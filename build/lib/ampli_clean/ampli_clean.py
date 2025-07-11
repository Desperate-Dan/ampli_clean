#!/usr/bin/env python
from Bio import SeqIO
import itertools
import os
import argparse
import pysam
import sys
from collections import defaultdict
import re
import gzip

#Minimap2 is a prerequisite in the env where this is run... could add a "which minimap2" check to see if it's installed!


def read_parser(input_files, min_len=False, max_len=False, log=False):
    #Deal with the read file input...currently need one read file for input to minimap2
    #This section either reads all input lines and filters them, or if no filtering is needed skips that and zcats the files together

    if min_len or max_len:
        read_bin = gzip.open("binned_reads.fastq.gz", "wt")
        counter = 1
        read_count = 0    
        for file in input_files:
            logger("Binning read file %s" % counter, log)
            read_file = gzip.open(file, "rt")
            for record in SeqIO.parse(read_file, "fastq"):
                #Filtering happens here, possible that this could be sped up with SeqIO.index...    
                if min_len and max_len:
                    if min_len < len(record) < max_len:
                        SeqIO.write(record, read_bin, "fastq")
                        read_count += 1
                elif min_len:
                    if min_len < len(record):
                        SeqIO.write(record, read_bin, "fastq")
                        read_count += 1
                elif max_len:
                    if len(record) < max_len:
                        SeqIO.write(record, read_bin, "fastq")
                        read_count += 1
            read_file.close()
            counter += 1
        logger("%s reads have been binned after filtering" % read_count, log)    
        read_bin.close()
        
        return read_bin

    else:
        #Simple but gross looking command to run zcat if filtering is not needed, probably a more succinct way of calling this
        os.system("zcat %s > ./binned_reads.fastq.gz" % str(input_files).replace(",","").lstrip("[").rstrip("]"))


def mini_mapper(output_name, input_ref, secondary=False, log=False):
    #Map, sort and index ready for cleaning
    #Extract the ref names
    with open(input_ref, "r") as refs:
        ref_names = []
        for line in refs:
            if re.match(">", line):
                if re.match("\s*", line):
                    ref_names.append(line.lstrip(">").rstrip("\n").split()[0])    
                else:
                    ref_names.append(line.lstrip(">").rstrip("\n"))
        logger("Found %s references: %s" % (len(ref_names),ref_names), log)
        if secondary:
            sec = "--secondary=no"
        else:
            sec = ""
    #Echo the minimap command
    logger("minimap2 -a %s -x map-ont -o %s.bam %s binned_reads.fastq.gz" % (sec, output_name, input_ref), log)
    os.system("minimap2 -a %s -x map-ont -o %s.bam %s binned_reads.fastq.gz" % (sec, output_name, input_ref))
    
    pysam.sort("-o", "%s.sorted.bam" % output_name, "%s.bam" % output_name, catch_stdout=False)
    pysam.index("%s.sorted.bam" % output_name, catch_stdout=False)

    return ref_names

def bed_file_reader(input_bed):
    #Read a bed file in format CHR, start_pos, end_pos, primer_id, etc and return an amplicon dictionary
    #Currently does not handle "alt" primers
    bed = open(input_bed)
    primer_list = []
    primer_pos_dict = {}
    for line in bed:
        #Take each line of the bed file, strip \n, split on tab, keep first 4 elements to deal with 4 or 6 element bed files
        primer_list.append(line.rstrip("\n").split("\t")[:4])
    for line in primer_list:
        bed_ref = line[0]
        if re.match(".*_LEFT", line[-1]):
            #This is a bit dodgy as it currently would need the LEFT primers to be first in the bed file... Maybe that is a safe assumption?
            #Could add a try: statement to deal with this?
            primer_pos_dict[line[3].rstrip("_LEFT")] = {}
            primer_pos_dict[line[3].rstrip("_LEFT")].update({'LEFT_START': int(line[1]),
                                                      'LEFT_END': int(line[2])})
        elif re.match(".*_RIGHT", line[-1]):
            primer_pos_dict[line[3].rstrip("_RIGHT")].update({'RIGHT_START': int(line[1]),
                                                      'RIGHT_END': int(line[2])})

    return primer_pos_dict, bed_ref


def ampli_clean(input_file,ref_names,output_name, bed_dict, wobble, all_vs_all=False, log=False):
    #This code iterates through the bam file and checks if each read starts and ends within a pair of primer positions.
    
    #Open the alignment file - expects BAM at the moment.
    aln_file = pysam.AlignmentFile(input_file,'rb')

    read_count_dict = map_stats(ref_names,aln_file, log)

    #This checks if the ref with the most mapping reads is greater than 0 and exits at this point - necessary for clean nextflow execution
    if read_count_dict[max(read_count_dict, key=read_count_dict.get)] == 0:
        aln_file.close()
        logger("%s has no mapped reads... Exiting..." % input_file, log)
        sys.exit(0)

    if not all_vs_all or len(bed_dict) == 1 :
        ref_names = []
        ref_names.append(max(read_count_dict, key=read_count_dict.get))
    
    for ref in ref_names:
        primer_position_dict = bed_dict[ref] #Should add an error check here for if the bed ref and the fasta ref names match - perhaps could ignore this step if there is only one bed file given?
        #For each ref mapped against set up some files
        logger("cleaning %s..." % ref, log)
        outfile = pysam.AlignmentFile("%s.%s.clean.bam" % (output_name, ref), "wb", template=aln_file)
        #Counter for cleaned reads per ref
        counter = 0
        #Fetch allows you to pull out reads matching to a single ref - could make a crude counter if we only want X number of refs returned
        for read in aln_file.fetch(ref):
            for primer in primer_position_dict:
                #This would be "strict" where a sequence fragment needs to start and end in a primer site. Could also add a more permissive mode?
                if (primer_position_dict[primer]["LEFT_START"] - wobble < read.reference_start < primer_position_dict[primer]["LEFT_END"] + wobble):
                    if (primer_position_dict[primer]["RIGHT_START"] - wobble < read.reference_end < primer_position_dict[primer]["RIGHT_END"] + wobble):
                        counter += 1
                        outfile.write(read)
        logger("%s has ~%s mapped reads after cleaning..." % (ref, counter), log)

        #Need to remove empty outfiles after cleaning
        #It seems that longhot in fieldbioinformatics will fail if passed an empty file or a file with exactly one mapped read so filtering based on that
        if counter <= 1:
            logger("One or less reads left after cleaning...removing output...")
            os.remove("%s.%s.clean.bam" % (output_name, ref))
            sys.exit(0)

        outfile.close()
        

    aln_file.close()
    
    return ref_names


def sam_sort_index(output_name,ref_name):
    #Sort and index your bam file. Works on the assumption that the outfile is the correct name for your bam, may need to change.
    clean_bam_sorted = pysam.sort("-o", "%s.%s.clean.sorted.bam" % (output_name, ref_name), "%s.%s.clean.bam" % (output_name, ref_name))
    clean_bam_index = pysam.index("%s.%s.clean.sorted.bam" % (output_name, ref_name))    

    return clean_bam_sorted, clean_bam_index


def bam_to_fq(output_name,ref_name):
    #Return your cleaned bam to fq for sticking into fieldbioinf
    fq = pysam.fastq("-0", "%s.%s.fastq.gz" % (output_name, ref_name), "%s.%s.clean.bam" % (output_name, ref_name))

    return fq


def map_stats(ref_names, aln_file, log):
    read_count_dict = {}
    for ref in ref_names:
        read_count_dict[ref] = aln_file.count(ref)
        logger("%s has ~%s mapped reads..." % (ref, read_count_dict[ref]), log)
    
    return read_count_dict


def logger(msg,log=False):
    if log:
        with open("log.txt", "a") as log_file:
            log_file.write(msg + "\n")
    else:
        os.system("echo %s" % msg)


def skip_clean(ref_names, input_bam, output_name, log=False, out_fastq=False, read_cutoff=50):
    with pysam.AlignmentFile(input_bam, 'rb') as aln_file:
        #Get a dict of counts of reads for each ref using map_stats
        read_count_dict = map_stats(ref_names, aln_file, log)
        #Parse that and kick out any that are below the cut off
        above_cutoff_dict = {}
        for ref in read_count_dict:
            if read_count_dict[ref] >= read_cutoff:
                above_cutoff_dict[ref] = read_count_dict[ref]
        if not above_cutoff_dict:
            logger("No refs have more than %s reads mapping to them...exiting..." % read_cutoff)
            sys.exit(0)
        #Pull out the reads that correspond to each remaining ref
        for ref in above_cutoff_dict:
            with pysam.AlignmentFile("%s.%s.bam" % (output_name, ref), "wb", template=aln_file) as outfile:
                for read in aln_file.fetch(ref):
                    outfile.write(read)
            if out_fastq:
                pysam.fastq("-0", "%s.%s.fastq.gz" % (output_name, ref), "%s.%s.bam" % (output_name, ref))
        sys.exit(0)


def runner(args):
    #Gathers the fastq.gz files and filters if required
    read_parser(args.input_reads,args.min_len,args.max_len,args.log)
    #Does the mapping and gets the reference names from the input ref file
    ref_names = mini_mapper(args.output_name,args.input_ref,args.sec,args.log)
    if args.map_only:
        skip_clean(ref_names, "%s.sorted.bam" % args.output_name, args.output_name, args.log, args.out_fastq, args.read_cutoff)
    #Parses the input bed files to get primer positions per amplicon and adds them to a dictionary
    bed_pos_dict = {}
    for bed_file in args.input_bed:
        primer_position_dict, bed_ref = bed_file_reader(bed_file)
        bed_pos_dict[bed_ref] = primer_position_dict
    #Runs the mapping stats then cleans the bam file using the bed_pos_dict
    ref_names = ampli_clean("%s.sorted.bam" % args.output_name,ref_names,args.output_name,bed_pos_dict,args.wobble,args.all_vs_all, args.log)
    #Optional: Sorts and indexes the cleaned bam file
    if args.out_sort:
        for ref in ref_names:
            sam_sort_index(args.output_name,ref)
    #Optional: Pulls the fastqs out of the cleaned bam; necessary for input into further pipelines
    if args.out_fastq:
        for ref in ref_names:
            bam_to_fq(args.output_name,ref)


def main():
    parser = argparse.ArgumentParser()
    '''
    Parse the command line arguments.
    '''
    parser = argparse.ArgumentParser(description='Creates a "clean" bam file containing only reads that start and end near primer sites.')

    parser.add_argument('-w', '--wobble', dest='wobble', type=int, default=10,
                            help='Number of bases around the primer binding sites that reads can start or end and still be considered amplicon spanning. Default = 10')
    parser.add_argument('-o', '--output-name', dest='output_name', default="clean",
                            help='Prefix for the output. Default = "clean"')
    parser.add_argument('-s', dest = 'out_sort', action='store_true',
                            help='Output sorted and indexed bam')
    parser.add_argument('--fastq', dest = 'out_fastq', action='store_true',
                            help='Output cleaned fastq file')
    parser.add_argument('--all', dest= 'all_vs_all', action='store_true',
                            help='Cleans each references mapped reads using the appropriate input bed file rather than ust the one with the most mapped reads')
    parser.add_argument('--secondary', dest = 'sec', action='store_false',
                            help='Allow minimpa2 to output secondary alignments. Default = False')
    parser.add_argument('--min', dest = 'min_len', type=int,
                            help='Filter reads when binning by minimum read length')
    parser.add_argument('--max', dest='max_len', type=int,
                            help='Filter reads when binning by maximum read length')
    parser.add_argument('--log', dest='log', action='store_true',
                            help='Write messages to log file')
    parser.add_argument('--bam', dest = 'input_bam',
                            help='Path to the BAM file you want to clean (currently non-functional)')
    parser.add_argument('-n', '--ref-name', dest='ref_name',
                            help='Name of ref the bam files were aligned to (currently non-functional)')
    parser.add_argument('--map-only', dest='map_only', action='store_true',
                            help='Skip the bam cleaning with bedfile step, only map against the references provided. The output flags still work.')
    parser.add_argument('--map-cutoff', dest='read_cutoff', type=int, default=50,
                            help='Cutoff for number of mapped reads when using --map-only. Default = 50')

    
    required_group = parser.add_argument_group('required arguments')
    required_group.add_argument('-r', '--refs', dest = 'input_ref',
                            help='Path to input ref fasta')
    required_group.add_argument('-f', dest = 'input_reads', nargs='+',
                            help='Path to input fastq, currently expects them to be gzipped')
    required_group.add_argument('-b', '--bed', dest = 'input_bed', nargs='+',
                            help='Path to the bed file you want to get positions from')


    args = parser.parse_args()
    runner(args)    
    
          
if __name__ == "__main__":
    main()

