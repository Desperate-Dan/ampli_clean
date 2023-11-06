#!python
from Bio import SeqIO
import itertools
import os
import argparse
import pysam
import sys
from collections import defaultdict
import re


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
        if re.match(".*_LEFT", line[-1]):
            #This is a bit dodgy as it currently would need the LEFT primers to be first in the bed file... Maybe that is a safe assumption?
            primer_pos_dict[line[3].rstrip("_LEFT")] = {}
            primer_pos_dict[line[3].rstrip("_LEFT")].update({'LEFT_START': int(line[1]),
                                                      'LEFT_END': int(line[2])})
        elif re.match(".*_RIGHT", line[-1]):
            primer_pos_dict[line[3].rstrip("_RIGHT")].update({'RIGHT_START': int(line[1]),
                                                      'RIGHT_END': int(line[2])})

    return (primer_pos_dict)


def rsv_clean_main(input_file,ref_name,output_name, primer_position_dict, wobble):
    #
    
    #Open the alignment file - expects BAM at the moment
    aln_file = pysam.AlignmentFile(input_file,'rb')
    outfile = pysam.AlignmentFile(output_name + ".bam", "wb", template=aln_file)

    #New code
    #start_pos_list = []
    for read in aln_file.fetch(ref_name):
        for primer in primer_position_dict:
            if (primer_position_dict[primer]["LEFT_START"] - wobble < read.reference_start < primer_position_dict[primer]["LEFT_END"] + wobble):
                if (primer_position_dict[primer]["RIGHT_START"] - wobble < read.reference_end < primer_position_dict[primer]["RIGHT_END"] + wobble):
                    outfile.write(read)
    #                start_pos_list.append(read.reference_start)
    #print(start_pos_list)

    aln_file.close()
    outfile.close()
    return outfile


def sam_sort_index(output_name):
    #Sort and index your bam file. Works on the assumption that the outfile is the correct name for your bam, may need to change.

    clean_bam_sorted = pysam.sort("-o", "%s.sorted.bam" % output_name, "%s.bam" % output_name)
    clean_bam_index = pysam.index("%s.sorted.bam" % output_name)    
    
    return clean_bam_sorted, clean_bam_index


def bam_to_fq(output_name):
    #Return your cleaned bam to fq for sticking into fieldbioinf

    fq = pysam.fastq("-0", "%s.fastq.gz" % output_name, "%s.sorted.bam" % output_name)

    return fq


def runner(args):
    primer_position_dict = bed_file_reader(args.input_bed)
    rsv_clean_main(args.input_file,args.ref_name,args.output_name,primer_position_dict,args.wobble)
    sam_sort_index(args.output_name)
    bam_to_fq(args.output_name)


def main():
    parser = argparse.ArgumentParser()
    '''
    Parse the command line arguments.
    '''
    parser = argparse.ArgumentParser(description='Creates a "clean" bam file containing only reads that start and end near primer sites.')

    parser.add_argument('-w', '--wobble', dest='wobble', default=10,
                            help='If coverage is below this it will be masked')
    parser.add_argument('-r', '--ref-name', dest='ref_name', default="RSVA",
                            help='Name of ref the bam files were aligned to. Default = "RSVA"')
    parser.add_argument('-o', '--output-name', dest='output_name', default="clean",
                            help='Prefix for the output. Default = "clean"')
    
    required_group = parser.add_argument_group('required arguments')
    required_group.add_argument('-i', '--input', dest = 'input_file',
                            help='Path to the BAM file you want to clean')
    required_group.add_argument('-b', '--bed', dest='input_bed',
                            help='Path to the bed file you want to get positions from')


    args = parser.parse_args()
    runner(args)    
    
          
if __name__ == "__main__":
    main()

