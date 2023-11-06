#!/usr/bin/env python
from Bio import SeqIO
import itertools
import os
import argparse
import pysam
import sys

#Main body of the function

def runner(args):
    #Open the alignment file - expects BAM at the moment
    aln_file = pysam.AlignmentFile("barcode21.sorted.bam",'rb')
    
    #Make a dictionary of length = to the ref sequence
    coverage_dict = {}
    for i in range(0 ,aln_file.get_reference_length(args.ref_name)):
        coverage_dict[i] = 0
    
    #Populate that dictionary with the read depth at each position
    for pileup_column in aln_file.pileup(args.ref_name, truncate=False, min_base_quality=0):
        coverage_dict[pileup_column.pos] = pileup_column.n
    
    #Create a list of lists containing the runs of positions below the "depth" value
    mask_pos_list = []
    position_list = []
    for pos in coverage_dict:
        if coverage_dict[pos] < int(1):
            position_list.append(pos)
        else:
            if position_list:
                mask_pos_list.append(position_list)
            position_list = []    
    
    #Write a "bcftools consensus" friendly file of regions to mask
    mask_file = open(args.output_name + '.tsv', 'w')
    for mask_region in mask_pos_list:
<<<<<<< HEAD:build/lib/RSV_clean/rsv_clean.py
=======
        #bcftools expects your mask file to be one based so need to add one to all co-ordinates
>>>>>>> 2962088b9f7b8a81a876aec5cdd9be852c77ec24:maskara/maskara.py
        mask_file.write("%s\t%d\t%d\n" % (args.ref_name, mask_region[0] + 1, mask_region[-1] + 1 ))
    mask_file.close()

    aln_file.close()
    return coverage_dict, mask_file




def main():
    parser = argparse.ArgumentParser()
    '''
    Parse the command line arguments.
    '''
    parser = argparse.ArgumentParser(description='Creates a coverage mask to apply to your lovely consensus fasta.')

    
    optional_group = parser.add_argument_group('Optional')
    optional_group.add_argument('-d', '--depth', dest='depth', default="20",
                            help='If coverage is below this it will be masked. Default = 20')
    optional_group.add_argument('-r', '--ref-name', dest='ref_name', default="MN908947.3",
                            help='Name of ref the bam files were aligned to. Default = "MN908947.3"')
    optional_group.add_argument('-o', '--output-name', dest='output_name', default="depth_mask",
                            help='Prefix for the output. Default = "depth_mask"')
    
    parser.add_argument('input_file',
                            help='Path to the BAM file you want to create a mask for')


    args = parser.parse_args()
    runner(args)    
    
    
    
      
if __name__ == "__main__":
    main()

