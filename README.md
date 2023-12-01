# ampli_clean
A tool to clean up bam files from amplicon based sequencing methods so only sequences that start and/or end in primer site pairs will remain. Will hopefully be useful to deal with off mapping when using multiple different amplicon schemes simultaneously.

## Installation
`git clone https://github.com/Desperate-Dan/ampli_clean.git`

`cd ampli_clean`

`pip install .`

To check installation:

`ampli_clean -h`

## Usage
ampli_clean accepts either a single fastq file or multiple by globbing. If multiple files are provided it will then bin these together for cleaning. Multiple reference files can be given to ampli_clean and they will be competitively mapped to. Multiple bed files can be provided, but the first column of these bed files must match the names of their corresponding ref for mapping. Default behaviour chooses the reference with most mapped reads and then cleans the reads using the matching provided bed file. Default output returns a cleaned bam file but the fastq reads from this file can be output, ready to be provided to your consensus generation method of choice. This tool assumes your reads of interest will be amplicon spanning (ie. start and end near a pair of primer binding sites) and removes reads that do not start or end within the left and right primer sites of an amplicon in the bed file provided. The `--wobble=n` setting allows reads that start or end `n` bases up or downstream of the primer binding sites to be included. 

It is important to note that if a primer for a different scheme happens to bind at a site only slighlty offset from a "true" binding site from the scheme of interest, it is possible that those reads will not be correctly filtered out. As such we recommend using consensus generating tools that check both amplicons at overlapping regions for the presence of a SNP, such as the `--strict` flag in the `fieldbioinformatics` pipeline.


```
usage: ampli_clean [-h] [-w WOBBLE] [-o OUTPUT_NAME] [-s] [--fastq] [--all] [--secondary]
                   [--min MIN_LEN] [--max MAX_LEN] [--log] [-r INPUT_REF]
                   [-f INPUT_READS [INPUT_READS ...]] [-b INPUT_BED [INPUT_BED ...]]

Creates a "clean" bam file containing only reads that start and end near primer sites.

optional arguments:
  -h, --help            show this help message and exit
  -w WOBBLE, --wobble WOBBLE
                        Number of bases around the primer binding sites that reads can start or end and still be considered amplicon spanning. Default = 10 
  -o OUTPUT_NAME, --output-name OUTPUT_NAME
                        Prefix for the output. Default = "clean"
  -s                    Output sorted and indexed bam
  --fastq               Output cleaned fastq file
  --all                 Cleans each references mapped reads using the appropriate input bed file rather
                        than just the one with the most mapped reads
  --secondary           Allow minimpa2 to output secondary alignments. Default = False
  --min MIN_LEN         Filter reads when binning by minimum read length
  --max MAX_LEN         Filter reads when binning by maximum read length
  --log                 Write messages to log file

required arguments:
  -r INPUT_REF, --refs INPUT_REF
                        Path to input ref fasta
  -f INPUT_READS [INPUT_READS ...]
                        Path to input fastq, currently expects them to be gzipped
  -b INPUT_BED [INPUT_BED ...], --bed INPUT_BED [INPUT_BED ...]
                        Path to the bed file you want to get positions from
```
