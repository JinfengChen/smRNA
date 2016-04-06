#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
from numpy import *
import re
import os
import argparse

def usage():
    test="name"
    message='''
python smRNApipe_mirDeep.py

INPUT: 
SAM file of smRNA mapped to reference (adaptor cliped and quality filtered, produced by SmallRNA.pl)

STEP:
1. reads collapse
2. clean rRNA/tRNA/snRNA/snoRNA
3. mirDeep2 to predict mirRNA
4. summary siRNA and prepare for downstream analysis

OUTPUT: 
A summary file of smRNA (mirRNA and siRNA)
For mirRNA:
tag_id	mirRNA_name	NA	chromosome:start-end	read_count	TPM(transcripts per million reads)
For siRNA: (This siRNA could be used to downstream analysis with TE)
tag_id	gene	UTR/exon/intron	chromosome:start-end	read_count	TPM(transcripts per million reads)
tag_id  RIRE2_I	LTR/Gypsy chromosome:start-end    read_count      TPM(transcripts per million reads)
For other smRNA:
tag_id	tRNA/rRNA/snRNA	NA	chromosome:start-end    read_count      TPM(transcripts per million reads)

    '''
    print message

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)


if __name__ == '__main__':
    main()

