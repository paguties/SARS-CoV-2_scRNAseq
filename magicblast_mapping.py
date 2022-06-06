#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
   Written by Pablo Gutierrez
   '''

__author__ = "Pablo Gutierrez"
__credits__ = ["Pablo Gutierrez"]
__license__ = "MIT"
__version__ = "1.0.1"
__maintainer__ = "Pablo Gutierrez"
__email__ = "paguties@unal.edu.co"
__status__ = "Development"
__date__ = '2021/11/02'

# run as: python magicblast_mapping.py xxx.fq.gz

# ============================
# Imports
# ============================
import os
import sys
import gzip
import numpy as np


# ============================
# MAGICBLAST parameters
# ============================

query = sys.argv[1]
# query = "test_1e3_2.fq.gz"
output = "magicblast_" + query.split(".")[0] + ".tbl.gz"
unambiguous_output = "unambiguous_mappings_" + query.split(".")[0] + ".tbl"

if os.path.exists(unambiguous_output):
    os.system("rm " + unambiguous_output )
    print(f"Removed file: {unambiguous_output}")
if os.path.exists(unambiguous_output+".gz"):
    os.system("rm " + unambiguous_output+".gz" )
    print(f"Removed file: {unambiguous_output}.gz")


BLAST_command = ("magicblast -db DB/reference -no_unaligned  -reftype transcriptome "
                 "-outfmt tabular -infmt fastq "
                 "-num_threads 6 " +
                 "-score 95  "
                 "-query " + query +
                 " | cut -f 1,2,15 | grep 'plus$' | cut -f 1,2 | tr '(|)' '\t' "
                 " | cut -f 1,3 | uniq | gzip > " + output)


os.system(BLAST_command)

# ============================
# Removal of ambiguous mappings
# ============================

line1 = str()
line2 = str()
line3 = str()


with gzip.open(output, 'rt') as file:
    line1 = next(file)
    line2 = next(file)
    for line in file:
        line3 = line
        if line2.split()[0] != line1.split()[0] and line2.split()[0] != line3.split()[0]:
            with open(unambiguous_output, 'a') as output_file:
                output_file.write(line2)
            # print(line2)
        line1 = line2[:]
        line2 = line3[:]

if line3.split()[0] != line2.split()[0]:
    # print(line3)
    with open(unambiguous_output, 'a') as output_file:
        output_file.write(line2)

# os.system("rm " + output )
os.system("gzip " + unambiguous_output )


# ============================
# End of script
# ============================
