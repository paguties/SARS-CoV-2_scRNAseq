#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
   Written by Pablo Gutierrez
   Grupo de Biotecnologia microbiana
   Universidad Nacional de Colombia sede Medellin
   '''

__author__ = "Pablo Gutierrez"
__copyright__ = "Copyright 2020, Grupo de Biotecnologia Microbiana Unalmed"
__credits__ = ["Pablo Gutierrez"]
__license__ = "MIT"
__version__ = "1.0.1"
__maintainer__ = "Pablo Gutierrez"
__email__ = "paguties@unal.edu.co"
__status__ = "Development"
__date__ = '2021/11/02'

import gzip
import os
import sys
import time

# ============================
# Input files
# ============================

# technical_reads = 'SRR12508049_1.fastq.gz'
# biological_reads = 'SRR12508049_2.fastq.gz'
# technical_reads = 'test_1.fq.gz'
# biological_reads = 'test_2.fq.gz'

technical_reads  = sys.argv[1]
biological_reads  = sys.argv[2]
# query = "test_1e3_2.fq.gz"

# ============================
# Output files and removal of old files
# ============================
print(f"\n==========================================")
print(f"Removing old files")
print(f"==========================================")
start_time = time.time()

filtered_fasta = "filtered_reads_" + biological_reads.split(".")[0] + ".fa"

if os.path.exists(filtered_fasta):
    os.system("rm " + filtered_fasta)
    print(f"Removed file: {filtered_fasta}")

if os.path.exists(filtered_fasta + ".gz"):
    os.system("rm " + filtered_fasta + ".gz")
    print(f"Removed file: {filtered_fasta}.gz")

blast_mapping = "blast_mapping_" + biological_reads.split(".")[0] + ".tbl.gz"

unambiguous_mappings = "unambiguous_mapping_" + biological_reads.split(".")[0] + ".tbl"

if os.path.exists(unambiguous_mappings):
    os.system("rm " + unambiguous_mappings)
    print(f"Removed file: {unambiguous_mappings}")

if os.path.exists(unambiguous_mappings + ".gz"):
    os.system("rm " + unambiguous_mappings + ".gz")
    print(f"Removed file: {unambiguous_mappings}.gz")

sorted_mappings = "sorted_" + biological_reads.split(".")[0] + ".tbl.gz"

if os.path.exists(sorted_mappings):
    os.system("rm " + sorted_mappings)
    print(f"Removed file: {sorted_mappings}")

curated_mappings = "curated_mappings_" + biological_reads.split(".")[0] + ".tbl"

if os.path.exists(curated_mappings):
    os.system("rm " + curated_mappings)
    print(f"Removed file: {curated_mappings}")

if os.path.exists(curated_mappings + ".gz"):
    os.system("rm " + curated_mappings + ".gz")
    print(f"Removed file: {curated_mappings}.gz")

print(f'Total time: {round(time.time() - start_time, 2)} seconds')

# ============================
# Creation of whitelist set
# ============================
print(f"\n==========================================")
print(f"Reading whitelist")
print(f"==========================================")
start_time = time.time()
whitelist_set = set()
with gzip.open('DB/3M-february-2018.txt.gz', 'rt') as file:
    for line in file:
        whitelist_set.add(line.split()[0])

print(f'Whitelist set created')
print(f'Total time: {round(time.time() - start_time, 2)} seconds')

# ============================
# Removal of bad quality reads
# ============================
print(f"\n==========================================")
print(f"Removal of bad quality reads")
print(f"==========================================")
start_time = time.time()

print(f'Technical reads: {technical_reads}')
print(f'Biological reads: {biological_reads}')
print(f'Generating {filtered_fasta} file')
with gzip.open(technical_reads, 'rt') as f1, gzip.open(biological_reads, 'rt') as f2:
    for line1 in f1:
        temp1 = line1.rstrip()
        BC = temp1[:16]
        UMI = temp1[16:]
        temp2 = next(f2).rstrip()
        if BC in whitelist_set and 'N' not in temp1:
            with open(filtered_fasta, 'a') as output_file:
                output_file.write(f">{BC}_{UMI}\n{temp2}\n")
print(f'Total time: {round(time.time() - start_time, 2)} seconds')
print(f'Compressing {filtered_fasta} file')
start_time = time.time()
os.system("gzip " + filtered_fasta)
print(f'Total time: {round(time.time() - start_time, 2)} seconds')

# ============================
# Mapping of reads
# ============================
print(f"\n==========================================")
print(f"Read mapping step")
print(f"==========================================")
start_time = time.time()
BLAST_command = ("magicblast -db DB/reference -no_unaligned  -reftype transcriptome "
                 "-outfmt tabular "
                 "-num_threads 6 " +
                 # "-score 95  "
                 "-score 85  "
                 "-query " + filtered_fasta + ".gz" +
                 " | cut -f 1,2,15 | grep 'plus$' | cut -f 1,2 | tr '(|)' '\t' "
                 " | cut -f 1,3 | uniq | gzip > " + blast_mapping)

print(f'Running magicblast on {filtered_fasta} file')
os.system(BLAST_command)
print(f'Read mapping step finished')
print(f'Total time: {round(time.time() - start_time, 2)} seconds')

# ============================
# Removal of ambiguous mappings
# ============================

print(f"\n==========================================")
print(f"Removal of ambiguous mappings")
print(f"==========================================")
start_time = time.time()

line1 = str()
line2 = str()
line3 = str()

with gzip.open(blast_mapping, 'rt') as file:
    line1 = next(file)
    line2 = next(file)
    for line in file:
        line3 = line
        if line2.split()[0] != line1.split()[0] and line2.split()[0] != line3.split()[0]:
            with open(unambiguous_mappings, 'a') as output_file:
                output_file.write(line2)
            # print(line2)
        line1 = line2[:]
        line2 = line3[:]

if line3.split()[0] != line2.split()[0]:
    # print(line3)
    with open(unambiguous_mappings, 'a') as output_file:
        output_file.write(line2)

print(f'Deleting: {blast_mapping}')
os.system("rm " + blast_mapping)
print(f'Sorting {unambiguous_mappings}')
sorting_command = os.system("sort -k 1 --radixsort --mmap " + unambiguous_mappings + " | gzip > " + sorted_mappings)
print(f'Saved {sorted_mappings} file')
os.system("rm " + unambiguous_mappings)
print(f'Removed {unambiguous_mappings}')
print(f'Total time: {round(time.time() - start_time, 2)} seconds')
# ============================
# Removal of chimeric reads
# ============================
print(f"\n==========================================")
print(f"Removal of chimeric reads")
print(f"==========================================")
start_time = time.time()

line1 = str()
line2 = str()
line3 = str()

with gzip.open(sorted_mappings, 'rt') as file:
    line1 = next(file)
    line2 = next(file)
    for line in file:
        line3 = line
        if line2.split()[0] != line1.split()[0] and line2.split()[0] != line3.split()[0]:
            with open(curated_mappings, 'a') as output_file:
                output_file.write(line2)
        line1 = line2[:]
        line2 = line3[:]

if line3.split()[0] != line2.split()[0]:
    with open(curated_mappings, 'a') as output_file:
        output_file.write(line2)

print(f'Saved {curated_mappings} file')
os.system("rm " + sorted_mappings)
print(f'Removed {sorted_mappings}')
print(f'compressing {curated_mappings}')
os.system("gzip " + curated_mappings)
print(f'Total time: {round(time.time() - start_time, 2)} seconds')

print('\nEnd of execution')

# ============================
# End of script
# ============================
