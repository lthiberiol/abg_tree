#!/usr/bin/env python
#coding: utf-8

import os, re
from sys import argv, exit

fasta_file = argv[1]

if not os.path.isfile(fasta_file):
    exit('%s does not exist' %fasta_file)

file_name_no_ext = os.path.splitext(fasta_file)[0]

os.system('mafft --reorder --auto --quiet %s > %s.aln' %(fasta_file, file_name_no_ext))
os.system('hmmbuild %s.hmm %s.aln' %(file_name_no_ext, file_name_no_ext))