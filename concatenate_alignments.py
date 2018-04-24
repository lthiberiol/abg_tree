#!/usr/bin/env python
#coding: utf-8

from Bio import AlignIO, SeqIO, Align, Alphabet
import pandas as pd
import os, re, sys
from copy import deepcopy

aln_alphabet = Alphabet.Gapped(Alphabet.IUPAC.ambiguous_dna)

aln_folder    = '/work/abg_tree/concatenated_trees/3rd_try/alignments'
output_folder = '/work/abg_tree/concatenated_trees/3rd_try'
genomes    = {}
for aln in os.listdir(aln_folder):
    alignment    = AlignIO.read('%s/%s' %(aln_folder, aln), 'fasta')
    genomes[aln] = set()
    for entry in alignment:
        if re.match('GC[AF]_', entry.name):
            genome, gene = entry.name.split('|')
        else:
            genome, gene = entry.name.split('_')

        if genome in genomes[aln]:
            sys.exit('\t**Error, duplicated genome in %s: %s' %(aln, genome))

        genomes[aln].add(genome)

genome_union = set.union(*genomes.values())

missing_genes = {} # just to keep track of the number of missing marker genes in each genome
concatenation = {}
for genome in genome_union:
    missing_genes[genome]             = 0
    concatenation[genome]             = Align.SeqRecord( Align.Seq('', aln_alphabet) )
    concatenation[genome].name        = genome
    concatenation[genome].id          = genome
    concatenation[genome].description = genome

#
# fill the handles with the marker sequences from each genome
total_genes      = 0.0 # keep track of the number of genes added to the concatenation
current_position = 1
partitions       = open('%s/concatenated_partitions' %output_folder, 'wb')
for aln in os.listdir(aln_folder):

    tmp_aln      = AlignIO.read( '%s/%s' %(aln_folder, aln), 'fasta' )
    aln_length   = tmp_aln.get_alignment_length() # get the expected size of the alignment so you can compare if all have the same size
    total_genes += aln_length

    for entry in tmp_aln:
        # if this alignment has a different size from the rest, something is reaaaaaly wrong!
        if len(entry) != aln_length:
            sys.exit('\t**Error, block "%s" has a different length than the rest of the MSA: %s' %(entry.name, aln))

        if re.match('GC[AF]_', entry.name):
            genome, gene = entry.name.split('|')
        else:
            genome, gene = entry.name.split('_')

        concatenation[genome] += deepcopy(entry.seq)

    partitions.write('LG, %s = %i-%i\n' %(aln.replace('.fasta.aln', ''), current_position, current_position+aln_length-1) )
    current_position += aln_length

    #
    # add gaps for those genomes missing this gene (same size as the expected alignment)
    for genome in genome_union.difference(genomes[aln]):
        concatenation[genome] += Align.Seq( '-' * aln_length, aln_alphabet )
        missing_genes[genome] += aln_length
partitions.close()

#
# remove genomes missing more than 20% of the marker genes
for genome, num_missing_genes in missing_genes.items():
    if num_missing_genes/total_genes > 0.1:
        print '\t\t**%s: excluded from analysis for missing %.2f from concatenated alignment!' %(genome, (num_missing_genes/total_genes)*100)
        concatenation.pop( genome )

AlignIO.write( Align.MultipleSeqAlignment( concatenation.values() ), '%s/concatenated_alignment.aln' %output_folder, 'fasta' )