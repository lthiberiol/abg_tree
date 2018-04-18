#!/usr/bin/env python
#coding: utf-8

import pandas as pd
import ete3
import re
import os
os.chdir('/work/abg_tree')

genome_lineages = pd.read_table( '/work/lbi_backup/fournierLab/lineages_from_genbank_summary.tab', sep='\t', index_col=0, dtype=str )
genome_lineages.index = genome_lineages.index.astype( str )
genome_lineages['taxid'] = genome_lineages.index

genbank_summary                     = pd.read_table( '/work/lbi_backup/assembly_summary_genbank.txt', dtype={'taxid':str, 'infraspecific_name':str} )
genbank_summary['refseq_category']  = genbank_summary['refseq_category'].str.lower()
genbank_summary['assembly_level']   = genbank_summary['assembly_level'].str.lower()
genbank_summary['genome_rep']       = genbank_summary['genome_rep'].str.lower()
genbank_summary = genbank_summary.merge( genome_lineages[['taxid', 'class', 'phylum', 'order']], how='right', on='taxid' )
genbank_summary.set_index( 'assembly_accession', inplace=True )

refseq_summary                     = pd.read_table( '/work/lbi_backup/assembly_summary_refseq.txt', dtype={'taxid':str, 'infraspecific_name':str} )
refseq_summary['refseq_category']  = refseq_summary['refseq_category'].str.lower()
refseq_summary['assembly_level']   = refseq_summary['assembly_level'].str.lower()
refseq_summary['genome_rep']       = refseq_summary['genome_rep'].str.lower()
refseq_summary = refseq_summary.merge( genome_lineages[['taxid', 'class', 'phylum', 'order']], how='left', on='taxid' )
refseq_summary.set_index( 'assembly_accession', inplace=True )

assembly_summary = refseq_summary.append(genbank_summary)

for tree_file in os.listdir('phylogeny/'):
    if not tree_file.endswith('aln.treefile'):
        continue

    tree = ete3.Tree('phylogeny/%s' %tree_file)
    taxa = {'class':set(), 'phylum':set(), 'order':set()}
    out  = open('phylogeny/%s.figTree' %tree_file, 'wb')
    out.write("#NEXUS\nbegin taxa;\n\tdimensions ntax=%i;\n\ttaxlabels\n" %len(tree))
    for leaf in tree.get_leaves():
        if leaf.name.startswith('GCF_') or leaf.name.startswith('GCA_'):
            genome, gene = re.search('^(GC[FA]_\d+\.\d+)_(\S+)$', leaf.name).groups()
        else:
            genome, gene = re.search('^(\S+)_(\S+)$', leaf.name).groups()

        if genome not in assembly_summary.index:
            out.write('\t%s\n' %(leaf.name))
            continue

        out.write('\t%s ' %(leaf.name))
        comment = []
        for rank in ['organism_name', 'class', 'phylum', 'order']:
            comment.append('tax_%s="%s"' %(rank, assembly_summary.loc[genome, rank]))

        out.write('[&%s]\n' %' '.join(comment))
    out.write(';\nend;\n')

    out.write('begin trees;\n\ttree tree_1 = [&R] %s\nend;' %tree.write())
    out.close()
