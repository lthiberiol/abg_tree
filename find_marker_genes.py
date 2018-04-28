#!/usr/bin/env python
#encoding: utf-8

############################################################
#                                                          #
# Script to assess taxonomy available in Silva and define  #
#     genomes to sampled...                                #
#                                                          #
#                                       L. Thibério Rangel #
#                                     lthiberiol@gmail.com #
#                                                          #
############################################################

import sys

if __name__ != '__main__':
    sys.exit()

import os
import itertools
from popen2 import popen2


class cd:
    """
    Context manager for changing the current working directory
    """
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)

def qsub_hmmsearch(xargs_input_file, num_nodes=1, num_threads=1, memmory_size='1GB', time_limit='170:00:00', working_dir='~', genome_dir='.', profile_dir='.'):
    # Open a pipe to the qsub command.
    qsub_out, qsub_in = popen2('sbatch')

    # Customize your options here
    job_string = """#!/bin/bash
#SBATCH -p sched_mit_g4nier               #greg's partition on the cluster
#SBATCH -N {num_nodes}                    #number of nodes
#SBATCH -c {num_threads}                  #number of cores
#SBATCH --mem={memmory_size}              #max amount of memory
#SBATCH --time={time_limit}               #wall time
#SBATCH -J hmmsearch                      #job name
#SBATCH --output=hmmsearch.out            #output name
#SBATCH --error=hmmsearch.err             #error file name
#SBATCH --mail-user=lthiberiol@gmail.com  #if you want emails upon start/finish

module load engaging/hmmer/3.1b2
cd {working_dir}

cat {input_file} | xargs -n 2 -P {simultaneous_processes} sh -c 'hmmsearch --acc --noali --cpu 1 -o $1_-_$2.hmm_out {profile_dir}/$1 {genome_dir}/$2' sh
""".format(input_file=xargs_input_file, num_nodes=num_nodes, num_threads=num_threads, memmory_size=memmory_size,
           time_limit=time_limit, working_dir=working_dir, simultaneous_processes=num_threads*num_nodes,
           genome_dir=genome_dir, profile_dir=profile_dir)

    # Send job_string to qsub
    qsub_in.write(job_string)
    qsub_in.close()

    # Print your job and the system response to the screen as it's submitted
    return qsub_out.read()

profile_search_folder = '/home/thiberio/abg_tree/profile_search'
profile_folder        = '/home/thiberio/abg_tree/hmm_models'
genome_folder         = '/home/thiberio/abg_tree/faas'

if not os.path.isdir(profile_folder) or not os.path.isdir(genome_folder):
    exit('Source directories do not exists...')

print 'Running maker HMM against genomes!'
if not os.path.isdir(profile_search_folder):
    os.mkdir(profile_search_folder)
else:
    os.system('rm -rf %s/*' %profile_search_folder)

with cd(profile_search_folder):
    out = open('tmp_xargs_list', 'wb')
    for line in itertools.product(os.listdir(profile_folder), os.listdir(genome_folder)):
        out.write('%s\n' %' '.join(line))
    out.close()

print qsub_hmmsearch('tmp_xargs_list', num_nodes=2, num_threads=20, memmory_size='10GB', working_dir=profile_search_folder, genome_dir=genome_folder, profile_dir=profile_folder)
