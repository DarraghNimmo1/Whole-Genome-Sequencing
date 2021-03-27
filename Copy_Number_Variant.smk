#############################################################################
##Copy number variant analyses using GATK in the absence of a matched normal.
##Darragh Nimmo, Trinity College Dublin.
##March 2021
##############################################################################

##############################################################################
##Imports and variables.
##############################################################################

import os

import functools

configfile: 'config.yaml'

GENOME = config['genome']

IN_BAM = config['in_bam']

directory_function = functools.partial(os.path.join, config['results'])
INTERVALS_DIR = directory_function('Intervals')
COUNTS_DIR  = directory_function('Counts')

################################################################################
##Rules
################################################################################

rule get_intervals:
  input:
      fasta = GENOME
  output:
      list = os.path.join(INTERVALS_DIR,  '{sample}.interval_list')
  shell:
      """
      gatk PreprocessIntervals -R {input.fasta} --bin-length 1000 --padding 0 -O {output}
      """
      
rule read_counts_sample:
  input:
      bam = os.path.join(IN_BAM, '{sample}.bam')
      intervals = rules.get_intervals.output.list
  output:
      counts = os.path.join(COUNTS_DIR, {sample}.counts.hdf5)
  shell:
      """
      gatk CollectReadCounts -I {input.bam} -L {input.intervals} --interval-merging-rule OVERLAPPING_ONLY -O {output.counts}
      """
   
rule   
    
    
