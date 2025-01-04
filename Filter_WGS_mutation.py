
#############################################################################
## A script to perform the filtering of WGS variants for downstream analysis
##############################################################################

#Remove qual score variants
#Remove variants in hartwig germline PON
#gnomad filter file needs to be fixed up


#packages

import os
import pybedtools as pybed
import pandas as pd
import re
import numpy as np
import shutil
import glob
import csv
import pysam

pd.set_option('display.max_columns', None)

##########################################################################################################################################
# Hartwig filters
##########################################################################################################################################

######################################
#Filter out radiation therapy samples
######################################

pd.set_option('display.max_columns', None)
Metadata_file = pd.read_table("/home/darragh/Hartwig/metadata.tsv",header=0)
Metadata_file = Metadata_file[["sampleId", "primaryTumorLocation", "hasRadiotherapyPreTreatment"]]
Metadata_file.columns = ["patient", "primaryTumorLocation", "hasRadiotherapyPreTreatment" ]
Metadata_file = Metadata_file.loc[Metadata_file['primaryTumorLocation'] == "Breast"]
Metadata_file = Metadata_file.loc[Metadata_file['hasRadiotherapyPreTreatment'] != "Yes"]

for i in Metadata_file['patient']:
  shutil.move( os.path.join("/home/darragh/Hartwig/Data/breast", i+".purple.somatic.vcf"), os.path.join('/home/darragh/Hartwig/Data/breast/no_radiation', i+".purple.somatic.vcf"))


##############################
#Filter out duplicates
##############################
#Done manually based on visual inspection of file names


###############################
#Convert vcfs to a bed format that can be parsed - also filter out any variants not pass 

#ls *.purple.somatic.vcf |while read line ; do egrep -v "^#" $line | awk '{ if($7 == "PASS") print "chr"$1"\t"$2-1"\t"$2"\t"$4"\t"$5"\t"$7}' > $line.bed ; done

###############################
#Filter samples with >20% of mutations in highly repetitive regions
###############################


#Filter out the repeat file to only keep repetitions of 1-5bp with >5 repetitions
repeat_file = pd.read_table("/home/darragh/perf/PERF/first_run/repeating_dna.hg19.txt",header=None)
repeat_file.columns = ["chr", "start", "end", "repeating_nucs", "repeat_length", "strand", "n_motif_repeats", "actual_repeat"]
#repeat_file = repeat_file[["chr", "start", "end","repeat_length","n_motif_repeats" ]]
repeat_file['repeat_size'] = repeat_file['repeat_length']/repeat_file['n_motif_repeats']
repeat_file = repeat_file.loc[repeat_file['n_motif_repeats'] > 5]
repeat_file = repeat_file.loc[repeat_file['repeat_size'] < 5]
repeat_file = repeat_file[["chr", "start", "end", "repeating_nucs","repeat_length","n_motif_repeats"]]
repeat_file.to_csv("/home/darragh/Hartwig/Repetitive_dna_filter_vcf.txt", sep='\t', index=False, header=None)
repeat_file = pybed.BedTool.from_dataframe(repeat_file)
#ls *.vcf |cut -f 2 -d '/' |while read line ; do egrep -v "^#" "$line" | grep "PASS" | awk '{print "chr"$1"\t"$2"\t"$2"\t"$4"\t"$5}'  > "$line".not_pon.bed ; done
vcf_files = glob.glob('/home/darragh/Hartwig/Data/breast/no_radiation/*.somatic.vcf.bed')

for file in vcf_files:
  mutation_file = pybed.BedTool(file)
  mutation_file_replicative = mutation_file.intersect(repeat_file, wa=True, u=True)
  if len(mutation_file_replicative)/len(mutation_file) > 0.2:
    print(os.path.basename(file))

#When all of the filenames are printed to the screen, manually move them into a new filtered folder


##############################
#Filter out Kataegis variants
##############################

import os
import pybedtools as pybed
import pandas as pd
import re
import numpy as np
import shutil
import glob
import csv
import pysam

vcf_files = glob.glob('/home/darragh/Hartwig/Data/breast/no_radiation/*.somatic.vcf.bed')

for file in vcf_files:
  vcf = pybed.BedTool(file)
  vcf = vcf.cluster(d=1000)
  vcf_df = pd.read_table(vcf.fn,header=None)
  #vcf_df = vcf_df.iloc[:, [0,1,3,4,6,11]]
  vcf_df.columns = ["chr", "start", "end", "ref", "alt", "annotation", "cluster"]
  #vcf_df['chr'] = 'chr'+str(vcf_df['chr'])
  vcf_df =  vcf_df.groupby('cluster').filter(lambda x : len(x)<6)
  output_file = os.path.join( "/home/darragh/Hartwig/Data/breast/no_radiation", os.path.basename(file)).split(".")[0] + ".not_pon.not.kataegis.bed"
  vcf_df.to_csv(output_file, sep='\t', index=False, header=False)


####################################################
#Create mutation matrix with all remaining mutations
#####################################################

import os
import pybedtools as pybed
import pandas as pd
import re
import numpy as np
import shutil
import glob
import csv
import pysam

def parse_unfiltered_vcfs(file):
  mutation_file=set()
  for line in open(file).read().splitlines():
    row=line.split('\t')
    Chr= row[0]
    Start=int(row[1])
    End=int(row[2])
    Ref = row[3]
    Alt = row[4]
    Start=str(Start)
    End=str(End)
    Name = os.path.basename(file)
    mutation_file.add('\t'.join([Chr,Start,End,Ref,Alt,Name]))
  return mutation_file

def run_hartwig():
  if not os.path.exists('/home/darragh/Hartwig/Data/breast/no_radiation/HARTWIG_merged_mutation_file_july_5'):
    Hartwig_file_list = glob.glob('/home/darragh/Hartwig/Data/breast/no_radiation/*.not_pon.not.kataegis.bed')
    with open("/home/darragh/Hartwig/Data/breast/no_radiation/HARTWIG_merged_mutation_file_july_5","w") as f:
        for file in Hartwig_file_list:
            try:
                mutation_file = parse_unfiltered_vcfs(file)
            except OSError:
                pass
            wr = csv.writer(f,delimiter="\n")
            wr.writerow(mutation_file)

    return

run_hartwig()

##############################
# Filter out variants in 1000G
##############################
#  
#egrep "^#" ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf >1000g_header
#egrep -v "^#" ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf | awk '{print "chr"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}' > 1000G_no_header
# awk '{ print $}
#cat 1000g_header 1000G_no_header > 1000G_hg19_chr_fixed.vcf
#bedtools subtract -A -a /home/darragh/Hartwig/Data/breast/no_radiation/HARTWIG_merged_mutation_file_july_5 -b 1000G_hg19_chr_fixed.bed > /home/darragh/Hartwig/Data/breast/no_radiation/HARTWIG_merged_mutation_file_july_5_not_1000G

###############################
# Filter out variants in gnomAD
################################
#egrep -v "^#" af-only-gnomad.hg19.vcf | awk '{print $1"\t"$2-1"\t"$2"\t"$4"\t"$5}' > af-only-gnomad.hg19.bed
#bedtools subtract -A -a /home/darragh/Hartwig/Data/breast/no_radiation/HARTWIG_merged_mutation_file_july_5_not_1000G -b af-only-gnomad.hg19.bed > /home/darragh/Hartwig/Data/breast/no_radiation/HARTWIG_merged_mutation_file_july_5_not_1000G_not_gnomad


#################################################################
# Filter out mutations in highly repetitive regions of the genome
#################################################################

#bedtools subtract -A -a /home/darragh/Hartwig/Data/breast/no_radiation/HARTWIG_merged_mutation_file_july_5_not_1000G_not_gnomad -b /home/darragh/Hartwig/Repetitive_dna_filter_vcf.txt > /home/darragh/Hartwig/Data/breast/no_radiation/HARTWIG_merged_mutation_file_july_5_not_1000G_not_gnomad_not_repetitive


###################################################
#Filter out exon mutations
###################################################

#bedtools subtract -A -a /home/darragh/Hartwig/Data/breast/no_radiation/HARTWIG_merged_mutation_file_july_5_not_1000G_not_gnomad_not_repetitive -b /home/darragh/Annotation_reg_elements/hg19_exons.gtf > /home/darragh/Hartwig/Data/breast/no_radiation/HARTWIG_merged_mutation_file_july_5_not_1000G_not_gnomad_not_repetitive_non_coding


#########################################################
#Change to 1 based counting for the activedriverwgs model
#########################################################



################
#Additional filters we could try
################

#Filter out regions not covered by 30x in gnomad
#Filter mismatch repair-related mutations in AAGTTT or complementary nucleotide contexts
#positions with >4 possible matches of the surrounding 36mer to the reference genome (up to 2 mismatches allowed) based on the CRC GEM alignability track
#mutations detected by <2 independent callers / centers in PCAWG (9, 50) (Broad, DKFZ, MuSE, Sanger)


####################
#Filter SCAN-B data
####################

import os
import pybedtools as pybed
import pandas as pd
import re
import numpy as np
import shutil
import glob
import csv
import pysam

#Create individual vcfs for the patients
SCAN_B = pd.read_table("/home/darragh/ICGC/staafVariants.bed",header=0)
SCAN_B.columns = ["chr", "start", "end", "ref", "alt", "gene", "sample"]
SCAN_B["chr"] = 'chr' + SCAN_B["chr"].astype(str)
SCAN_B["start"] = SCAN_B["start"]-1

for sample, df_sample in SCAN_B.groupby('sample'):
    df_sample = df_sample.drop_duplicates()
    df_sample.to_csv(os.path.join("/home/darragh/ICGC/SCANB",sample+".vcf" ), sep='\t', index=False, header=False)

###############
#Filter out samples with >20% of their mutations in highly repetitive regions - print their name and remove manually
###############

repeat_file = pybed.BedTool("/home/darragh/Hartwig/Repetitive_dna_filter_vcf.txt")  


vcf_files = glob.glob('/home/darragh/ICGC/SCANB/*.vcf')

for file in vcf_files:
  mutation_file = pybed.BedTool(file)
  mutation_file_replicative = mutation_file.intersect(repeat_file, wa=True, u=True)
  if len(mutation_file_replicative)/len(mutation_file) > 0.2:
    print(os.path.basename(file))

#None to remove


##############
#Filter out kataegis mutations
##############

vcf_files = glob.glob('/home/darragh/ICGC/SCANB/*.vcf')

for file in vcf_files:
  vcf = pybed.BedTool(file)
  vcf = vcf.cluster(d=1000)
  vcf_df = pd.read_table(vcf.fn,header=None)
  #vcf_df = vcf_df.iloc[:, [0,1,3,4,6,11]]
  vcf_df.columns = ["chr", "start", "end", "ref", "alt", "gene", "sample", "cluster"]
  vcf_df =  vcf_df.groupby('cluster').filter(lambda x : len(x)<6)
  output_file = os.path.join( "/home/darragh/ICGC/SCANB/not_kataegis", os.path.basename(file)).split(".")[0] + ".not.kataegis.bed"
  vcf_df.to_csv(output_file, sep='\t', index=False, header=False)


####################################################
#Create mutation matrix with all remaining mutations
#####################################################

import os
import pybedtools as pybed
import pandas as pd
import re
import numpy as np
import shutil
import glob
import csv
import pysam
pd.set_option('display.max_columns', None)

def parse_unfiltered_vcfs(file):
  mutation_file=set()
  for line in open(file).read().splitlines():
    row=line.split('\t')
    Chr= row[0]
    Start=int(row[1])
    End=int(row[2])
    Ref = row[3]
    Alt = row[4]
    Start=str(Start)
    End=str(End)
    Name = os.path.basename(file)
    mutation_file.add('\t'.join([Chr,Start,End,Ref,Alt,Name]))
  return mutation_file

def run_SCANB():
  if not os.path.exists('/home/darragh/ICGC/SCANB/not_kataegis/SCANB_merged_mutation_file_july_5'):
    SCANB_file_list = glob.glob('/home/darragh/ICGC/SCANB/not_kataegis/*.not.kataegis.bed')
    with open("/home/darragh/ICGC/SCANB/not_kataegis/SCANB_merged_mutation_file_july_5","w") as f:
        for file in SCANB_file_list:
            try:
                mutation_file = parse_unfiltered_vcfs(file)
            except OSError:
                pass
            wr = csv.writer(f,delimiter="\n")
            wr.writerow(mutation_file)

    return

run_SCANB()

##############################
# Filter out variants in 1000G
##############################
#bedtools subtract -A -a /home/darragh/ICGC/SCANB/not_kataegis/SCANB_merged_mutation_file_july_5 -b 1000G_hg19_chr_fixed.bed > /home/darragh/ICGC/SCANB/not_kataegis/SCANB_merged_mutation_file_july_5_not_1000G

###############################
# Filter out variants in gnomAD
################################
#bedtools subtract -A -a /home/darragh/ICGC/SCANB/not_kataegis/SCANB_merged_mutation_file_july_5_not_1000G -b af-only-gnomad.hg19.bed > /home/darragh/ICGC/SCANB/not_kataegis/SCANB_merged_mutation_file_july_5_not_1000G_not_gnomad

##########################################
# Filter out variants in repetitve regions
##########################################
#bedtools subtract -A -a /home/darragh/Hartwig/Data/breast/no_radiation/HARTWIG_merged_mutation_file_july_5_not_1000G_not_gnomad -b /home/darragh/Hartwig/Repetitive_dna_filter_vcf.txt > /home/darragh/Hartwig/Data/breast/no_radiation/HARTWIG_merged_mutation_file_july_5_not_1000G_not_gnomad_not_repetitive


###################################################
#Filter out exon mutations
###################################################

#bedtools subtract -A -a /home/darragh/Hartwig/Data/breast/no_radiation/HARTWIG_merged_mutation_file_july_5_not_1000G_not_gnomad_not_repetitive -b /home/darragh/Annotation_reg_elements/hg19_exons.gtf > /home/darragh/Hartwig/Data/breast/no_radiation/HARTWIG_merged_mutation_file_july_5_not_1000G_not_gnomad_not_repetitive_non_coding







##########################################################################################################################################
# ICGC filters
##########################################################################################################################################

pd.set_option('display.max_columns', None)

#Keep just the whitelisted samples
ICGC_mutation_file = pd.read_table("/home/darragh/ICGC/simple_somatic_mutation.controlled.BRCA_ALL.tsv",header=0)
ICGC_mutation_file = ICGC_mutation_file[["chromosome", "chromosome_start", "chromosome_end", "mutated_from_allele", "mutated_to_allele", "icgc_donor_id", "icgc_sample_id"]]
white_list = pd.read_table("/home/darragh/ICGC/ICGC_white_list_breast.txt",header=None)
white_list.columns = ["icgc_sample_id"]
ICGC_mutation_file_white_list = ICGC_mutation_file.merge(white_list, on='icgc_sample_id', how='inner')

#Remove Kataegis mutations.

#First create individual vcfs for each of the whitelisted patients
for sample, df_sample in ICGC_mutation_file_white_list.groupby('icgc_sample_id'):
    df_sample = df_sample.drop_duplicates()
    df_sample.to_csv(os.path.join("/home/darragh/ICGC/vcfs",sample+".vcf" ), sep='\t', index=False)


###############################
#Filter samples with >20% of mutations in highly repetitive regions
###############################

#Filter out the repeat file to only keep repetitions of 1-5bp with >5 repetitions
#repeat_file.to_csv("/home/darragh/Hartwig/Repetitive_dna_filter_vcf.txt", sep='\t', index=False)
repeat_file = pybed.BedTool("/home/darragh/Hartwig/Repetitive_dna_filter_vcf.txt")
#ls *.vcf |cut -f 2 -d '/' |while read line ; do cat "$line" | sed '1d' | awk '{print "chr"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$7 }'  > "$line".filtered.vcf ; done
vcf_files = glob.glob('/home/darragh/ICGC/vcfs/*.filtered.vcf')

for file in vcf_files:
  mutation_file = pybed.BedTool(file)
  mutation_file = mutation_file.sort()
  mutation_file_replicative = mutation_file.intersect(repeat_file, wa=True, u=True)
  if len(mutation_file_replicative)/len(mutation_file) > 0.2:
    print(os.path.basename(file))

#When all of the filenames are printed to the screen, manually move them into a new filtered folder


##############################
#Filter out Kataegis variants
##############################

import os
import pybedtools as pybed
import pandas as pd
import re
import numpy as np
import shutil
import glob
import csv
import pysam

#ls *.filtered.vcf |cut -f 2 -d '/' |while read line ; do cat "$line" | bedtools sort -i stdin  > "$line".filtered.sorted.vcf ; done
vcf_files = glob.glob('/home/darragh/ICGC/vcfs/*.filtered.sorted.vcf')

for file in vcf_files:
  vcf = pybed.BedTool(file)
  vcf = vcf.cluster(d=1000)
  vcf_df = pd.read_table(vcf.fn,header=None)
  #vcf_df = vcf_df.iloc[:, [0,1,3,4,6,11]]
  vcf_df.columns = ["chr", "start", "end", "ref", "alt", "cluster"]
  #vcf_df['chr'] = 'chr'+str(vcf_df['chr'])
  vcf_df =  vcf_df.groupby('cluster').filter(lambda x : len(x)<6)
  output_file = os.path.join( "/home/darragh/ICGC/vcfs", os.path.basename(file)).split(".")[0] + ".not.kataegis.bed"
  vcf_df.to_csv(output_file, sep='\t', index=False, header=False)

  #come back to this



####################################################################################################################################
#staaf filters
####################################################################################################################################

import os
import pybedtools as pybed
import pandas as pd
import re
import numpy as np
import shutil
import glob
import csv
import pysam

pd.set_option('display.max_columns', None)

STAAF_variants = pd.read_table("/home/darragh/ICGC/staafVariants.bed", header=None)
STAAF_variants.columns = ["chr", "start", "end", "ref", "alt", "null", "patient"]
STAAF_variants["start"] = STAAF_variants["start"]-1

for sample, df_sample in STAAF_variants.groupby('patient'):
    df_sample = df_sample.drop_duplicates()
    df_sample.to_csv(os.path.join("/home/darragh/ICGC/STAAF",sample+".vcf" ), sep='\t', index=False)


###############################
#Filter samples with >20% of mutations in highly repetitive regions
###############################

#Filter out the repeat file to only keep repetitions of 1-5bp with >5 repetitions
#repeat_file.to_csv("/home/darragh/Hartwig/Repetitive_dna_filter_vcf.txt", sep='\t', index=False)
repeat_file = pybed.BedTool("/home/darragh/Hartwig/Repetitive_dna_filter_vcf.txt")
repeat_file = repeat_file.iloc[1: , :]
#ls *.vcf |cut -f 2 -d '/' |while read line ; do cat "$line" | sed '1d' | awk '{print "chr"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$7 }'  > "$line".filtered.vcf ; done
vcf_files = glob.glob('/home/darragh/ICGC/STAAF/*.vcf')

for file in vcf_files:
  mutation_file = pybed.BedTool(file)
  mutation_file = mutation_file.sort()
  mutation_file_replicative = mutation_file.intersect(repeat_file, wa=True, u=True)
  if len(mutation_file_replicative)/len(mutation_file) > 0.2:
    print(os.path.basename(file))