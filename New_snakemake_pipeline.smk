#Write a snakemake script to analyse WGS data for PDAC project.

################################################
#WGS analysis snakemake script
#Darragh Nimmo
#April 2022
#################################################

import os

import functools

configfile: 'config.yaml'
  
FILE_NAMES = config['file_names']

SAMPLE = config['sample']

READS_DIR = config['reads']

GENOME = config['genome']

PLATFORM = config['platform']

PICARD = config['picard']

REGIONS = config['regions']

GNOMAD = config['gnomad']

DBSNP = config['dbsnp']


directory_function = functools.partial(os.path.join, config['results'])
BAM_DIR = directory_function('Bam')
VCF_DIR = directory_function('VCF')


rule all:
    input:
        expand()
        
#Running these steps for all seperate read group fastq files

rule fastqtosam:
    input:
        forward_read = os.path.join(READS_DIR, '{file_names}_1.fq.gz'),
        reverse_read = os.path.join(READS_DIR, '{file_names}_2.fq.gz')
    output:
        bam = os.path.join(BAM_DIR, '{file_names}_fastqtosam.bam')
    params:
        picard = PICARD,
        platform = PLATFORM,
        sample_name = SAMPLE,
        read_group = SAMPLE + '_' + {file_name}[11:18],
        library_name = {file_name}[15:18]
    shell:
        """
        java -Xmx377G -jar {params.picard} FastqToSam FASTQ={input.forward_read} FASTQ2={input.reverse_read} OUTPUT={output.bam} \
        READ_GROUP_NAME={params.read_group} SAMPLE_NAME={params.sample_name} LIBRARY_NAME={params.library_name} PLATFORM={params.platfrom}
        """
        
rule markilluminaadapters:
    input:
        rules.fastqtosam.output.bam
    output:
        bam = os.path.join(BAM_DIR, '{file_names}_markilluminaadapters.bam')
        txt = os.path.join(BAM_DIR, '{file_names}_markilluminaadapters.txt')
    params:
        picard = PICARD
    shell:
        """
        java -Xmx377G -jar {params.picard} MarkIlluminaAdapters I={input.bam} O={output.bam} M={output.txt} 
        """

        
rule alignment:
    input:
        bam = rules.markilluminaadapters.output.bam,
        unmapped_bam = rules.fastqtosam.output.bam 
    output:
        bam = os.path.join(BAM_DIR, '{file_names}_piped.bam')
        
    params:
        picard = PICARD
        genome = GENOME
    shell:
        """
        java -Xmx377G -jar {params.picard} SamToFastq I={input.bam} FASTQ=/dev/stdout CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true \
        |  bwa mem -M -t 7 -p {params.genome} /dev/stdin \
        | java -Xmx377G -jar {params.picard} MergeBamAlignment ALIGNED_BAM=/dev/stdin UNMAPPED_BAM={output.unmapped_bam} OUTPUT={output.bam} R={params.genome} \
        CREATE_INDEX=true ADD_MATE_CIGAR=true CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 \
        PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS --SORT_ORDER queryname
        """
        
rule subset_alignments:
    input:
        bam = rules.alignment.output.bam
    output:
        bam = os.path.join(BAM_DIR, '{file_names}_subsetted.bam')
    params:
        regions = REGIONS
    shell:
        """
        bedtools intersect -a {input.bam} -b {params.regions} > {output.bam}
        """
        
rule mark_duplicates:
    input:
        bam = expand(os.path.join(BAM_DIR, '{file_names}_subsetted.bam'), file_names = FILE_NAMES)
    output:
        bam = os.path.join(BAM_DIR, '{sample}_marked_duplicates.bam'),
        txt = os.path.join(BAM_DIR, '{sample}_marked_duplicates.txt')
    params:
        picard = PICARD
    shell:
        """
        gatk MarkDuplicates -I {input.bam}  -O {output.bam} -M {output.txt}  --REMOVE_DUPLICATES true 
        """
 
rule sortsam:
    input:
        bam = rules.mark_duplicates.output.bam
    output:
        bam = os.path.join(BAM_DIR, '{sample}_marked_duplicates.sorted.bam')
    params:
        picard = PICARD
    shell:
        """
        java -Xmx377G -jar {params.picard} SortSam I={input.bam} O={output.bam} SORT_ORDER=coordinate
        """
        
rule baserecalibrator:
    input:
        bam = rules.sortsam.output.bam
    output:
        table = os.path.join(BAM_DIR, '{sample}_recal_data.table')
    params:
        gnomad = GNOMAD,
        dbsnp = DBSNP
        genome = GENOME
    shell:
      
        """
        gatk BaseRecalibrator -I {input.bam} -R {params.genome} --known-sites {params.dbsnp} --known-sites {params.gnomad} -O {output.table}
        """
        
rule applybqsr:
    input:
        bam = rules.sortsam.output.bam,
        table = rules.baserecalibrator.output.bam
    output:
        bam = os.path.join(BAM_DIR, '{sample}_BQSR.bam')
    params:
        genome = GENOME
    shell:
        """
        gatk ApplyBQSR -R {params.genome} -I {input.bam} --bqsr-recal-file {input.table} -O {output.bam}
        """
        
rule mutect2:
    input:
        bam = rules.applybqsr.output.bam
    output:
        vcf = os.path.join(VCF_DIR, '{sample}_MUTECT2.vcf.gz')    
    params:
        genome = GENOME,
        gnomad = GNOMAD,
        pon = PON
    shell:
        """
        gatk Mutect2 -R {params.genome} -I {input.bam} --germline-resource {input --panel-of-normals /home/darragh/GATK/resources/1000g_pon.hg38.vcf.gz -O P1A_MUTECT2.vcf.gz --af-of-alleles-not-in-resource 5e-8 --f1r2-tar-gz f1r2.tar.gz --native-pair-hmm-threads 48
        """
        

     
        
