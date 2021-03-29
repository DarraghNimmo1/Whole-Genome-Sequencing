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

import pysam

configfile: 'config.yaml'

GENOME = config['genome']

IN_BAM = config['in_bam']

THOUSAND_GENOMES_GENOME = config['thousand_genomes_genome']

THOUSAND_GENOMES_DIR = config['thousand_genomes_dir']

hg38_DICTIONARY = config['hg38_dictionary']

directory_function = functools.partial(os.path.join, config['results'])
INTERVALS_DIR = directory_function('Intervals')
COUNTS_DIR  = directory_function('Counts')
PON_DIR = directory_function('PON')
DENOISED_DIR = directory_function('Denoised')
PLOT_DENOISED_DIR = directory_function('Plot_Denoised')
ALLELIC_COUNTS_DIR = directory_function('Allelic_Counts')
MODEL_SEGMENTS_DIR = directory_function('Model_Segments')
CALL_COPY_RATIOS_DIR = directory_function('Call_Copy_Ratios')

################################################################################
##Rules
################################################################################

rule all:
    input:
        

rule get_intervals:
    input:
        fasta = GENOME
    output:
        list = os.path.join(INTERVALS_DIR,  '{sample}.interval_list') #sample will be specified at the command line
    message:
        "Making an interval list. These are the regions which will be used to search for copy number gains or losses."
  shell:
      """
      gatk PreprocessIntervals -R {input.fasta} --bin-length 1000 --padding 0 -O {output}
      """

rule annotate_intervals_GC:
    input:
        intervals = rules.get_intervals.output.list
        genome = GENOME
    output:
        annotated_interval = os.path.join(PON_DIR, '{GC_annotated_intervals.tsv')
    message:
        "Annotating hg38 genome intervals with GC percentages"
    shell:
        """
        gatk AnnotateIntervals -R {input.genome} -L {input.intervals} --interval-merging-rule OVERLAPPING_ONLY -O {output.annotated_interval}
        """
## A seperate pon is created here for male and female 1000 genome samples. This workflow is senstive to the sex of individuals in the pon and it is important to use a sex 
## matched pon to denoise read counts.
rule convert_bam_index:
    input:
        cram = os.path.join(THOUSAND_GENOMES_DIR, '{males}.final.cram')
        genome = THOUSAND_GENOMES_GENOME
    output:
        bam = os.path.join(THOUSAND_GENOMES_DIR, '{males}.bam')
        index = os.path.join(THOUSAND_GENOMES_DIR, '{males}.bam.bai')
    message:
        "Converting male thousand genomes crams to bams and indexing. 50 to do, might take a while..."
    shell:
        """
        samtools view -@ 48 -b -T {input.genome} {input.cram} > {output.bam}; samtools index -@ 48 {output.bam}  
        """

rule collect_read_counts:   #do a rule all for this!! use expand and wildcards 
    input:
        bam = rules.convert_bam_index.output.bam
        intervals = rules.get_intervals.output.list
    output:
        counts = os.path.join(PON_DIR, '{males}.hdf5')
    message:
        "Collecting read counts for the male thousand genome bams."
    shell:
        """
        gatk CollectReadCounts -I {input.bam} -L {input.intervals} --interval-merging-rule OVERLAPPING_ONLY -O {output.counts}
        """        

rule create_pon:
    input:
        hdf5s = expand(os.path.join(PON_DIR, '{males}.hdf5'))
        GC_percentages = rules.annotate_intervals_GC.output.annotated_interval
    output:
        PON = os.path.join(PON_DIR, 'pon.hdf5')
    message:
        "Creating a PON from the thousand genome samples"
    shell:
        """
        gatk CreateReadCountPanelOfNormals -I {input.hdf5s} --annotated-intervals {input.GC_percentages} -O {output.PON}
        """

#rule convert_bam_index_females:
#    input:
#        cram = os.path.join(THOUSAND_GENOMES_DIR, '{females}.final.cram')
#        genome = THOUSAND_GENOMES_GENOME
#    output:
#        bam = os.path.join(THOUSAND_GENOMES_DIR, '{females}.bam')
#        index = os.path.join(THOUSAND_GENOMES_DIR, '{females}.bam.bai')
#    message:
#        "Converting female thousand genomes crams to bams and indexing. 50 to do, might take a while..."
#    shell:
#        """
#        samtools view -@ 48 -b -T {input.genome} {input.cram} > {output.bam}; samtools index -@ 48 {output.bam}
#        """

#rule collect_read_counts_females: #do a rule all for this!! use expand and wildcards
#    input:
#        bam = rules.convert_bam_index_females.output.bam
#        intervals = rules.get_intervals.output.list
#    output:
#        counts = os.path.join(PON_DIR, '{females}.hdf5')
#    message:
#        "Collecting read counts for the female thousand genome bams."
#    shell:
#        """
#        gatk CollectReadCounts -I {input.bam} -L {input.intervals} --interval-merging-rule OVERLAPPING_ONLY -O {output.counts}
#        """

#rule create_female_pon:
#    input:
#        hdf5s = expand(os.path.join(PON_DIR, '{females}.hdf5'))
#        GC_percentages = rules.annotate_intervals_GC.output.annotated_interval
#    output:
#        PON = os.path.join(PON_DIR, 'female.pon.hdf5')
#    message:
#        "Creating a PON for female thousand genome samples"
#    shell:
#        """
#        gatk CreateReadCountPanelOfNormals -I {input.hdf5s} --annotated-intervals {input.GC_percentages} -O {output.PON}
#        """

rule read_counts_sample:
    input:
        bam = os.path.join(IN_BAM, '{sample}.bam')
        intervals = rules.get_intervals.output.list
    output:
        counts = os.path.join(COUNTS_DIR, {sample}.counts.hdf5)
    message:
        "Collecting read counts across the interval list in the sample bam."
    run:
        bamfile = pysam.AlignmentFile({input.bam}, "rb")
        shell("gatk CollectReadCounts -I {input.bam} -L {input.intervals} --interval-merging-rule OVERLAPPING_ONLY -O {output.counts}")

rule denoise_read_counts:
    input:
        counts = rules.read_counts_sample.output.counts
        pon = rules.create_pon.output.PON
    output:
        standard_copy_ratios = os.path.join(DENOISED_DIR, '{sample}.standardizedCR.tsv')
        denoised_copy_ratios = os.path.join(DENOISED_DIR, '{sample}.denoisedCR.tsv')
    message:
        "Denoising read counts using pre-made panel of normals"
    shell:
        """
        gatk DenoiseReadCounts -I {input.counts} --count-panel-of-normals {input.pon} --standardized-copy-ratios {output.standard_copy_ratios} --denoised-copy-ratios {output.denoised_copy_ratios}
        """

rule plot_denoised_ratios:
    input:
        standard_copy_ratios = rules.read_counts_sample.output.standard_copy_ratios
        denoised_copy_ratios = rules.read_counts_sample.output.denoised_copy_ratios
        dict = hg38_DICTIONARY
    output:
        plot_1 = os.path.join(PLOT_DENOISED_DIR, '{sample}.denoised.png')
        plot_2 = os.path.join(PLOT_DENOISED_DIR, '{sample}.denoisedLimit4.png')
        stndrd_mad = os.path.join(PLOT_DENOISED_DIR, '{sample}.standardizedMAD.txt')
        denoised_mad = os.path.join(PLOT_DENOISED_DIR, '{sample}.denoisedMAD.txt')
        delta_mad = os.path.join(PLOT_DENOISED_DIR, '{sample}.deltaMAD.txt')
        scaled_delta_mad = os.path.join(PLOT_DENOISED_DIR, '{sample}.scaledDeltaMAD.txt')
    params:
        prefix = {sample}
        out_dir = PLOT_DENOISED_DIR
    shell:
        """
        gatk PlotDenoisedCopyRatios --standardized-copy-ratios {input.standard_copy_ratios} --denoised-copy-ratios {input.denoised_copy_ratios} --sequence-dictionary {input.hg38_DICTIONARY} --output-prefix {params.prefix} -O {params.out_dir}
        """
        
rule collect_allelic_counts:
    input:
        bam = os.path.join(IN_BAM, '{sample}.bam')
        genome = GENOME
        sites = rules.get_intervals.output.list
    output:
        counts = os.path.join(ALLELIC_COUNTS_DIR, '{sample}.allelicCounts.tsv')
    shell:
        """
        gatk CollectAllelicCounts -I {input.bam} -R {input.genome} -L {input.sites} -O {output.counts}
        """

rule model_segments:
    input:
        denoised_copy_ratios = rules.denoise_read_counts.output.denoised_copy_ratios,
        allelic_counts = rules.collect_allelic_counts.output.counts
    output:
        seg_file = os.path.join(MODEL_SEGMENTS_DIR, '{sample}.cr.seg')
    params:
        dir = MODEL_SEGMENTS_DIR
        prefix = '{sample}'
    shell:
        """
        gatk ModelSegments --denoised-copy-ratios {input.denoised_copy_ratios} --allelic-counts {input.counts} --output {params.dir} --output-prefix {params.prefix} --minimum-total-allele-count-case 20
        """

rule call_copy_ratio_segments:
    input:
        segments = rules.model_segments.output.seg_file
    output:
        called_segs = os.path.join(CALL_COPY_RATIOS_DIR, '{sample}_called.seg')
    shell:
        """
        gatk CallCopyRatioSegments --input hcc1143_T_clean.cr.seg --output sandbox/hcc1143_T_clean.called.seg
        """

##To do: plot copy ratio calls.

        
