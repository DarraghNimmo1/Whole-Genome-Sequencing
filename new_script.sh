zcat V300068040_L01_405_1.fq.gz V300068040_L01_406_1.fq.gz V300068040_L01_407_1.fq.gz V300068040_L01_408_1.fq.gz V300073024_L01_405_1.fq.gz V300073024_L01_406_1.fq.gz V300073024_L01_407_1.fq.gz V300073024_L01_408_1.fq.gz > P1A_1.fq

zcat V300068040_L01_405_1.fq.gz V300068040_L01_406_1.fq.gz V300068040_L01_407_1.fq.gz V300068040_L01_408_1.fq.gz V300073024_L01_405_1.fq.gz V300073024_L01_406_1.fq.gz V300073024_L01_407_1.fq.gz V300073024_L01_408_1.fq.gz > P1A_2.fq

java -Xmx8G -jar /home/darragh/picard.jar FastqToSam FASTQ=P1A_1.fq FASTQ2=P1A_2.fq OUTPUT=P1A_fastqtosam.bam READ_GROUP_NAME=V300 SAMPLE_NAME=P1A LIBRARY_NAME=V300 PLATFORM=illumina

java -Xmx8G -jar /home/darragh/picard.jar MarkIlluminaAdapters I=P1A_fastqtosam.bam O=P1A_markilluminaadapters.bam M=P1A_markilluminaadapters_metrics.txt TMP_DIR=/home/darragh/PDAC_project/data/P1A/tmp/ 

java -Xmx8G -jar /home/darragh/picard.jar SamToFastq I=P1A_markilluminaadapters.bam FASTQ=/dev/stdout CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true |  bwa mem -M -t 7 -p /index/hg38/hg38.fa /dev/stdin | java -Xmx16G -jar /home/darragh/picard.jar MergeBamAlignment ALIGNED_BAM=/dev/stdin UNMAPPED_BAM=P1A_fastqtosam.bam OUTPUT=P1A_piped.bam R=/index/hg38/hg38.fa CREATE_INDEX=true ADD_MATE_CIGAR=true CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS --SORT_ORDER queryname

gatk MarkDuplicates -I P1A_piped.bam  -O P1A_marked_duplicates.bam --REMOVE_DUPLICATES true 

java -jar /home/darragh/picard.jar SortSam I=P1A_marked_duplicates.bam O=P1A_marked_duplicates.sorted.bam SORT_ORDER=coordinate

BaseRecalibrator -I P1A_marked_duplicates.sorted.bam -R /index/hg38/hg38.fa --known-sites /home/darragh/GATK/resources/dbSNP_hg38_00-All_named.vcf.gz -O P1A_recal_data.table

gatk ApplyBQSR -R /index/hg38/hg38.fa -I P1A_marked_duplicates.sorted.bam --bqsr-recal-file P1A_recal_data.table -O P1A_BQSR.bam

gatk Mutect2 -R /index/hg38/hg38.fa -I P1A_BQSR.bam --germline-resource /home/darragh/GATK/resources/af-only-gnomad.hg38.vcf.gz --panel-of-normals /home/darragh/GATK/resources/1000g_pon.hg38.vcf.gz -O P1A_MUTECT2.vcf.gz --af-of-alleles-not-in-resource 5e-8 --f1r2-tar-gz f1r2.tar.gz --native-pair-hmm-threads 48

gatk GetPileupSummaries -I P1A_BQSR.bam -V /home/darragh/GATK/resources/af-only-gnomad.hg38.vcf.gz -L /home/darragh/GATK/resources/af-only-gnomad.hg38.vcf.gz -O P1A_pileups.table
