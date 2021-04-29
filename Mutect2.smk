split = [1,2,3,4,5,6,7,8,9,10,11,12,13]
SAMPLE = "P4A"


rule all:
        input:
            SAMPLE+"_funcotated.vcf"


rule mutect2_1:
    input:
        fasta="/index/hg38/hg38.fa",
        bam=SAMPLE+".bam",
        gnomad="/home/darragh/P51A/VQSR/af-only-gnomad.hg38.vcf.gz",
        hg38="/home/darragh/P51A/VQSR/1000g_pon.hg38.vcf.gz",
        intervals="/home/darragh/P51A_Somatic/test/file/{split}_of_13/scattered.bed"
    output:
        a="{split}_unfiltered.vcf",
        c="{split}_unfiltered.vcf.stats",
        b="{split}_f1r2.tar.gz"
    shell:
        """
        gatk Mutect2 -R {input.fasta} -I {input.bam} --germline-resource {input.gnomad} --panel-of-normals {input.hg38} --f1r2-tar-gz {output.b} --genotype-germline-sites true -O {output.a} --intervals {input.intervals}
        """


rule GatherVCFs_1:
    input:
        expand("{directory}_unfiltered.vcf", directory=split),
    output:
        "P4A_gathered_mutect.vcf",
    run:
        INPUTS = " ".join(["--INPUT {}".format(x) for x in input])
        shell("gatk MergeVcfs {INPUTS} -O P4A_gathered_mutect.vcf --CREATE_INDEX true".format(INPUTS=INPUTS))
rule MergeMutectStats:
    input:
        one = "1_unfiltered.vcf.stats",
        two = "2_unfiltered.vcf.stats",
        three = "3_unfiltered.vcf.stats",
        four = "4_unfiltered.vcf.stats",
        five = "5_unfiltered.vcf.stats",
        six = "6_unfiltered.vcf.stats",
        seven = "7_unfiltered.vcf.stats",
        eight = "8_unfiltered.vcf.stats",
        nine = "9_unfiltered.vcf.stats",
        ten = "10_unfiltered.vcf.stats",
        eleven = "11_unfiltered.vcf.stats",
        twelve = "12_unfiltered.vcf.stats",
        thirteen = "13_unfiltered.vcf.stats",

    output:
        stats = SAMPLE+"_gathered_mutect.vcf.stats",
    shell:
        """
        gatk MergeMutectStats -stats {input.one} -stats {input.two} -stats {input.three} -stats {input.four} -stats {input.five} -stats {input.six} -stats {input.seven} -stats {input.eight} -stats {input.nine} -stats {input.ten} -stats {input.eleven} -stats {input.twelve} -stats {input.thirteen} -O {output.stats}
        """
rule learn_orien_model:
    input:
        one = "1_f1r2.tar.gz",
        two = "2_f1r2.tar.gz",
        three = "3_f1r2.tar.gz",
        four = "4_f1r2.tar.gz",
        five  = "5_f1r2.tar.gz",
        six = "6_f1r2.tar.gz",
        seven = "7_f1r2.tar.gz",
        eight = "8_f1r2.tar.gz",
        nine = "9_f1r2.tar.gz",
        ten = "10_f1r2.tar.gz",
        eleven = "11_f1r2.tar.gz",
        twelve = "12_f1r2.tar.gz",
        thirteen = "13_f1r2.tar.gz",
    output:
        model = SAMPLE+"_all_read_orien_model.tar.gz"
    shell:
        """
        gatk LearnReadOrientationModel -I {input.one} -I {input.two} -I {input.three} -I {input.four} -I {input.five} -I {input.six} -I {input.seven} -I {input.eight} -I {input.nine} -I {input.ten} -I {input.eleven} -I {input.twelve} -I {input.thirteen} -O {output.model}
        """

rule get_pile_up_summaries:
    input:
        bam=SAMPLE+".bam",
        gnomad="/home/darragh/P51A_Somatic/BWA/test/af-only-gnomad.hg38.vcf.gz",
        int="/home/darragh/P51A_Somatic/BWA/test/wgs_calling_regions.hg38.interval_list"
    output:
        table = SAMPLE+"_pileups.table"
    shell:
        "gatk GetPileupSummaries -I {input.bam} -V {input.gnomad} -L {input.int} -O {output.table}"

rule calculate_contamination:
    input:
        rules.get_pile_up_summaries.output.table
    output:
        contam= SAMPLE+"_calculatecontamination.table"
    shell:
        "gatk CalculateContamination -I {input} -O {output.contam}"

rule filter_mutect_calls:
    input:
        variants=rules.GatherVCFs_1.output,
        model=rules.learn_orien_model.output.model,
        stats=rules.MergeMutectStats.output.stats,
        fasta="/index/hg38/hg38.fa",
        contam=rules.calculate_contamination.output.contam,
        intervals="/home/darragh/P51A_Somatic/test/file/{directory}_of_13/scattered.bed",
    output:
        filtered = "variants_{directory}_filtered.vcf"
    shell:
         "gatk FilterMutectCalls -V {input.variants} -R {input.fasta} --contamination-table {input.contam} --ob-priors {input.model} -O {output.filtered} --intervals {input.intervals}"

rule GatherVCFs_2:
        input:
                expand("variants_{directory}_filtered.vcf", directory=split)
        output:
                filtered = SAMPLE+"_gathered_filtered.vcf"
        run:
                INPUTS = " ".join(["--INPUT {}".format(x) for x in input])
                shell("gatk MergeVcfs {INPUTS} -O P4A_gathered_filtered.vcf --CREATE_INDEX true".format(INPUTS=INPUTS))

rule funcotate:
    input:
        vcf = rules.GatherVCFs_2.output.filtered,
        fasta = "/index/hg38/hg38.fa"
    output:
        funcotated = SAMPLE+"_funcotated.vcf"
    params:
        data = "/home/darragh/funcotator_dataSources.v1.7.20200521s/use"
    shell:
        """
        gatk Funcotator -R {input.fasta} -V {input.vcf} -O {output.funcotated}  --output-file-format VCF --data-sources-path {params.data} --ref-version hg38 --exclude-field Gencode_34_ncbiBuild --exclude-field Gencode_34_chromosome --exclude-field Gencode_34_start --exclude-field Gencode_34_end --exclude-field Gencode_34_refAllele --exclude-field Gencode_34_tumorSeqAllele1 --exclude-field Gencode_34_tumorSeqAllele2 --exclude-field Gencode_34_annotationTranscript --exclude-field Gencode_34_transcriptStrand --exclude-field Gencode_34_transcriptExon --exclude-field Gencode_34_transcriptPos --exclude-field Gencode_34_codonChange --exclude-field Gencode_34_gcContent --exclude-field Gencode_34_referenceContext --exclude-field Gencode_34_otherTranscripts
        """
