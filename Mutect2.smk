
rule mutect2_1:
        input:
                fasta="/index/hg38/hg38.fa",
                bam="GatheredBamFiles.bam",
                gnomad="/home/darragh/P51A/VQSR/af-only-gnomad.hg38.vcf.gz",
                hg38="/home/darragh/P51A/VQSR/1000g_pon.hg38.vcf.gz",
                intervals="/home/darragh/P51A_Somatic/test/file/{RUN1}_of_13/scattered.bed"
        output:
                a=temp("first_{RUN1}_unfiltered.vcf"),
                c=temp("first_{RUN1}_unfiltered.vcf.stats"),
                b=temp("first_{RUN1}_f1r2.tar.gz")
        shell:
                "gatk Mutect2 -R {input.fasta} -I {input.bam} --germline-resource {input.gnomad} --panel-of-normals {input.hg38} --f1r2-tar-gz {output.b} --geno
type-germline-sites true -O {output.a} --intervals {input.intervals}"

rule mutect2_2:
        input:
                fasta="/index/hg38/hg38.fa",
                bam="GatheredBamFiles.bam",
                gnomad="/home/darragh/P51A/VQSR/af-only-gnomad.hg38.vcf.gz",
                hg38="/home/darragh/P51A/VQSR/1000g_pon.hg38.vcf.gz",
                intervals="/home/darragh/P51A_Somatic/test/file/{RUN2}_of_13/scattered.bed",
        output:
                a=temp("second_{RUN2}_unfiltered.vcf"),
                c=temp("second_{RUN2}_unfiltered.vcf.stats"),
                b=temp("second_{RUN2}_f1r2.tar.gz")
        shell:
                "gatk Mutect2 -R {input.fasta} -I {input.bam} --germline-resource {input.gnomad} --panel-of-normals {input.hg38} --f1r2-tar-gz {output.b} -O {output.a} --genotype-germline-sites true  --intervals {input.intervals}"

rule mutect2_3:
        input:
                fasta="/index/hg38/hg38.fa",
                bam="GatheredBamFiles.bam",
                gnomad="/home/darragh/P51A/VQSR/af-only-gnomad.hg38.vcf.gz",
                hg38="/home/darragh/P51A/VQSR/1000g_pon.hg38.vcf.gz",
                intervals="/home/darragh/P51A_Somatic/test/file/{RUN3}_of_13/scattered.bed",
        output:
                a=temp("third_{RUN3}_unfiltered.vcf"),
                c=temp("third_{RUN3}_unfiltered.vcf.stats"),
                b=temp("third_{RUN3}_f1r2.tar.gz"),

        shell:
                "gatk Mutect2 -R {input.fasta} -I {input.bam} --germline-resource {input.gnomad} --panel-of-normals {input.hg38} --f1r2-tar-gz {output.b} --genotype-germline-sites true -O {output.a} --intervals {input.intervals}"


rule GatherVCFs_1:
        input:
                run_1=expand("first_{directory}_unfiltered.vcf", directory=RUN1),
                run_2=expand("second_{directory}_unfiltered.vcf", directory=RUN2),
                run_3=expand("third_{directory}_unfiltered.vcf", directory=RUN3),
        output:
                "GatheredVariants.vcf",
        run:
                INPUTS = " ".join(["--INPUT {}".format(x) for x in input])
                shell("gatk MergeVcfs {INPUTS} -O GatheredVariants.vcf --CREATE_INDEX true".format(INPUTS=INPUTS))
              
rule MergeMutectStats:
    input:
        chr1="first_1_unfiltered.vcf.stats",
        chr2="first_2_unfiltered.vcf.stats",
        chr3="first_3_unfiltered.vcf.stats",
        chr4="second_4_unfiltered.vcf.stats",
        chr5="second_5_unfiltered.vcf.stats",
        chr6="second_6_unfiltered.vcf.stats",
        chr7="second_7_unfiltered.vcf.stats",
        chr8="third_8_unfiltered.vcf.stats",
        chr9="third_9_unfiltered.vcf.stats",
        chr10="third_10_unfiltered.vcf.stats",
        chr11="third_11_unfiltered.vcf.stats",
        chr12="third_12_unfiltered.vcf.stats",
        chr13="third_13_unfiltered.vcf.stats",
    output:
        "GatheredVariants.vcf.stats",
    shell:
        "gatk MergeMutectStats -stats {input.chr1} -stats {input.chr2} -stats {input.chr3} -stats {input.chr4} -stats {input.chr5} -stats {input.chr6} -stats {input.chr7} -stats {input.chr8} -stats {input.chr9} -stats {input.chr10} -stats {input.chr11} -stats {input.chr12} -stats {input.chr13} -O {output}"

