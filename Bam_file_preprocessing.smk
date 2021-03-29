RUN1=[1,2,3]

RUN2=[4,5,6,7]

RUN3=[8,9,10,11,12,13]

RUN4 = [1,2,3,4,5,6,7,8,9,10,11,12,13]

rule all:
    input:
        expand("/home/darragh/P51A_Somatic/test/file/{directory}_of_13/scattered.bed", directory=RUN1),
        expand("/home/darragh/P51A_Somatic/test/file/{directory}_of_13/scattered.bed", directory=RUN2),
        expand("/home/darragh/P51A_Somatic/test/file/{directory}_of_13/scattered.bed", directory=RUN3),
        expand("/home/darragh/P51A_Somatic/test/file/{directory}_of_13/scattered.bed", directory=RUN3),

rule BwaMem:
        input:
                fasta = "/index/hg38/hg38.fa",
                read1 = "V300_R1.fq",
                read2 = "V300_R2.fq"
        output:
                temp("V300_mapped.bam")
        shell:
                "bwa mem -t 24 {input.fasta} {input.read1} {input.read2}"
                " | samtools view -@ 48 -b | samtools sort -@ 48 -o {output}"


rule add_readgroups:
    input:
        "V300_mapped.bam"
    output:
        temp("V300_mapped_RG.bam")
    shell:
        "java -jar /home/darragh/picard.jar AddOrReplaceReadGroups CREATE_INDEX=true"
        " I={input} O={output} RGID=4 RGLB=lib1 RGPL=BGI RGPU=unit1 RGSM=20 "

rule MarkDupSpark:
    input:
        bam="V300_mapped_RG.bam"
    output:
        bam=temp("markedDuplicates.bam"),
        metrics=temp("markedDuplicates.txt"),
    shell:
        "java -jar /home/darragh/picard.jar MarkDuplicates -I {input.bam} -O {output.bam} -M {output.metrics} --REMOVE_DUPLICATES true"

rule index_md_bam:
        input:
                "markedDuplicates.bam"
        output:
                temp("markedDuplicates.bai")
        shell:
                "java -jar /home/darragh/picard.jar BuildBamIndex I={input}"

rule BaseRecalibrator_1:
        input:
                bam="markedDuplicates.bam",
                index="markedDuplicates.bai",
                fasta="/index/hg38/hg38.fa",
                dbsnp="/home/darragh/P51A/VQSR/Homo_sapiens_assembly38.dbsnp138.vcf",
                gnomad="/home/darragh/P51A_Somatic/BWA/test/af-only-gnomad.hg38.vcf.gz",
                intervals="/home/darragh/P51A_Somatic/test/file/{RUN1}_of_13/scattered.bed"
        output:
                temp("BQSR1_{RUN1}.table"),
        threads:5
        shell:
                "gatk BaseRecalibrator --reference {input.fasta} --input {input.bam} -O {output} --known-sites {input.dbsnp} --known-sites {input.gnomad} --intervals {input.intervals}"
            
 rule BaseRecalibrator_2:
        input:
                bam="markedDuplicates.bam",
                index="markedDuplicates.bai",
                fasta="/index/hg38/hg38.fa",
                dbsnp="/home/darragh/P51A/VQSR/Homo_sapiens_assembly38.dbsnp138.vcf",
                gnomad="/home/darragh/P51A_Somatic/BWA/test/af-only-gnomad.hg38.vcf.gz",
                intervals="/home/darragh/P51A_Somatic/test/file/{RUN2}_of_13/scattered.bed"
        output:
                temp("BQSR2_{RUN2}.table"),
        threads:5

        shell:
                "gatk BaseRecalibrator --reference {input.fasta} --input {input.bam} -O {output} --known-sites {input.dbsnp} --known-sites {input.gnomad} --intervals {input.intervals}"

rule BaseRecalibrator_3:
        input:
                bam="markedDuplicates.bam",
                index="markedDuplicates.bai",
                fasta="/index/hg38/hg38.fa",
                dbsnp="/home/darragh/P51A/VQSR/Homo_sapiens_assembly38.dbsnp138.vcf",
                gnomad="/home/darragh/P51A_Somatic/BWA/test/af-only-gnomad.hg38.vcf.gz",
                intervals="/home/darragh/P51A_Somatic/test/file/{RUN3}_of_13/scattered.bed"
        output:
                temp("BQSR3_{RUN3}.table"),
        threads:5

        shell:
                "gatk BaseRecalibrator --reference {input.fasta} --input {input.bam} -O {output} --known-sites {input.dbsnp} --known-sites {input.gnomad} --intervals {input.intervals}"
rule GatherBQSR:
        input:
            expand("BQSR1_{directory}.table", directory=RUN1),
            expand("BQSR2_{directory}.table", directory=RUN2),
            expand("BQSR3_{directory}.table", directory=RUN3)
        output:
            "GatheredBQSR.table"
        run:
            INPUTS = " ".join(["--input {}".format(x) for x in input])
            shell("gatk GatherBQSRReports {INPUTS} -O GatheredBQSR.table".format(INPUTS=INPUTS))

rule apply_bqsr:
    input:
        table="GatheredBQSR.table",
        fasta="/index/hg38/hg38.fa",
        bam="markedDuplicates.bam",
        index="markedDuplicates.bai",

    output:
        bam = "GatheredBamFiles.bam",

    threads:5
    shell:
        "gatk ApplyBQSR -R {input.fasta} -I {input.bam} --bqsr-recal-file {input.table} -O {output} --create-output-bam-index true "

         
