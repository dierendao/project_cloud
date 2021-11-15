configfile: "config/config.yaml",

SRA,FRR=glob_wildcards("results/subset/{sra}_{frr}.fastq")

rule all:
        input:
                expand("results/fastqc_init/{sample}_fastqc.zip", sample = config["samples"]),
                expand("results/fastqc_init/{sample}_fastqc.html", sample = config["samples"]),
                expand("results/trimmedReads/{sra}_1P.{extension}", sra=SRA, frr=FRR, extension=["fastq"]),
                expand("results/trimmedReads/{sra}_2P.{extension}", sra=SRA, frr=FRR, extension=["fastq"]),
                expand("results/trimmedReads/{sra}_1U.{extension}", sra=SRA, frr=FRR, extension=["fastq"]),
                expand("results/trimmedReads/{sra}_2U.{extension}", sra=SRA, frr=FRR, extension=["fastq"]),
                expand("results/fastqc_post/{sra}_fastqc.{extension}", sra=SRA, extension=["zip","html"]),
                expand("results/qc/{sra}_fastqc.{extension}", sra=SRA, extension=["html"]),
                expand("results/bowtie2/{sra}trim.mapped.sorted.bam", sra=SRA),
                expand("results/picard/{sra}_remove.dup.bam", sra=SRA),
                expand("results/picard/{sra}Picard.metrics.txt", sra=SRA),
                expand("results/GC_remove/{sra}_gc_corrected.bam", sra=SRA),
                expand("results/GC_remove/{sra}_freq.txt", sra=SRA),
                expand("results/npz/{sra}_deeptool.npz", sra=SRA),
                expand("results/expl_data_correlation/{sra}_scatterplot_PearsonCorr_bigwigScores.png", sra=SRA),
                expand("results/expl_data_correlation/{sra}_scatterplot_PearsonCorr_bigwigScores.tab", sra=SRA),
                expand("results/expl_data_correlation/{sra}_heatmap_SpearmanCorr_readCounts.png", sra=SRA),
                expand("results/expl_data_correlation/{sra}_heatmap_SpearmanCorr_readCounts.tab", sra=SRA)

rule unzip:
        input:
                lambda wildcards: config["samples"][wildcards.sample]
        output:
                "results/subset/{sample}.fastq"
        shell:
                """
                mkdir -p results/subset
                gunzip -c {input} > {output}
                """
rule fastqc_init:
        input:
                "results/subset/{sample}.fastq"
        output:
                "results/fastqc_init/{sample}_fastqc.html",
                "results/fastqc_init/{sample}_fastqc.zip"
        conda:
                "env.yaml"
        threads: 2
        shell:
                """
                mkdir -p results/fastqc_init
                fastqc {input} -o "results/fastqc_init"
                """

rule trimmomatic:
        input:
                read1="results/subset/{sra}_1.fastq",
                read2="results/subset/{sra}_2.fastq"
        output:
                forwardPaired="results/trimmedReads/{sra}_1P.fastq",
                forwardUnpaired="results/trimmedReads/{sra}_1U.fastq",
                reversePaired="results/trimmedReads/{sra}_2P.fastq",
                reverseUnpaired="results/trimmedReads/{sra}_2U.fastq"
        threads:
                8
        params:
                basename="results/trimmedReads/{sra}.fastq",
                log="results/trimmedReads/{sra}.log"
        shell:
                """
                trimmomatic PE -threads {threads} {input.read1} {input.read2} \
                -baseout {params.basename} \
                ILLUMINACLIP:data/mydatalocal/atacseq/NexteraPE-PE.fa:2:30:10:>
                LEADING:3 TRAILING:3 MINLEN:33 SLIDINGWINDOW:4:15 2>{params.lo>
                """

rule fastqc_post:
        input:
                forwardPaired=rules.trimmomatic.output.forwardPaired
                forwardUnpaired=rules.trimmomatic.output.forwardUnpaired
                reversePaired=rules.trimmomatic.output.reversePaired
                reverseUnpaired=rules.trimmomatic.output.reverseUnpaired
        output:
                html="results/fastqc_post/{sra}_fastqc.html",
                zip="results/fastqc_post/{sra}_fastqc.zip"
        log:
                "results/fastqc_post/{sra}_fastqc.log"
        conda:
                "env.yaml"
        threads: 8
        shell:
                """
                mkdir -p results/fastqc_post
                fastqc {input} -o "results/fastqc_post 2> {log}" 
                """

rule multiqc:
        input:
                rules.fastqc_post.output
        output:
                "results/qc/{sra}_multiqc.html"
        log:
                "results/qc/{sra}_multiqc.out"
        conda:
                "env.yaml"
        threads:
                8
        shell:
                "multiqc --threads {threads}"
                " -o {output}"
                " -i {input}"
                " 2> {log}"

rule bowtie2:
        input:
                read1=rules.trimmomatic.output.forwardPaired
                read2=rules.trimmomatic.output.reversePaired
        output:
                bam="results/bowtie2/{sra}trim.mapped.sorted.bam"
        log:
                "results/bowtie2/{sra}Log.final.out"
        params:
                index="/data/mydatalocal/atacseq/databank/bowtie2/all"
        conda:
                "env.yaml"
        threads: 8
        shell:
                "bowtie2 --threads {threads} --very-sensitive "
                "-x {params.index} -1 {input.read1} -2 {input.read2} "
                "| samtools view -Sbh -o {snakemake.output} 2> {log}"

rule picard:
        input:
                rules.bowtie2.output.bam
        output:
                bam="results/picard/{sra}_remove.dup.bam",
                metrics= "results/picard/{sra}Picard.metrics.txt"
        log:
                "results/bowtie2/{sra}Picard.remove.out"
        params:
                extra="REMOVE_DUPLICATES=true"
        conda:
                "env.yaml"
        threads: 8
        shell:
                "java -jar picard.jar MarkDuplicates --threads {threads} " 
                "{extra} "  # User defined parmeters
                "{input} "  # Input bam(s)
                "OUTPUT={output.bam} "  # Output bam
                "METRICS_FILE={output.metrics} "  # Output metrics
                "2> {log}" 

rule  GC_remove:
        input: 
                rules.picard.output.bam
        output:
                bam="results/GC_remove/{sra}_gc_corrected.bam"
                freq="results/GC_remove/{sra}_freq.txt"
        log:
                "results/GC_remove/{sra}_gc_corrected.out"
        params:
                genome="/home/users/shared/databanks/bio/ncbi/genomes/Mus_musculus/Mus_musculus_GRCm38.p6/Mus_musculus_2020-7-9/2bit/all.2bit",
                freq="results/GC_remove/{sra}_freq.txt"
        conda:
                "env.yaml"
        threads:
                8
        shell:
                #Part I: Detecting GC Biases
                "computeGCBias --threads {threads} "
                "-b {input} "
                "--effectiveGenomeSize 2652783500 "
                "-g {genome} "
                "--GCbiasFrequenciesFile {output.freq} "
                #Part II: Correcting GC Biases
                "correctGCBias " 
                "-b {input} "
                "--effectiveGenomeSize 2652783500 " 
                "-g {genome} "
                --GCbiasFrequenciesFile {output.freq}
                "-o {output.bam} "
                " 2> {log}"

rule deeptools:
        input:
                rules.GC_remove.output.bam
        output:
                "results/npz/{sra}_deeptool.npz"
        log:
               "results/npz/{sra}_deeptool.out" 
        conda:
               "env.yaml"
        threads:
                8
        shell:
                "multiBamSummary bins --threads {threads}"
                "-b {input} "
                "-o {output} "

rule deeptools_correlation:
        input:
                rules.deeptools.output
        output:
                png1="results/expl_data_correlation/{sra}_scatterplot_PearsonCorr_bigwigScores.png",
                tab1="results/expl_data_correlation/{sra}_scatterplot_PearsonCorr_bigwigScores.tab",
                png2="results/expl_data_correlation/{sra}_heatmap_SpearmanCorr_readCounts.png",
                tab2="results/expl_data_correlation/{sra}_heatmap_SpearmanCorr_readCounts.tab"
        params:
                title1="Pearson Correlation of Average Scores Per Transcript",
                title2="Spearman Correlation of Read Counts"
        threads:
                8
        conda:
               "env.yaml"        
        shell:
                "plotCorrelation --threads {threads} "
                "-in {input} "
                "--corMethod pearson --skipZeros "
                "--plotTitle {title1} "
                "--whatToPlot scatterplot "
                "-o {png1} "
                "--outFileCorMatrix {tab1} "
                #Plot Spearman
                "plotCorrelation "
                "-in {input} "
                "--corMethod spearman --skipZeros "
                "--plotTitle {title2} "
                "--whatToPlot heatmap --colorMap RdYlBu --plotNumbers "
                "-o {png2} "
                "--outFileCorMatrix {tab2}"