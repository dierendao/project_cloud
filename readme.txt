# Title
ATACSeq data set from "STATegra, a comprehensive multi-omics dataset of B-cell differentiation in mouse"
David Gomez-Cabrero et al. 2019
https://doi.org/10.1038/s41597-019-0202-7

# Where
Datasets collected from https://www.ebi.ac.uk/ena
# When
2021-07-26
# Who
Nadia Goué
#Modified by
Babacar Ndao
# How
with enaDataGet tool https://www.ebi.ac.uk/about/news/service-news/new-tools-download-data-ena
included in inhouse script ena_array_atac.slurm

# What
## Raw dataset:
SRR4785152  50k-Rep1-0h-sample.0h   GSM2367179  0.7G
SRR4785153  50k-Rep2-0h-sample.0h   GSM2367180  0.7G
SRR4785154  50k-Rep3-0h-sample.0h   GSM2367181  0.7G

SRR4785341  50k-24h-R1-sample.24h.2 GSM2367368  0.6G
SRR4785342  50k-24h-R2-sample.24h.2 GSM2367369  0.7G
SRR4785343  50k-24h-R3-sample.24h.2 GSM2367370  0.6G

## Experimental Design

Cells collected for 0 and 24hours post-treatment with tamoxifen
3 biological replicates of ~50,000 cells

## Subset
Create data subset to test tools
zcat dataset | head -n 4000000 | gzip > ss_dataset.fastq.gz
Exemple:
zcat 50k_0h_R1_1.fastq.gz | head -n 4000000 | gzip > ss_50k_0h_R1_1.fastq.gz


## Sequencing Tech

### Reference
Gomez-Cabrero et al. refer to Buenrostro et al. "ATAC-seq: A Method for Assaying Chromatin Accessibility Genome-Wide"
doi:10.1002/0471142727.mb2129s109
Illumina with Nextera-based sequencing primers

### Collect adapters for read cleaning
https://github.com/timflutre/trimmomatic/tree/master/adapters/NexteraPE-PE.fa
2015-03-5
v0.33

>PrefixNX/1
AGATGTGTATAAGAGACAG
>PrefixNX/2
AGATGTGTATAAGAGACAG
>Trans1
TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
>Trans1_rc
CTGTCTCTTATACACATCTGACGCTGCCGACGA
>Trans2
GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG
>Trans2_rc
CTGTCTCTTATACACATCTCCGAGCCCACGAGAC