# RNAseq_pipeline
This is latest version of RNA-seq pipeline used in Genomics Core of Cedars Sinai Medical Center. 

# Mapping_QC_Auto_v3
Mapping/QC pipeline for RNA-seq data version 3.4


# DESCRIPTION

This pipeline integrats the Mapping, gene counts/tpm by RSEM and RseQC

This pipeline is compatible for reads of "single-end" and "paired-end" which is specified by the option "-t".

# USAGE

In the folder with only "fastq.gz" files:

nohup perl /common/genomics-core/apps/mapping_qc_auto/Mapping_QC_Auto_v3.4.pl -t <SE|PE> -o <Human_mRNA|Mouse_mRNA|Human_totalRNA|Mouse_totalRNA|Rat> -p <project_ID> -n <1,2,3,...> -qc -gb > projectid.log.txt >2&1&

example:
nohup perl /common/genomics-core/apps/mapping_qc_auto/Mapping_QC_Auto_v3.4.pl -t SE -o Mouse_mRNA -p KF-9131--03--04--2020 -n 20,21,22 -qc > KF-9131--03--04--2020.log.txt 2>&1&


# REQUIREMENT

- Perl 5

- perl module: Getopt::Long

# OPTIONS

Running options:

-t or --type [required if no samplesheet supplied]
        The sequencing type is "single end (SE)" or "paired end (PE)"

-o or --organism
        the reference genoem: Human or Mouse

-p or --project
        project ID for this run

-qc or --qualitycontrol
                whether run quality contorl (RSeQC) or not after mapping

-gb or --genebody
                whether run genebody test to check 3' or 5' bias

-n or --nodes
                which nodes the jobs will be ran on the all.q.The number of jobs shouldn't smaller than the number of nodes used

-h or --help
        Help information


# OUTPUT FILES

bam: folder saving all of bam/bai file after mapping

fastq: folder saving all of .fastq files

final_results: stat for mapping and counts and tpm matrixs.

genes_isoforms_results: folder saving *.genes.results and *.isoforms.results for each sample

log: folder saving all of log files during mapping

node_log: folder saving all or *.e* and *.o* files from clusters

others: useless files

RseQC_results: folder saving all of QC results for each sample (not generated if -qc deactive)
