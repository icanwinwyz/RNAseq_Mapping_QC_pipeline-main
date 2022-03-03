export PATH=/home/wangyiz/genomics/anaconda2/bin:$PATH
export PATH=/common/genomics-core/anaconda2/bin:$PATH

module load samtools 
module load R


multiqc -n $1.multiqc -f -c /common/genomics-core/apps/multiqc/multiqc_config_RNAseq_HPC.yaml -e rsem ./  

