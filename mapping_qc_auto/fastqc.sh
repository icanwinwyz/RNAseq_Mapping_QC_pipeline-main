export PATH=/home/wangyiz/genomics/anaconda2/bin:$PATH
export PATH=/common/genomics-core/anaconda2/bin:$PATH

module load samtools 
module load R

FASTQ=".fastq.gz"
INFILE=$1
echo $INFILE

fastqc $INFILE$FASTQ
