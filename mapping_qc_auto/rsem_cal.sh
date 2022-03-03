#########################################################################
# File Name: rsem_cal.sh
# Author: Yizhou Wang
# mail: yizhou@cshs.org
# Created Time: Fri 02 Sep 2016 02:05:25 PM PDT
#########################################################################
#!/bin/bash

export PATH=/home/wangyiz/genomics/anaconda2/bin:$PATH
export PATH=/common/genomics-core/anaconda2/bin:$PATH
module load samtools
BAM="Aligned.toTranscriptome.out.bam"
COR="Aligned.sortedByCoord.out.bam"
INFILE=$3
echo $INFILE

if [ "$1" = "Human_mRNA" ]; then

		GENOMEPATH="/common/genomics-core/reference/STAR/primary_GRCh38_23"
		REFPATH="/common/genomics-core/reference/RseQC/GRCh38_v23/QC_reference_files.txt"
		RSEMPATH="/common/genomics-core/reference/STAR/primary_GRCh38_23/GRCh38_primary"

elif [ "$1" = "Mouse_mRNA" ]; then

		GENOMEPATH="/common/genomics-core/reference/STAR/primary_GRCm38_mm8"
		REFPATH="/common/genomics-core/reference/RseQC/GRCm38_m8/QC_reference_files.txt"
		RSEMPATH="/common/genomics-core/reference/STAR/primary_GRCm38_mm8/GRCm38_m8"

elif [ "$1" = "Human_totalRNA" ]; then

                GENOMEPATH="/common/genomics-core/reference/STAR/GRCh38_23_totalRNA"
                REFPATH="/common/genomics-core/reference/RseQC/GRCh38_v23/QC_reference_files.txt"
                RSEMPATH="/common/genomics-core/reference/STAR/GRCh38_23_totalRNA/GRCh38_23_totalRNA_ERCC92"

elif [ "$1" = "Mouse_totalRNA" ]; then

                GENOMEPATH="/common/genomics-core/reference/STAR/GRCm38_M8_totalRNA"
                REFPATH="/common/genomics-core/reference/RseQC/GRCm38_m8/QC_reference_files.txt"
                RSEMPATH="/common/genomics-core/reference/STAR/GRCm38_M8_totalRNA/GRCm38_M8_totalRNA"

elif [ "$1" = "Rat" ]; then

	       	GENOMEPATH="/common/genomics-core/reference/STAR/Rnor_6.0_Rat"
		REFPATH="/common/genomics-core/reference/RseQC/Rnor6.0/QC_reference_files.txt"
		RSEMPATH="/common/genomics-core/reference/STAR/Rnor_6.0_Rat/Rnor6"

elif [ "$1" = "Breunig" ]; then

GENOMEPATH="/common/genomics-core/reference/STAR/Breunig_mouse_custom"
        REFPATH="/common/genomics-core/reference/RseQC/GRCm38_m8/QC_reference_files.txt"
	RSEMPATH="/common/genomics-core/reference/STAR/Breunig_mouse_custom/Breunig_mouse"

elif [ "$1" = "Human_mRNA_hg19" ]; then

                GENOMEPATH="/common/genomics-core/reference/STAR/GRCh37"
                REFPATH="/common/genomics-core/reference/RseQC/GRCh37/QC_reference_files.txt"
		RSEMPATH="/common/genomics-core/reference/STAR/GRCh37/GRCh37"
elif [ "$1" = "test" ]; then
	GENOMEPATH="/common/genomics-core/reference/Pipeline_test/GRCh38/STAR/"
RSEMPATH="/common/genomics-core/reference/Pipeline_test/GRCh38/STAR/chr19"


else
       	echo "Please perfomr alignment manually since no reference genome available"

fi

if [ "$2" = "SE" ]; then
	echo $INFILE

	rsem-calculate-expression --no-bam-output --append-names --bam -p 5 --time $INFILE$BAM $RSEMPATH $INFILE

	samtools index $INFILE$COR

	echo "Done for RSEM calculation!"

elif [ "$2" = "PE" ]; then
	echo $INFILE

	rsem-calculate-expression --no-bam-output --paired-end --append-names --bam -p 5 --time $INFILE$BAM $RSEMPATH $INFILE

	samtools index $INFILE$COR

	echo "Done for RSEM calculation!"

fi


