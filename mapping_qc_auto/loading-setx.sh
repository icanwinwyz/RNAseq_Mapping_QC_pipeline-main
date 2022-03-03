#########################################################################
# File Name: loading.sh
# Author: Yizhou Wang
# mail: yizhou@cshs.org
# Created Time: Thu 01 Sep 2016 11:40:43 AM PDT
#########################################################################
#!/bin/bash

#set -x

export PATH=/home/wangyiz/genomics/anaconda2/bin:$PATH
export PATH=/common/genomics-core/anaconda2/bin:$PATH

if [ "$1" = "Human_mRNA" ]; then

		GENOMEPATH="/common/genomics-core/reference/STAR/primary_GRCh38_23"
		REFPATH="/common/genomics-core/reference/RseQC/GRCh38_v23/QC_reference_files.txt"

elif [ "$1" = "Mouse_mRNA" ]; then

	       	GENOMEPATH="/common/genomics-core/reference/STAR/primary_GRCm38_mm8"
		REFPATH="/common/genomics-core/reference/RseQC/GRCm38_m8/QC_reference_files.txt"

elif [ "$1" = "Human_totalRNA" ]; then

                GENOMEPATH="/common/genomics-core/reference/STAR/GRCh38_23_totalRNA"
                REFPATH="/common/genomics-core/reference/RseQC/GRCh38_v23/QC_reference_files.txt"

elif [ "$1" = "Mouse_totalRNA" ]; then

                GENOMEPATH="/common/genomics-core/reference/STAR/GRCm38_M8_totalRNA"
                REFPATH="/common/genomics-core/reference/RseQC/GRCm38_m8/QC_reference_files.txt"


elif [ "$1" = "Rat" ]; then

       		GENOMEPATH="/common/genomics-core/reference/STAR/Rnor_6.0_Rat"
		REFPATH="/common/genomics-core/reference/RseQC/Rnor6.0/QC_reference_files.txt"

elif [ "$1" = "test" ]; then
		GENOMEPATH="/common/genomics-core/reference/Pipeline_test/GRCh38/STAR/"
       		REFPATH="/common/genomics-core/reference/Pipeline_test/GRCh38/RSEQC/QC_reference_files.txt"

elif [ "$1" = "Breunig" ]; then

        	GENOMEPATH="/common/genomics-core/reference/STAR/Breunig_mouse_custom"
        	REFPATH="/common/genomics-core/reference/RseQC/GRCm38_m8/QC_reference_files.txt"

elif [ "$1" = "Human_mRNA_hg19" ]; then

                GENOMEPATH="/common/genomics-core/reference/STAR/GRCh37"
                REFPATH="/common/genomics-core/reference/RseQC/GRCh37/QC_reference_files.txt"

else
       	echo "Please perfomr alignment manually since no reference genome available"

fi


/common/genomics-core/anaconda2/bin/STAR --genomeDir $GENOMEPATH --genomeLoad LoadAndExit --outFileNamePrefix $2

echo $GENOMEPATH
echo "Done for ref genome loading!"

#s=0
#for (( i=1; i<=1000000; i=i+1 ))
#do
#	s=$(($s+$i))
#done
#echo $s
#echo 'genome loading is successful run!!'


