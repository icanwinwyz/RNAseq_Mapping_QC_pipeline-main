#########################################################################
# File Name: geneBody.sh
# Author: Yizhou Wang
# mail: yizhou@cshs.org
# Created Time: Tue 06 Sep 2016 11:07:40 AM PDT
#########################################################################
#!/bin/bash
export PATH=/home/wangyiz/genomics/anaconda2/bin:$PATH
export PATH=/common/genomics-core/anaconda2/bin:$PATH
if [ "$1" = "Human_mRNA" ]; then

	HOUSEDEEP="/common/genomics-core/reference/RseQC/GRCh38_v23/hg38.HouseKeepingGenes.bed"

if [ "$1" = "Human_totalRNA" ]; then

	HOUSEDEEP="/common/genomics-core/reference/RseQC/GRCh38_v23/hg38.HouseKeepingGenes.bed"

elif [ "$1" = "Mouse_mRNA" ]; then

	HOUSEDEEP="/common/genomics-core/reference/RseQC/GRCm38_m8/mm10.HouseKeepingGenes.bed"

elif [ "$1" = "Mouse_totalRNA" ]; then

	HOUSEDEEP="/common/genomics-core/reference/RseQC/GRCm38_m8/mm10.HouseKeepingGenes.bed"

elif [ "$1" = "Rat" ]; then

	HOUSEDEEP="/common/genomics-core/reference/RseQC/Rnor6.0/housekeeping.bed"

elif [ "$1" = "other" ]; then

	HOUSEDEEP="/common/genomics-core/reference/RseQC/GRCm38_m8/mm10.HouseKeepingGenes.bed"

elif [ "$1" = "Breunig" ]; then

	HOUSEDEEP="/common/genomics-core/reference/RseQC/GRCm38_m8/mm10.HouseKeepingGenes.bed"

else
       	echo "Please perfomr alignment manually since no reference genome available"

fi

ls *Aligned.sortedByCoord.out.bam > sortedBam.txt

geneBody_coverage.py -r $HOUSEDEEP -i sortedBam.txt -o $2


