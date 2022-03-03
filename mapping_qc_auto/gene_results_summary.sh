#########################################################################
# File Name: gene_results_summary.sh
# Author: Yizhou Wang
# mail: yizhou@cshs.org
# Created Time: Tue 06 Sep 2016 10:05:30 AM PDT
#########################################################################
#!/bin/bash

export PATH=/home/wangyiz/genomics/anaconda2/bin:$PATH
export PATH=/common/genomics-core/anaconda2/bin:$PATH
#/common/genomics-core/anaconda2/bin/rsem-generate-data-matrix *genes.results > $1.COUNTS.txt

#/common/genomics-core/anaconda2/bin/rsem-generate-data-tpm *genes.results > $1.TPM.txt

#/common/genomics-core/anaconda2/bin/rsem-generate-data-fpkm *genes.results > $1.FPKM.txt

#/common/genomics-core/anaconda2/bin/Rscript /common/genomics-core/apps/sequencing_data_distri/format.R $1.TPM.txt $1.COUNTS.txt $1.FPKM.txt $1


/common/genomics-core/anaconda2/bin/rsem-generate-data-matrix *genes.results > $1.COUNTS.txt

/common/genomics-core/anaconda2/bin/rsem-generate-data-tpm *genes.results > $1.TPM.txt

/common/genomics-core/anaconda2/bin/rsem-generate-data-fpkm *genes.results > $1.FPKM.txt

/hpc/apps/R/3.4.1/bin/Rscript /common/genomics-core/apps/sequencing_data_distri/format.R $1.TPM.txt $1.COUNTS.txt $1.FPKM.txt $1


rm *.COUNTS.txt
rm *.TPM.txt
rm *.FPKM.txt
