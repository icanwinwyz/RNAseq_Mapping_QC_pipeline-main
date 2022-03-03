#########################################################################
# File Name: rseqc_summary.sh
# Author: Yizhou Wang
# mail: yizhou@cshs.org
# Created Time: Tue 06 Sep 2016 05:32:19 PM PDT
#########################################################################
#!/bin/bash

export PATH=/common/genomics-core/anaconda2/bin:$PATH

module load R/3.4.1

Rscript /common/genomics-core/bin/Rse_QC_summary.R ./RseQC_results/*/*.final_results_summary.txt $1



