#########################################################################
# File Name: orgnize.sh
# Author: Yizhou Wang
# mail: yizhou@cshs.org
# Created Time: Tue 06 Sep 2016 06:04:53 PM PDT
#########################################################################
#!/bin/bash

mv *_fastqc.* fastqc
mv *.fastq* fastq
mv *.csv *QC.txt *QC.pdf *.geneBodyCoverage.*.pdf *.multiqc.html final_results/
mv *.genes.results genes_isoforms_results
mv *.isoforms.results genes_isoforms_results
mv *.bam *.bai bam
mv *Log* *.time log
mv *.tab others
mv -f *_STARtmp *.stat *.r *.geneBodyCoverage.txt input*_*.txt others
rm log.txt sortedBam.txt Aligned.out.sam
mv *.e* *.o* node_log

echo "Subject: Mapping/QC is done" | sendmail -v yizhou.wang@cshs.org
echo "Subject: Mapping/QC is done" | sendmail -v di.wu@cshs.org

