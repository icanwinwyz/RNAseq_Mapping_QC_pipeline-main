#!/usr/bin/perl -w
use strict;

for(my $i=0;$i<5;$i++){
	print "qsub -q all.q -l h=csclprd1-c",$i,"v -cwd ./remove_ref_clean.sh Human_mRNA Human_mRNA_node",$i,"\n";
	print "qsub -q all.q -l h=csclprd1-c",$i,"v -cwd ./remove_ref_clean.sh Human_totalRNA Human_total_node",$i,"\n";
	print "qsub -q all.q -l h=csclprd1-c",$i,"v -cwd ./remove_ref_clean.sh Mouse_mRNA Mouse_mRNA_node",$i,"\n";
	print "qsub -q all.q -l h=csclprd1-c",$i,"v -cwd ./remove_ref_clean.sh Mouse_totalRNA Mouse_totalRNA_node",$i,"\n";
	print "qsub -q all.q -l h=csclprd1-c",$i,"v -cwd ./remove_ref_clean.sh Rat Rat_node",$i,"\n";
	print "qsub -q all.q -l h=csclprd1-c",$i,"v -cwd ./remove_ref_clean.sh Breunig Breunig_node",$i,"\n";
#	print "qsub -q all.q -l h=csclprd1-c",$i,"v -cwd ./remove_ref_clean.sh other node",$i,"\n";
}


for(my $i=11;$i<36;$i++){
	print "qsub -q all.q -l h=csclprd1-c",$i,"v -cwd ./remove_ref_clean.sh Human_mRNA Human_mRNA_node",$i,"\n";
	print "qsub -q all.q -l h=csclprd1-c",$i,"v -cwd ./remove_ref_clean.sh Human_totalRNA Human_totalRNA_node",$i,"\n";
	print "qsub -q all.q -l h=csclprd1-c",$i,"v -cwd ./remove_ref_clean.sh Mouse_mRNA Mouse_mRNA_node",$i,"\n";
	print "qsub -q all.q -l h=csclprd1-c",$i,"v -cwd ./remove_ref_clean.sh Mouse_totalRNA Mouse_totalRNA_node",$i,"\n";
	print "qsub -q all.q -l h=csclprd1-c",$i,"v -cwd ./remove_ref_clean.sh Rat Rat_node",$i,"\n";
	print "qsub -q all.q -l h=csclprd1-c",$i,"v -cwd ./remove_ref_clean.sh Breunig Breunig_node",$i,"\n";
#	print "qsub -q all.q -l h=csclprd1-c",$i,"v -cwd ./remove_ref_clean.sh other node",$i,"\n";
}
