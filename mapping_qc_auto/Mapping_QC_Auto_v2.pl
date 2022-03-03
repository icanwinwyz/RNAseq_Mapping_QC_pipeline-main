#!/usr/bin/perl -w
use strict;
use Data::Dumper qw(Dumper);
use Getopt::Long;
use Term::ANSIColor   qw(:constants);
$Term::ANSIColor::AUTORESET=1; 

=head1 NAME

RNA-seq Mapping(STAR) QC RSEM pipeline v1.1

=head1 DESCRIPTION

This pipeline integrats the Mapping, gene counts/tpm by RSEM and RseQC

This pipeline is compatible for reads of "single-end" and "paired-end" which is specified by the option "-t".

=head1 USAGE

If the "SampleSheet(inputfile.txt)" supplied:

Mapping_auto.pl -e <your_email_address> 

If no "SampleSheet(inputfile.txt)" supplied:

Mapping_auto.pl -e <your_email_address> -t <SE|PE> -o <Human|Mouse> -p <project_ID> -n <1238> -qc

=head1 REQUIREMENT

- Perl 5

- perl module: Getopt::Long

=head1 OPTIONS

Running options:

-e or --email [optional]
       Provide your email address if you would like to be notified after jobs completed.	

-t or --type [required if no samplesheet supplied]
        The sequencing type is "single end (SE)" or "paired end (PE)"

-o or --organism
       	the reference genoem: Human or Mouse 

-p or --project
       	project ID for this run 

-qc or --qualitycontrol
		whether run quality contorl (RSeQC) or not after mapping

-n or --nodes
		which nodes the jobs will be ran on (not recommand using csclp1-0-0.local).The number of jobs shouldn't smaller than the number of nodes used

-h or --help
       	Help information


=head1 OUTPUT FILES

bam: folder saving all of bam/bai file after mapping

fastq: folder saving all of .fastq files

final_results: stat for mapping and counts and tpm matrixs.

genes_isoforms_results: folder saving *.genes.results and *.isoforms.results for each sample

log: folder saving all of log files during mapping 

node_log: folder saving all or *.e* and *.o* files from clusters

others: useless files

RseQC_results: folder saving all of QC results for each sample (not generated if -qc deactive)

=head1 AUTHOR

If you have any questions, please email: yizhou.wang@cshs.org

Genomics CORE, 09/09/2016

=cut


my ($email,$type,$org,$proj,$nodes);
my $qc = 0;
my $genebody = 0;

GetOptions(
	'email|e=s' => \$email,	
	'type|t=s' => \$type,
	'organism|o=s' => \$org,
	'project|p=s' => \$proj,
	'nodes|n=s' => \$nodes,
	'qualitycontrol|qc!' => \$qc,
	'genebody|gb!' => \$genebody,
	'help|h' => sub{exec('perldoc',$0);
	exit(0);},
);

my ($proj_ID,$organism,$seq,$file,@data,$number);

if( defined $type && defined $org && defined $proj){
	if($type eq "PE"){
		system("ls *fastq.gz > input_file_fastq_tmp.txt");
		open IN, "input_file_fastq_tmp.txt" or die $!;
		my @data = <IN>;
		foreach my $line (@data){
			chomp $line;
			if ($line =~ /R1/){
				my @a = split("R1",$line);
				substr($a[0],-1) = "";
				my $new_name = join(".",$a[0],"R1","fastq.gz");
				if($new_name eq $line){
					next;
				}else{
					my $cmd = "mv $line $new_name";
					system($cmd);
				}
			}elsif($line =~ /R2/){
				my @a = split("R2",$line);
				substr($a[0],-1) = "";
				my $new_name = join(".",$a[0],"R2","fastq.gz");
				if($new_name eq $line){
					next;
				}else{
					my $cmd = "mv $line $new_name";
					system($cmd);
				}
				#	}else{
				#	die "please provide paired-end reads with right format when using the "PE" argument";
			}else{
				die "please provide paired-end reads!\n";
			}
		}

#		system("rm input_file_fastq_tmp.txt");
		system("ls *.fastq.gz|sed 's/.R[1\|2].*fastq.gz//g'|sort -u > input_fastq.txt");
	}elsif($type eq "SE"){
		#	system("ls *.fastq|sed 's/.[fastq|fq]//g' > input_fastq.txt");
		my $temp_cmd = "ls *.fastq.gz|sed 's/\\\(.fastq.gz\\|.fq.gz\\)//g' > input_fastq.txt";
		system($temp_cmd);
		#	system("ls *.fastq|sed 's/\(.fastq\|.fq\)//g' > input_fastq.txt");
	}
	$proj_ID = $proj;
	$organism = $org;
	$seq = $type;
	$file = "input_fastq.txt";
	open IN, $file;
	@data = <IN>;
	$number = scalar(@data);
	if($type eq "PE"){
		foreach my $line (@data){
			chomp $line;
			my $pair1 = $line.".R1".".fastq";
			my $pair2 = $line.".R2".".fastq";
#			print $pair1,"\t",$pair2,"\n";

			if (-e $pair1 && -e $pair2){
				next;
			}else{
				die "no pair mates found for $line!\n";
			}
		}
	}
	close IN;
}elsif( -e "inputfile.txt"){
	print "inputfile exists\n";
	$file = "inputfile.txt";
	open IN, $file or die $!;
	@data = <IN>;
	my $proj_ID_cmd = "awk \'{print \$2}\' $file|uniq";
	chomp($proj_ID = `$proj_ID_cmd`);
	$number = scalar(@data);
	my $comman_org = "awk \'\{print \$3\}\' inputfile\.txt\|uniq";
	chomp($organism = `$comman_org`);
	my $comman_seq = "awk \'\{print \$4\}\' inputfile\.txt\|uniq";
	chomp($seq = `$comman_seq`);
	#print $seq,"\n";
}else{
	die "FATAL ERROR: parameters are not defined or inputfile doesn't exist! See usage with -h\n";
}
print "The sequencing type is:  ";
print RED "$seq\n";
print "The organism is: ";
print RED "$organism\n";
print "The project name is: ";
print RED "$proj_ID\n";
if (defined $email){
	print "The complete notification will be sent to: ";
	print RED "$email\n";
}else{
	print GREEN "no email notification after job complete.\n";
}

my $cmd_path="/common/genomics-core/apps/sequencing_data_distri";

if ($qc ==1){
	print  "QC for samples?  ";
	print RED "YES!\n";
}else{
	print "QC for samples?  ";
	print GREEN "NO!\n";
}

my @nodes_free=split(",",$nodes);
my @sort;
foreach my $line(@nodes_free){
	chomp $line;
	my $a = "csclprd1-c".$line."v";
	push(@sort,$a);
}
print "job will run on nodes:\n";
print RED join("\t",@sort),"\n\n\n\n";

my $nodes_number = scalar(@sort);
#print $nodes_number,"\n";

my $number_check=int($number/2);

if ($nodes_number == 1){
	print "ref genome of $organism is loading to $sort[0]\n";
	####reference genome loading####
#	print "qsub -q all.q -l h=$sort[0] -N $sort[0] -cwd $cmd_path/loading.sh $organism\n";
	system("qsub -q all.q -l h=$sort[0] -N genome_load_$nodes_free[0] -cwd $cmd_path/loading.sh $organism");
	####Mapping by STAR####
	if($number == 1){
	my $mapping_comman1 = "qsub -q all.q -l h=$sort[0] -hold_jid genome_load_$nodes_free[0] -N mapping_$nodes_free[0] -t 1 -tc 1 -cwd $cmd_path/mapping_star.sh $organism $seq $file first";
	system($mapping_comman1);
	
#	my $mapping_check = "qsub -q all.q -l h=$sort[0] -hold_jid mapping_$nodes_free[0] -N check_mapping_$nodes_free[0] -cwd $cmd_path/check_mapping.sh";
#	system($mapping_check);
	
#	my $mapping_comman2 = "qsub -q all.q -l h=$sort[0] -hold_jid check_mapping_$nodes_free[0] -N mapping_retry_$nodes_free[0] -t 1 -tc 1 -cwd $cmd_path/mapping_star.sh $organism $seq mapping_retry\.txt second";
 #       system($mapping_comman2);
####Remove reference genome on nodes####
	my $remove_comman = "qsub -q all.q -l h=$sort[0] -hold_jid mapping_$nodes_free[0] -N remove_genome_$nodes_free[0] -cwd $cmd_path/remove_ref.sh $organism";
	system($remove_comman);
	my $resem_comman1 = "qsub -q all.q -l h=$sort[0] -hold_jid remove_genome_$nodes_free[0] -N rsem_$nodes_free[0] -t 1 -tc 1 -cwd $cmd_path/rsem_cal.sh $organism $seq $file first";
		system($resem_comman1);
		
#		my $rsem_check = "qsub -q all.q -l h=$sort[0] -hold_jid rsem_$nodes_free[0] -N check_rsem_$nodes_free[0] -cwd $cmd_path/check_rsem.sh";
#	system($rsem_check);
	
#		my $rsem_comman2 = "qsub -q all.q -l h=$sort[0] -hold_jid check_rsem_$nodes_free[0] -N rsem_retry_$nodes_free[0] -t 1 -tc 1 -cwd $cmd_path/rsem_cal.sh $organism $seq rsem_retry\.txt second";
 #       system($rsem_comman2);
####RseQC####
		my $rseqc_comman1 = "qsub -q all.q -l h=$sort[0] -hold_jid remove_genome_$nodes_free[0] -N rseqc_$nodes_free[0] -t 1 -tc 1 -cwd $cmd_path/rseqc.sh $organism $seq $file first";
		system($rseqc_comman1);
#		my $rseqc_check = "qsub -q all.q -l h=$sort[0] -hold_jid rseqc_$nodes_free[0] -N check_rseqc_$nodes_free[0] -cwd $cmd_path/check_rseqc.sh";
 #      		system($rseqc_check);
 #               my $rseqc_comman2 = "qsub -q all.q -l h=$sort[0] -hold_jid check_rseqc_$nodes_free[0] -N rseqc_retry_$nodes_free[0] -t 1 -tc 1 -cwd $cmd_path/rseqc.sh $organism $seq rseqc_retry\.txt second";
 #       	system($rseqc_comman2);
####RSEM_summary####
		my $rsem_summary = "qsub -q all.q -l h=$sort[0] -hold_jid rsem_$nodes_free[0] -N rsem_count_tmp_summary_$nodes_free[0] -cwd $cmd_path/gene_results_summary.sh $proj_ID";
		system($rsem_summary);
####Gene_body_coverage####
		my $rseqc_summary = "qsub -q all.q -hold_jid rseqc_$nodes_free[0] -N rseqc_summary_$nodes_free[0] -cwd $cmd_path/rseqc_summary.sh";
		system($rseqc_summary);
		system("qsub -q all.q -hold_jid rseqc_summary_$nodes_free[0],rsem_count_tmp_summary_$nodes_free[0] -N folder_build_$nodes_free[0] -cwd $cmd_path/mkdir.sh");
			my $cmd_email="qsub -q all.q -hold_jid folder_build_$nodes_free[0] -o ./node_log -e ./node_log -N organize_$nodes_free[0] -cwd $cmd_path/organize.sh";
			system($cmd_email);	
	}else{
	my $mapping_comman1 = "qsub -q all.q -l h=$sort[0] -hold_jid genome_load_$nodes_free[0] -N mapping_$nodes_free[0] -t 1\-$number -tc 2 -cwd $cmd_path/mapping_star.sh $organism $seq $file first";
	system($mapping_comman1);
	
	my $mapping_check = "qsub -q all.q -l h=$sort[0] -hold_jid mapping_$nodes_free[0] -N check_mapping_$nodes_free[0] -cwd $cmd_path/check_mapping.sh";
	system($mapping_check);
	
#	my $number_check=int($number/2);

	my $mapping_comman2 = "qsub -q all.q -l h=$sort[0] -hold_jid check_mapping_$nodes_free[0] -N mapping_retry_$nodes_free[0] -t 1\-$number_check -tc 1 -cwd $cmd_path/mapping_star.sh $organism $seq mapping_retry\.txt second";
        system($mapping_comman2);
####Remove reference genome on nodes####
	my $remove_comman = "qsub -q all.q -l h=$sort[0] -hold_jid mapping_retry_$nodes_free[0] -N remove_genome_$nodes_free[0] -cwd $cmd_path/remove_ref.sh $organism";
	system($remove_comman);
	
	if ($qc == 1){
####Rsem calculation####
		my $resem_comman1 = "qsub -q all.q -l h=$sort[0] -hold_jid remove_genome_$nodes_free[0] -N rsem_$nodes_free[0] -t 1\-$number -tc 2 -cwd $cmd_path/rsem_cal.sh $organism $seq $file first";
		system($resem_comman1);
		
		my $rsem_check = "qsub -q all.q -l h=$sort[0] -hold_jid rsem_$nodes_free[0] -N check_rsem_$nodes_free[0] -cwd $cmd_path/check_rsem.sh";
	system($rsem_check);
	
		my $rsem_comman2 = "qsub -q all.q -l h=$sort[0] -hold_jid check_rsem_$nodes_free[0] -N rsem_retry_$nodes_free[0] -t 1\-$number_check -tc 1 -cwd $cmd_path/rsem_cal.sh $organism $seq rsem_retry\.txt second";
        system($rsem_comman2);
####RseQC####
		my $rseqc_comman1 = "qsub -q all.q -l h=$sort[0] -hold_jid remove_genome_$nodes_free[0] -N rseqc_$nodes_free[0] -t 1\-$number -tc 2 -cwd $cmd_path/rseqc.sh $organism $seq $file first";
		system($rseqc_comman1);
		my $rseqc_check = "qsub -q all.q -l h=$sort[0] -hold_jid rseqc_$nodes_free[0] -N check_rseqc_$nodes_free[0] -cwd $cmd_path/check_rseqc.sh";
       		system($rseqc_check);
                my $rseqc_comman2 = "qsub -q all.q -l h=$sort[0] -hold_jid check_rseqc_$nodes_free[0] -N rseqc_retry_$nodes_free[0] -t 1\-$number_check -tc 1 -cwd $cmd_path/rseqc.sh $organism $seq rseqc_retry\.txt second";
        	system($rseqc_comman2);
####RSEM_summary####
		my $rsem_summary = "qsub -q all.q -l h=$sort[0] -hold_jid rsem_retry_$nodes_free[0] -N rsem_count_tmp_summary_$nodes_free[0] -cwd $cmd_path/gene_results_summary.sh $proj_ID";
		system($rsem_summary);
####Gene_body_coverage####
		my $rseqc_summary = "qsub -q all.q -hold_jid rseqc_retry_$nodes_free[0] -N rseqc_summary_$nodes_free[0] -cwd $cmd_path/rseqc_summary.sh";
		system($rseqc_summary);
		if ($genebody == 1){
			my $gene_body = "qsub -q all.q -l h=$sort[0] -hold_jid rsem_retry_$nodes_free[0],rseqc_retry_$nodes_free[0] -N geneBody_$nodes_free[0] -cwd $cmd_path/geneBody.sh $organism $proj_ID";
			system($gene_body);
			system("qsub -q all.q -hold_jid rseqc_summary_$nodes_free[0],geneBody_$nodes_free[0] -N folder_build_$nodes_free[0] -cwd $cmd_path/mkdir.sh");
			my $cmd_email="qsub -q all.q -hold_jid folder_build_$nodes_free[0] -o ./node_log -e ./node_log -N organize_$nodes_free[0] -cwd $cmd_path/organize.sh";
			system($cmd_email);	
		}else{
			system("qsub -q all.q -hold_jid rseqc_summary_$nodes_free[0],rsem_count_tmp_summary_$nodes_free[0] -N folder_build_$nodes_free[0] -cwd $cmd_path/mkdir.sh");
			my $cmd_email="qsub -q all.q -hold_jid folder_build_$nodes_free[0] -o ./node_log -e ./node_log -N organize_$nodes_free[0] -cwd $cmd_path/organize.sh";
			system($cmd_email);	
		}
			
	}else{
		my $resem_comman1 = "qsub -q all.q -l h=$sort[0] -hold_jid remove_genome_$nodes_free[0] -N rsem_$nodes_free[0] -t 1\-$number -tc 4 -cwd $cmd_path/rsem_cal.sh $organism $seq $file first";
		system($resem_comman1);
		
		my $rsem_check = "qsub -q all.q -l h=$sort[0] -hold_jid rsem_$nodes_free[0] -N check_rsem_$nodes_free[0] -cwd $cmd_path/check_rsem.sh";
	system($rsem_check);
	
		my $rsem_comman2 = "qsub -q all.q -l h=$sort[0] -hold_jid check_rsem_$nodes_free[0] -N rsem_retry_$nodes_free[0] -t 1\-$number_check -tc 2 -cwd $cmd_path/rsem_cal.sh $organism $seq rsem_retry\.txt second";
        system($rsem_comman2);
		my $rsem_summary = "qsub -q all.q -l h=$sort[0] -hold_jid rsem_retry_$nodes_free[0] -N rsem_count_tmp_summary_$nodes_free[0] -cwd $cmd_path/gene_results_summary.sh $proj_ID";
		system($rsem_summary);
		system("qsub -q all.q -hold_jid rsem_count_tmp_summary_$nodes_free[0] -N folder_build_$nodes_free[0] -cwd $cmd_path/mkdir_noqc.sh");
		if(defined $email){
			my $cmd_email="qsub -q all.q -hold_jid folder_build_$nodes_free[0] -m e -M $email -o ./node_log -e ./node_log -N organize_$nodes_free[0] -cwd $cmd_path/organize_noqc.sh";
			system($cmd_email);	
#		print $cmd_email,"\n";
		}else{
		my $cmd_email="qsub -q all.q -hold_jid folder_build_$nodes_free[0] -o ./node_log -e ./node_log -N organize_$nodes_free[0] -cwd $cmd_path/organize_noqc.sh";
		system($cmd_email);	
		}	
	}
	}
}


if ($nodes_number == 2){
	if($number < 2){
		die "USE TOO MANY NODES, PLEAE USE LESS NODES SMALL THAN 2!";
	}
	print "ref genome of $organism is loading to $sort[0] and $sort[1]\n";
	####reference genome loading####
	system("qsub -q all.q -l h=$sort[0] -N genome_load_$nodes_free[0] -cwd $cmd_path/loading.sh $organism");
	system("qsub -q all.q -l h=$sort[1] -hold_jid genome_load_$nodes_free[0] -N genome_load_$nodes_free[1] -cwd $cmd_path/loading.sh $organism");
	####Mapping by STAR####
	my $div = int($number/2);
	open OUT, ">inputfile1_temp.txt" or die $!;
	open OUTT, ">inputfile2_temp.txt" or die $!;
	print OUT @data[0..($div-1)];
	print OUTT @data[($div)..$#data];
	my $number1 = $#data-$div+1;
	my $mapping_comman1 = "qsub -q all.q -l h=$sort[0] -hold_jid genome_load_$nodes_free[0] -N mapping_$nodes_free[0] -t 1\-$div -tc 2 -cwd $cmd_path/mapping_star.sh $organism $seq inputfile1_temp\.txt first";
	system($mapping_comman1);
	my $mapping_comman2 = "qsub -q all.q -l h=$sort[1] -hold_jid genome_load_$nodes_free[1] -N mapping_$nodes_free[1] -t 1\-$number1 -tc 2 -cwd $cmd_path/mapping_star.sh $organism $seq inputfile2_temp\.txt first";
	system($mapping_comman2);
	my $mapping_check = "qsub -q all.q -l h=$sort[0] -hold_jid mapping_$nodes_free[0],mapping_$nodes_free[1] -N check_mapping_$nodes_free[0] -cwd $cmd_path/check_mapping.sh";
	system($mapping_check);
#	my $number_check=int($number/2);
	my $mapping_comman3 = "qsub -q all.q -l h=$sort[0] -hold_jid check_mapping_$nodes_free[0] -N mapping_retry_$nodes_free[0] -t 1\-$number_check -tc 1 -cwd $cmd_path/mapping_star.sh $organism $seq mapping_retry\.txt second";
	system($mapping_comman3);
	####Remove reference genome on nodes####
	my $remove_comman_1 = "qsub -q all.q -l h=$sort[0] -hold_jid mapping_retry_$nodes_free[0] -N remove_$nodes_free[0] -cwd $cmd_path/remove_ref.sh $organism";
	system($remove_comman_1);
	my $remove_comman_2 = "qsub -q all.q -l h=$sort[1] -hold_jid mapping_retry_$nodes_free[0],remove_$nodes_free[0] -N remove_$nodes_free[1] -cwd $cmd_path/remove_ref.sh $organism";
	system($remove_comman_2);
	
	my $resem_comman1 = "qsub -q all.q -l h=$sort[0] -hold_jid remove_$nodes_free[0] -N rsem_$nodes_free[0] -t 1\-$div -tc 2 -cwd $cmd_path/rsem_cal.sh $organism $seq inputfile1_temp\.txt first";
	system($resem_comman1);
	my $resem_comman2 = "qsub -q all.q -l h=$sort[1] -hold_jid remove_$nodes_free[1] -N rsem_$nodes_free[1] -t 1\-$number1 -tc 2 -cwd $cmd_path/rsem_cal.sh $organism $seq inputfile2_temp\.txt first";
	system($resem_comman2);
	my $rsem_check = "qsub -q all.q -l h=$sort[0] -hold_jid rsem_$nodes_free[0],rsem_$nodes_free[1] -N check_rsem_$nodes_free[0] -cwd $cmd_path/check_rsem.sh";
	system($rsem_check);
	my $rsem_comman3 = "qsub -q all.q -l h=$sort[0] -hold_jid check_rsem_$nodes_free[0] -N rsem_retry_$nodes_free[0] -t 1\-$number_check -tc 1 -cwd $cmd_path/rsem_cal.sh $organism $seq rsem_retry\.txt second";
	system($rsem_comman3);

	####RSEM_summary####
	my $rsem_summary = "qsub -q all.q -l h=$sort[0] -hold_jid rsem_retry_$nodes_free[0] -N rsem_count_tmp_summary -cwd $cmd_path/gene_results_summary.sh $proj_ID";
	system($rsem_summary);
	if($qc == 1){
		####RseQC####
		my $rseqc_comman1 = "qsub -q all.q -l h=$sort[0] -hold_jid rsem_$nodes_free[0] -N rseqc_$nodes_free[0] -t 1\-$div -tc 2 -cwd $cmd_path/rseqc.sh $organism $seq inputfile1_temp\.txt first";
		system($rseqc_comman1);
		my $rseqc_comman2 = "qsub -q all.q -l h=$sort[1] -hold_jid rsem_$nodes_free[1] -N rseqc_$nodes_free[1] -t 1\-$number1 -tc 2 -cwd $cmd_path/rseqc.sh $organism $seq inputfile2_temp\.txt first";
		system($rseqc_comman2);
		my $rseqc_check = "qsub -q all.q -l h=$sort[0] -hold_jid rseqc_$nodes_free[0],rseqc_$nodes_free[1] -N check_rseqc_$nodes_free[0] -cwd $cmd_path/check_rseqc.sh";
		system($rseqc_check);
		my $rseqc_comman3 = "qsub -q all.q -l h=$sort[0] -hold_jid check_rseqc_$nodes_free[0] -N rseqc_retry_$nodes_free[0] -t 1\-$number_check -tc 1 -cwd $cmd_path/rseqc.sh $organism $seq rseqc_retry\.txt second";
		system($rseqc_comman3);
		my $rseqc_summary = "qsub -q all.q -hold_jid rseqc_retry_$nodes_free[0] -N rseqc_summary_ -cwd $cmd_path/rseqc_summary.sh";
		system($rseqc_summary);
		if ($genebody == 1){
	####Gene_body_coverage####
			my $gene_body = "qsub -q all.q -l h=$sort[1] -hold_jid rsem_$nodes_free[0],rsem_$nodes_free[1] -N geneBody_ -cwd $cmd_path/geneBody.sh $organism $proj_ID";
			system($gene_body);
			system("qsub -q all.q -hold_jid rseqc_summary_,geneBody_ -N folder_build -cwd $cmd_path/mkdir.sh");
			my $cmd_email="qsub -q all.q -hold_jid folder_build -o ./node_log -e ./node_log -N organize -cwd $cmd_path/organize.sh";
			system($cmd_email);	
		}else{
			system("qsub -q all.q -hold_jid rseqc_summary_ -N folder_build -cwd $cmd_path/mkdir.sh");
			my $cmd_email="qsub -q all.q -hold_jid folder_build -o ./node_log -e ./node_log -N organize -cwd $cmd_path/organize.sh";
			system($cmd_email);	
		}
	}else{
		system("qsub -q all.q -hold_jid rsem_count_tmp_summary -N folder_build -cwd $cmd_path/mkdir_noqc.sh");
		if(defined $email){
			my $cmd_email="qsub -q all.q -hold_jid folder_build -m e -M $email -o ./node_log -e ./node_log -N organize -cwd $cmd_path/organize_noqc.sh";
			system($cmd_email);	
#		print $cmd_email,"\n";
		}else{
		my $cmd_email="qsub -q all.q -hold_jid folder_build -o ./node_log -e ./node_log -N organize -cwd $cmd_path/organize_noqc.sh";
		system($cmd_email);	
		}	
	}
}

if($nodes_number == 3){
	if($number < 3){
		die "USE TOO MANY NODES, PLEAE USE LESS NODES SMALL THAN 3!";
	}
	print "ref genome of $organism is loading to $sort[0] and $sort[1] and $sort[2]\n";
	####reference genome loading####
	system("qsub -q all.q -l h=$sort[0] -N genome_load_$nodes_free[0] -cwd $cmd_path/loading.sh $organism");
	system("qsub -q all.q -l h=$sort[1] -hold_jid genome_load_$nodes_free[0] -N genome_load_$nodes_free[1] -cwd $cmd_path/loading.sh $organism");
	system("qsub -q all.q -l h=$sort[2] -hold_jid genome_load_$nodes_free[1] -N genome_load_$nodes_free[2] -cwd $cmd_path/loading.sh $organism");
	####Mapping by STAR####
	my $div = int($number/3);
	open OUT, ">inputfile1_temp.txt" or die $!;
	open OUTT, ">inputfile2_temp.txt" or die $!;
	open OUTTT, ">inputfile3_temp.txt" or die $!;
	print OUT @data[0..($div-1)];
	print OUTT @data[$div..($div*2-1)];
	print OUTTT @data[($div*2)..$#data];
	my $number1 = $#data-$div*2+1;
	my $mapping_comman1 = "qsub -q all.q -l h=$sort[0] -hold_jid genome_load_$nodes_free[0] -N mapping_$nodes_free[0] -t 1\-$div -tc 3 -cwd $cmd_path/mapping_star.sh $organism $seq inputfile1_temp\.txt";
	system($mapping_comman1);
	my $mapping_comman2 = "qsub -q all.q -l h=$sort[1] -hold_jid genome_load_$nodes_free[1] -N mapping_$nodes_free[1] -t 1\-$div -tc 3 -cwd $cmd_path/mapping_star.sh $organism $seq inputfile2_temp\.txt";
	system($mapping_comman2);
	my $mapping_comman3 = "qsub -q all.q -l h=$sort[2] -hold_jid genome_load_$nodes_free[2] -N mapping_$nodes_free[2] -t 1\-$number1 -tc 3 -cwd $cmd_path/mapping_star.sh $organism $seq inputfile3_temp\.txt";
	system($mapping_comman3);
	my $mapping_check = "qsub -q all.q -l h=$sort[0] -hold_jid mapping_$nodes_free[0],mapping_$nodes_free[1],mapping_$nodes_free[2] -N check_mapping_$nodes_free[0] -cwd $cmd_path/check_mapping.sh";
	system($mapping_check);
#	my $number_check=int($number/2);
	my $mapping_comman4 = "qsub -q all.q -l h=$sort[0] -hold_jid check_mapping_$nodes_free[0] -N mapping_retry_$nodes_free[0] -t 1\-$number_check -tc 1 -cwd $cmd_path/mapping_star.sh $organism $seq mapping_retry\.txt second";
	system($mapping_comman4);
	####Remove reference genome on nodes####
	my $remove_comman_1 = "qsub -q all.q -l h=$sort[0] -hold_jid mapping_retry_$nodes_free[0] -N remove_$nodes_free[0] -cwd $cmd_path/remove_ref.sh $organism";
	system($remove_comman_1);
	my $remove_comman_2 = "qsub -q all.q -l h=$sort[1] -hold_jid remove_$nodes_free[0] -N remove_$nodes_free[1] -cwd $cmd_path/remove_ref.sh $organism";
	system($remove_comman_2);
	my $remove_comman_3 = "qsub -q all.q -l h=$sort[2] -hold_jid remove_$nodes_free[1] -N remove_$nodes_free[2] -cwd $cmd_path/remove_ref.sh $organism";
	system($remove_comman_3);
	####Rsem calculation####
	my $resem_comman1 = "qsub -q all.q -l h=$sort[0] -hold_jid remove_$nodes_free[0] -N rsem_$nodes_free[0] -t 1\-$div -tc 3 -cwd $cmd_path/rsem_cal.sh $organism $seq inputfile1_temp\.txt first";
	system($resem_comman1);
	my $resem_comman2 = "qsub -q all.q -l h=$sort[1] -hold_jid remove_$nodes_free[1] -N rsem_$nodes_free[1] -t 1\-$div -tc 3 -cwd $cmd_path/rsem_cal.sh $organism $seq inputfile2_temp\.txt first";
	system($resem_comman2);
	my $resem_comman3 = "qsub -q all.q -l h=$sort[2] -hold_jid remove_$nodes_free[2] -N rsem_$nodes_free[2] -t 1\-$number1 -tc 3 -cwd $cmd_path/rsem_cal.sh $organism $seq inputfile3_temp\.txt first";
	system($resem_comman3);
	my $rsem_check = "qsub -q all.q -l h=$sort[0] -hold_jid rsem_$nodes_free[0],rsem_$nodes_free[1],rsem_$nodes_free[2] -N check_rsem_$nodes_free[0] -cwd $cmd_path/check_rsem.sh";
	system($rsem_check);
	my $rsem_comman4 = "qsub -q all.q -l h=$sort[0] -hold_jid check_rsem_$nodes_free[0] -N rsem_retry_$nodes_free[0] -t 1\-$number_check -tc 1 -cwd $cmd_path/rsem_cal.sh $organism $seq rsem_retry\.txt second";
	system($rsem_comman4);
	####RSEM_summary####
	my $rsem_summary = "qsub -q all.q  -hold_jid rsem_retry_$nodes_free[0] -N rsem_count_tmp_summary -cwd $cmd_path/gene_results_summary.sh $proj_ID";
	system($rsem_summary);
	if($qc == 1){
		####RseQC####
		my $rseqc_comman1 = "qsub -q all.q -l h=$sort[0] -hold_jid rsem_$nodes_free[0] -N rseqc_$nodes_free[0] -t 1\-$div -tc 3 -cwd $cmd_path/rseqc.sh $organism $seq inputfile1_temp\.txt first";
		system($rseqc_comman1);
		my $rseqc_comman2 = "qsub -q all.q -l h=$sort[1] -hold_jid rsem_$nodes_free[1] -N rseqc_$nodes_free[1] -t 1\-$div -tc 3 -cwd $cmd_path/rseqc.sh $organism $seq inputfile2_temp\.txt first";
		system($rseqc_comman2);
		my $rseqc_comman3 = "qsub -q all.q -l h=$sort[2] -hold_jid rsem_$nodes_free[2] -N rseqc_$nodes_free[2] -t 1\-$number1 -tc 3 -cwd $cmd_path/rseqc.sh $organism $seq inputfile3_temp\.txt first";
		system($rseqc_comman3);
		my $rseqc_check = "qsub -q all.q -l h=$sort[0] -hold_jid rseqc_$nodes_free[0],rseqc_$nodes_free[1],rseqc_$nodes_free[2] -N check_rseqc_$nodes_free[0] -cwd $cmd_path/check_rseqc.sh";
		system($rseqc_check);
		my $rseqc_comman4 = "qsub -q all.q -l h=$sort[0] -hold_jid check_rseqc_$nodes_free[0] -N rseqc_retry_$nodes_free[0] -t 1\-$number_check -tc 1 -cwd $cmd_path/rseqc.sh $organism $seq rseqc_retry\.txt second";
		system($rseqc_comman4);
		my $rseqc_summary = "qsub -q all.q -hold_jid rseqc_retry_$nodes_free[0] -N rseqc_summary_ -cwd $cmd_path/rseqc_summary.sh";
		system($rseqc_summary);
		if ($genebody == 1){
	####Gene_body_coverage####
			my $gene_body = "qsub -q all.q -l h=$sort[1] -hold_jid rsem_retry_$nodes_free[0] -N geneBody_ -cwd $cmd_path/geneBody.sh $organism $proj_ID";
			system($gene_body);
			system("qsub -q all.q -hold_jid rseqc_summary_,geneBody_ -N folder_build -cwd $cmd_path/mkdir.sh");
			my $cmd_email="qsub -q all.q -hold_jid folder_build -o ./node_log -e ./node_log -N organize -cwd $cmd_path/organize.sh";
			system($cmd_email);	
		}else{
			system("qsub -q all.q -hold_jid rseqc_summary_ -N folder_build -cwd $cmd_path/mkdir.sh");
			my $cmd_email="qsub -q all.q -hold_jid folder_build -o ./node_log -e ./node_log -N organize -cwd $cmd_path/organize.sh";
			system($cmd_email);	
		}
	}else{
		system("qsub -q all.q -hold_jid rsem_count_tmp_summary -N folder_build -cwd $cmd_path/mkdir_noqc.sh");
		if(defined $email){
			my $cmd_email="qsub -q all.q -hold_jid folder_build -m e -M $email -o ./node_log -e ./node_log -N organize -cwd $cmd_path/organize_noqc.sh";
			system($cmd_email);	
#		print $cmd_email,"\n";
		}else{
		my $cmd_email="qsub -q all.q -hold_jid folder_build -o ./node_log -e ./node_log -N organize -cwd $cmd_path/organize_noqc.sh";
		system($cmd_email);	
		}	
	}
}


if($nodes_number == 4){
	if($number < 4){
		die "USE TOO MANY NODES, PLEAE USE LESS NODES SMALL THAN 4!";
	}
	####Reference genome loading####
	print "ref genome of $organism is loading to $sort[0] and $sort[1] and $sort[2] and $sort[3]\n";
	system("qsub -q all.q -l h=$sort[0] -N genome_load_$nodes_free[0] -cwd $cmd_path/loading.sh $organism");
	system("qsub -q all.q -l h=$sort[1] -hold_jid genome_load_$nodes_free[0] -N genome_load_$nodes_free[1] -cwd $cmd_path/loading.sh $organism");
	system("qsub -q all.q -l h=$sort[2] -hold_jid genome_load_$nodes_free[1] -N genome_load_$nodes_free[2] -cwd $cmd_path/loading.sh $organism");
	system("qsub -q all.q -l h=$sort[3] -hold_jid genome_load_$nodes_free[2] -N genome_load_$nodes_free[3] -cwd $cmd_path/loading.sh $organism");
	####Mapping by STAR####
	my $div = int($number/4);
	open OUT, ">inputfile1_temp.txt" or die $!;
	open OUTT, ">inputfile2_temp.txt" or die $!;
	open OUTTT, ">inputfile3_temp.txt" or die $!;
	open OUTTTT, ">inputfile4_temp.txt" or die $!;
	print OUT @data[0..($div-1)];
	print OUTT @data[$div..($div*2-1)];
	print OUTTT @data[($div*2)..($div*3-1)];
	print OUTTTT @data[($div*3)..$#data];
	my $number1 = $#data-$div*3+1;
	my $mapping_comman1 = "qsub -q all.q -l h=$sort[0] -hold_jid genome_load_$nodes_free[0] -N mapping_$nodes_free[0] -t 1\-$div -tc 4 -cwd $cmd_path/mapping_star.sh $organism $seq inputfile1_temp\.txt first";
	system($mapping_comman1);
	my $mapping_comman2 = "qsub -q all.q -l h=$sort[1] -hold_jid genome_load_$nodes_free[1] -N mapping_$nodes_free[1] -t 1\-$div -tc 4 -cwd $cmd_path/mapping_star.sh $organism $seq inputfile2_temp\.txt first";
	system($mapping_comman2);
	my $mapping_comman3 = "qsub -q all.q -l h=$sort[2] -hold_jid genome_load_$nodes_free[2] -N mapping_$nodes_free[2] -t 1\-$div -tc 4 -cwd $cmd_path/mapping_star.sh $organism $seq inputfile3_temp\.txt first";
	system($mapping_comman3);
	my $mapping_comman4 = "qsub -q all.q -l h=$sort[3] -hold_jid genome_load_$nodes_free[3] -N mapping_$nodes_free[3] -t 1\-$number1 -tc 4 -cwd $cmd_path/mapping_star.sh $organism $seq inputfile4_temp\.txt first";
	system($mapping_comman4);
	my $mapping_check = "qsub -q all.q -l h=$sort[0] -hold_jid mapping_$nodes_free[0],mapping_$nodes_free[1],mapping_$nodes_free[2],mapping_$nodes_free[3] -N check_mapping_$nodes_free[0] -cwd $cmd_path/check_mapping.sh";
	system($mapping_check);
#	my $number_check=int($number/2);
	my $mapping_comman5 = "qsub -q all.q -l h=$sort[0] -hold_jid check_mapping_$nodes_free[0] -N mapping_retry_$nodes_free[0] -t 1\-$number_check -tc 1 -cwd $cmd_path/mapping_star.sh $organism $seq mapping_retry\.txt second";
	system($mapping_comman5);
	####Remove reference genome on nodes####
	my $remove_comman_1 = "qsub -q all.q -l h=$sort[0] -hold_jid mapping_retry_$nodes_free[0] -N remove_$nodes_free[0] -cwd $cmd_path/remove_ref.sh $organism";
	system($remove_comman_1);
	my $remove_comman_2 = "qsub -q all.q -l h=$sort[1] -hold_jid remove_$nodes_free[0] -N remove_$nodes_free[1] -cwd $cmd_path/remove_ref.sh $organism";
	system($remove_comman_2);
	my $remove_comman_3 = "qsub -q all.q -l h=$sort[2] -hold_jid remove_$nodes_free[1] -N remove_$nodes_free[2] -cwd $cmd_path/remove_ref.sh $organism";
	system($remove_comman_3);
	my $remove_comman_4 = "qsub -q all.q -l h=$sort[3] -hold_jid remove_$nodes_free[2] -N remove_$nodes_free[3] -cwd $cmd_path/remove_ref.sh $organism";
	system($remove_comman_4);
	####Rsem calculation####
	my $resem_comman1 = "qsub -q all.q -l h=$sort[0] -hold_jid remove_$nodes_free[0] -N rsem_$nodes_free[0] -t 1\-$div -tc 4 -cwd $cmd_path/rsem_cal.sh $organism $seq inputfile1_temp\.txt first";
	system($resem_comman1);
	my $resem_comman2 = "qsub -q all.q -l h=$sort[1] -hold_jid remove_$nodes_free[1] -N rsem_$nodes_free[1] -t 1\-$div -tc 4 -cwd $cmd_path/rsem_cal.sh $organism $seq inputfile2_temp\.txt first";
	system($resem_comman2);
	my $resem_comman3 = "qsub -q all.q -l h=$sort[2] -hold_jid remove_$nodes_free[2] -N rsem_$nodes_free[2] -t 1\-$div -tc 4 -cwd $cmd_path/rsem_cal.sh $organism $seq inputfile3_temp\.txt first";
	system($resem_comman3);
	my $resem_comman4 = "qsub -q all.q -l h=$sort[3] -hold_jid remove_$nodes_free[3] -N rsem_$nodes_free[3] -t 1\-$number1 -tc 4 -cwd $cmd_path/rsem_cal.sh $organism $seq inputfile4_temp\.txt first";
	system($resem_comman4);
	my $rsem_check = "qsub -q all.q -l h=$sort[0] -hold_jid rsem_$nodes_free[0],rsem_$nodes_free[1],rsem_$nodes_free[2],rsem_$nodes_free[3] -N check_rsem_$nodes_free[0] -cwd $cmd_path/check_rsem.sh";
	system($rsem_check);
	my $rsem_comman5 = "qsub -q all.q -l h=$sort[0] -hold_jid check_rsem_$nodes_free[0] -N rsem_retry_$nodes_free[0] -t 1\-$number_check -tc 1 -cwd $cmd_path/rsem_cal.sh $organism $seq rsem_retry\.txt second";
	system($rsem_comman5);
	
####RSEM_summary####
	my $rsem_summary = "qsub -q all.q -hold_jid rsem_retry_$nodes_free[0] -N rsem_count_tmp_summary -cwd $cmd_path/gene_results_summary.sh $proj_ID";
	system($rsem_summary);
	if($qc == 1){
		####RseQC####
		my $rseqc_comman1 = "qsub -q all.q -l h=$sort[0] -hold_jid rsem_$nodes_free[0] -N rseqc_$nodes_free[0] -t 1\-$div -tc 4 -cwd $cmd_path/rseqc.sh $organism $seq inputfile1_temp\.txt first";
		system($rseqc_comman1);
		my $rseqc_comman2 = "qsub -q all.q -l h=$sort[1] -hold_jid rsem_$nodes_free[1] -N rseqc_$nodes_free[1] -t 1\-$div -tc 4 -cwd $cmd_path/rseqc.sh $organism $seq inputfile2_temp\.txt first";
		system($rseqc_comman2);
		my $rseqc_comman3 = "qsub -q all.q -l h=$sort[2] -hold_jid rsem_$nodes_free[2] -N rseqc_$nodes_free[2] -t 1\-$div -tc 4 -cwd $cmd_path/rseqc.sh $organism $seq inputfile3_temp\.txt first";
		system($rseqc_comman3);
		my $rseqc_comman4 = "qsub -q all.q -l h=$sort[3] -hold_jid rsem_$nodes_free[3] -N rseqc_$nodes_free[3] -t 1\-$number1 -tc 4 -cwd $cmd_path/rseqc.sh $organism $seq inputfile4_temp\.txt first";
		system($rseqc_comman4);
		my $rseqc_check = "qsub -q all.q -l h=$sort[0] -hold_jid rseqc_$nodes_free[0],rseqc_$nodes_free[1],rseqc_$nodes_free[2],rseqc_$nodes_free[3] -N check_rseqc_$nodes_free[0] -cwd $cmd_path/check_rseqc.sh";
		system($rseqc_check);
		my $rseqc_comman5 = "qsub -q all.q -l h=$sort[0] -hold_jid check_rseqc_$nodes_free[0] -N rseqc_retry_$nodes_free[0] -t 1\-$number_check -tc 1 -cwd $cmd_path/rseqc.sh $organism $seq rseqc_retry\.txt second";
		system($rseqc_comman5);
		my $rseqc_summary = "qsub -q all.q -hold_jid rseqc_retry_$nodes_free[0] -N rseqc_summary_ -cwd $cmd_path/rseqc_summary.sh";
		system($rseqc_summary);
		if ($genebody == 1){
	####Gene_body_coverage####
			my $gene_body = "qsub -q all.q -l h=$sort[1] -hold_jid rsem_retry_$nodes_free[0] -N geneBody_ -cwd $cmd_path/geneBody.sh $organism $proj_ID";
			system($gene_body);
			system("qsub -q all.q -hold_jid rseqc_summary_,geneBody_ -N folder_build -cwd $cmd_path/mkdir.sh");
			my $cmd_email="qsub -q all.q -hold_jid folder_build -o ./node_log -e ./node_log -N organize -cwd $cmd_path/organize.sh";
			system($cmd_email);	
		}else{
			system("qsub -q all.q -hold_jid rseqc_summary_ -N folder_build -cwd $cmd_path/mkdir.sh");
			my $cmd_email="qsub -q all.q -hold_jid folder_build -o ./node_log -e ./node_log -N organize -cwd $cmd_path/organize.sh";
			system($cmd_email);	
		}
	}else{
		system("qsub -q all.q -hold_jid rsem_count_tmp_summary -N folder_build -cwd $cmd_path/mkdir_noqc.sh");
		if(defined $email){
			my $cmd_email="qsub -q all.q -hold_jid folder_build -m e -M $email -o ./node_log -e ./node_log -N organize -cwd $cmd_path/organize_noqc.sh";
			system($cmd_email);	
#		print $cmd_email,"\n";
		}else{
		my $cmd_email="qsub -q all.q -hold_jid folder_build -o ./node_log -e ./node_log -N organize -cwd $cmd_path/organize_noqc.sh";
		system($cmd_email);	
		}	
	}
}






__END__


#my @progress = `qstat -f|grep large`;

#my @raw;

#my $useless1 = "all.q\@csclp1-0-1.local";

#foreach my $line (@progress){
#	if ($line =~ /$useless/){
#		next;
#	}else{
#	my @a = split(" ", $line);
#	my @b = split("@",$a[0]);
#	my @c = split('/',$line);
#	my $mem = $c[1]."_".$b[1];
#	push (@raw,$mem);
#	}
#}

#print "Jobs will be ran on @raw","\n\n\n\n";
#my @sort = sort @raw;

#foreach my $line (@sort){
#	$line =~ s/\d_//g;
	#push(@sort_modi,$line);
}


#print scalar(@data),"\n";


#my @node;
#my @mapping_node;
#my @remove_node;
#my @rsem_node;
#my @rseq;

#foreach my $line(@sort){
#	my @a = split("[.]",$line);
#	my @b = split("-",$a[0]);
#	my $genome_load="genome_load_".$b[2];
#	push(@node,$genome_load);
#	my $map= "mapping_".$b[2];
#	push(@mapping_node,$map);
#	my $remove= "removing_reference_".$b[2];
#	push(@remove_node,$remove);
#	my $rsem= "rsem_cal_".$b[2];
#	push(@rsem_node,$rsem);
#	my $rseqc= "rseqc_".$b[2];
#	push(@rseq,$rseqc);
#}
__END__
if($number <= 6 ){
	print "ref genome of $organism is loading to $sort[0]\n";
	####reference genome loading####
	system("qsub -q all.q -l h=$sort[0] -N $node[0] -cwd $cmd_path/loading.sh $organism");
	####Mapping by STAR####
	my $mapping_comman = "qsub -q all.q -l h=$sort[0] -hold_jid $node[0] -N $mapping_node[0] -t 1\-$number -tc 4 -cwd $cmd_path/mapping_star.sh $organism $seq $file";
	system($mapping_comman);
	####Remove reference genome on nodes####
	my $remove_comman = "qsub -q all.q -l h=$sort[0] -hold_jid $mapping_node[0] -N $remove_node[0] -cwd $cmd_path/remove_ref.sh $organism";
	system($remove_comman);
	####Rsem calculation####
	my $resem_comman = "qsub -q all.q -l h=$sort[0] -hold_jid $remove_node[0] -N $rsem_node[0] -t 1\-$number -tc 4 -cwd $cmd_path/rsem_cal.sh $organism $seq $file";
	system($resem_comman);
	####RseQC####
	my $rseqc_comman = "qsub -q all.q -l h=$sort[1] -hold_jid $remove_node[0] -N $rseq[1] -t 1\-$number -tc 2 -cwd $cmd_path/rseqc.sh $organism $seq $file";
	system($rseqc_comman);
	####RSEM_summary####
	my $rsem_summary = "qsub -q all.q -l h=$sort[2] -hold_jid $rsem_node[0] -N rsem_count_tmp_summary -cwd $cmd_path/gene_results_summary.sh $proj_ID";
	system($rsem_summary);
	####Gene_body_coverage####
	my $gene_body = "qsub -q all.q -l h=$sort[3] -hold_jid $rsem_node[0] -N geneBody_ -cwd $cmd_path/geneBody.sh $organism $proj_ID";
	system($gene_body);
	my $rseqc_summary = "qsub -q all.q -hold_jid $rseq[1] -N rseqc_summary_ -cwd $cmd_path/rseqc_summary.sh";
	system($rseqc_summary);
	system("qsub -q all.q -hold_jid rseqc_summary_,geneBody_ -N folder_build -cwd $cmd_path/mkdir.sh");
	if(defined $email){
	my $cmd_email="qsub -q all.q -hold_jid folder_build -m e -M $email -o ./node_log -e ./node_log -N organize -cwd $cmd_path/organize.sh";
	system($cmd_email);	
	print $cmd_email,"\n";
	}else{
	my $cmd_email="qsub -q all.q -hold_jid folder_build -o ./node_log -e ./node_log -N organize -cwd $cmd_path/organize.sh";
	system($cmd_email);	
	}	

}elsif($number >6 && $number <=15){
	print "ref genome of $organism is loading to $sort[0] and $sort[1]\n";
	####reference genome loading####
	system("qsub -q all.q -l h=$sort[0] -N $node[0] -cwd $cmd_path/loading.sh $organism");
	system("qsub -q all.q -l h=$sort[1] -hold_jid $node[0] -N $node[1] -cwd $cmd_path/loading.sh $organism");
	####Mapping by STAR####
	my $div = int($number/2);
	open OUT, ">inputfile1_temp.txt" or die $!;
	open OUTT, ">inputfile2_temp.txt" or die $!;
	print OUT @data[0..($div-1)];
	print OUTT @data[($div)..$#data];
	my $number1 = $#data-$div+1;
	my $mapping_comman1 = "qsub -q all.q -l h=$sort[0] -hold_jid $node[0] -N $mapping_node[0] -t 1\-$div -tc 4 -cwd $cmd_path/mapping_star.sh $organism $seq inputfile1_temp\.txt";
	system($mapping_comman1);
	my $mapping_comman2 = "qsub -q all.q -l h=$sort[1] -hold_jid $node[1] -N $mapping_node[1] -t 1\-$number1 -tc 4 -cwd $cmd_path/mapping_star.sh $organism $seq inputfile2_temp\.txt";
	system($mapping_comman2);
	####Remove reference genome on nodes####
	my $remove_comman_1 = "qsub -q all.q -l h=$sort[0] -hold_jid $mapping_node[0] -N $remove_node[0] -cwd $cmd_path/remove_ref.sh $organism";
	system($remove_comman_1);
	my $remove_comman_2 = "qsub -q all.q -l h=$sort[1] -hold_jid $mapping_node[1],$remove_node[0] -N $remove_node[1] -cwd $cmd_path/remove_ref.sh $organism";
	system($remove_comman_2);
	####Rsem calculation####
	my $resem_comman1 = "qsub -q all.q -l h=$sort[0] -hold_jid $remove_node[0] -N $rsem_node[0] -t 1\-$div -tc 4 -cwd $cmd_path/rsem_cal.sh $organism $seq inputfile1_temp\.txt";
	system($resem_comman1);
	my $resem_comman2 = "qsub -q all.q -l h=$sort[1] -hold_jid $remove_node[1] -N $rsem_node[1] -t 1\-$number1 -tc 4 -cwd $cmd_path/rsem_cal.sh $organism $seq inputfile2_temp\.txt";
	system($resem_comman2);
	####RseQC####
	my $rseqc_comman1 = "qsub -q all.q -l h=$sort[2] -hold_jid $remove_node[0] -N $rseq[2] -t 1\-$div -tc 2 -cwd $cmd_path/rseqc.sh $organism $seq inputfile1_temp\.txt";
	system($rseqc_comman1);
	my $rseqc_comman2 = "qsub -q all.q -l h=$sort[3] -hold_jid $remove_node[1] -N $rseq[3] -t 1\-$number1 -tc 2 -cwd $cmd_path/rseqc.sh $organism $seq inputfile2_temp\.txt";
	system($rseqc_comman2);
	####RSEM_summary####
	my $rsem_summary = "qsub -q all.q -l h=$sort[0] -hold_jid $rsem_node[0],$rsem_node[1] -N rsem_count_tmp_summary -cwd $cmd_path/gene_results_summary.sh $proj_ID";
	system($rsem_summary);
	####Gene_body_coverage####
	my $gene_body = "qsub -q all.q -l h=$sort[1] -hold_jid $rsem_node[0],$rsem_node[1] -N geneBody_ -cwd $cmd_path/geneBody.sh $organism $proj_ID";
	system($gene_body);
	my $rseqc_summary = "qsub -q all.q -hold_jid $rseq[2],$rseq[3] -N rseqc_summary_ -cwd $cmd_path/rseqc_summary.sh";
	system($rseqc_summary);
	system("qsub -q all.q -hold_jid rseqc_summary_,geneBody_ -N folder_build -cwd $cmd_path/mkdir.sh");
	if(defined $email){
	my $cmd_email="qsub -q all.q -hold_jid folder_build -m e -M $email -o ./node_log -e ./node_log -N organize -cwd $cmd_path/organize.sh";
	system($cmd_email);	
	}else{
	my $cmd_email="qsub -q all.q -hold_jid folder_build -o ./node_log -e ./node_log -N organize -cwd $cmd_path/organize.sh";
	system($cmd_email);	
	}	

}elsif($number > 15 && $number <= 40){
	print "ref genome of $organism is loading to $sort[0] and $sort[1] and $sort[2]\n";
	####reference genome loading####
	system("qsub -q all.q -l h=$sort[0] -N $node[0] -cwd $cmd_path/loading.sh $organism");
	system("qsub -q all.q -l h=$sort[1] -hold_jid $node[0] -N $node[1] -cwd $cmd_path/loading.sh $organism");
	system("qsub -q all.q -l h=$sort[2] -hold_jid $node[1] -N $node[2] -cwd $cmd_path/loading.sh $organism");
	####Mapping by STAR####
	my $div = int($number/3);
	open OUT, ">inputfile1_temp.txt" or die $!;
	open OUTT, ">inputfile2_temp.txt" or die $!;
	open OUTTT, ">inputfile3_temp.txt" or die $!;
	print OUT @data[0..($div-1)];
	print OUTT @data[$div..($div*2-1)];
	print OUTTT @data[($div*2)..$#data];
	my $number1 = $#data-$div*2+1;
	my $mapping_comman1 = "qsub -q all.q -l h=$sort[0] -hold_jid $node[0] -N $mapping_node[0] -t 1\-$div -tc 4 -cwd $cmd_path/mapping_star.sh $organism $seq inputfile1_temp\.txt";
#	print $mapping_comman1,"\n";
	system($mapping_comman1);
	my $mapping_comman2 = "qsub -q all.q -l h=$sort[1] -hold_jid $node[1] -N $mapping_node[1] -t 1\-$div -tc 4 -cwd $cmd_path/mapping_star.sh $organism $seq inputfile2_temp\.txt";
	system($mapping_comman2);
	my $mapping_comman3 = "qsub -q all.q -l h=$sort[2] -hold_jid $node[2] -N $mapping_node[2] -t 1\-$number1 -tc 4 -cwd $cmd_path/mapping_star.sh $organism $seq inputfile3_temp\.txt";
	system($mapping_comman3);
	####Remove reference genome on nodes####
	my $remove_comman_1 = "qsub -q all.q -l h=$sort[0] -hold_jid $mapping_node[0] -N $remove_node[0] -cwd $cmd_path/remove_ref.sh $organism";
	system($remove_comman_1);
	my $remove_comman_2 = "qsub -q all.q -l h=$sort[1] -hold_jid $mapping_node[1],$remove_node[0] -N $remove_node[1] -cwd $cmd_path/remove_ref.sh $organism";
	system($remove_comman_2);
	my $remove_comman_3 = "qsub -q all.q -l h=$sort[2] -hold_jid $mapping_node[2],$remove_node[1] -N $remove_node[2] -cwd $cmd_path/remove_ref.sh $organism";
	system($remove_comman_3);
	####Rsem calculation####
	my $resem_comman1 = "qsub -q all.q -l h=$sort[0] -hold_jid $remove_node[0] -N $rsem_node[0] -t 1\-$div -tc 4 -cwd $cmd_path/rsem_cal.sh $organism $seq inputfile1_temp\.txt";
	system($resem_comman1);
	my $resem_comman2 = "qsub -q all.q -l h=$sort[1] -hold_jid $remove_node[1] -N $rsem_node[1] -t 1\-$div -tc 4 -cwd $cmd_path/rsem_cal.sh $organism $seq inputfile2_temp\.txt";
	system($resem_comman2);
	my $resem_comman3 = "qsub -q all.q -l h=$sort[2] -hold_jid $remove_node[2] -N $rsem_node[2] -t 1\-$number1 -tc 4 -cwd $cmd_path/rsem_cal.sh $organism $seq inputfile3_temp\.txt";
	system($resem_comman3);
	####RseQC####
	my $rseqc_comman1 = "qsub -q all.q -l h=$sort[3] -hold_jid $remove_node[0] -N $rseq[3] -t 1\-$div -tc 2 -cwd $cmd_path/rseqc.sh $organism $seq inputfile1_temp\.txt";
	system($rseqc_comman1);
	my $rseqc_comman2 = "qsub -q all.q -l h=$sort[1] -hold_jid $remove_node[1],$rsem_node[1] -N $rseq[1] -t 1\-$div -tc 2 -cwd $cmd_path/rseqc.sh $organism $seq inputfile2_temp\.txt";
	system($rseqc_comman2);
	my $rseqc_comman3 = "qsub -q all.q -l h=$sort[0] -hold_jid $rsem_node[0] -N $rseq[0] -t 1\-$number1 -tc 2 -cwd $cmd_path/rseqc.sh $organism $seq inputfile3_temp\.txt";
	system($rseqc_comman3);
	####RSEM_summary####
	my $rsem_summary = "qsub -q all.q -l h=$sort[0] -hold_jid $rsem_node[0],$rsem_node[1],$rsem_node[2] -N rsem_count_tmp_summary -cwd $cmd_path/gene_results_summary.sh $proj_ID";
	system($rsem_summary);
	####Gene_body_coverage####
	my $gene_body = "qsub -q all.q -l h=$sort[1] -hold_jid $rsem_node[0],$rsem_node[1],$rsem_node[2] -N geneBody_ -cwd $cmd_path/geneBody.sh $organism $proj_ID";
	system($gene_body);
	my $rseqc_summary = "qsub -q all.q -hold_jid $rseq[3],$rseq[1],$rseq[0] -N rseqc_summary_ -cwd $cmd_path/rseqc_summary.sh";
	system($rseqc_summary);
	system("qsub -q all.q -hold_jid rseqc_summary_,geneBody_ -N folder_build -cwd $cmd_path/mkdir.sh");
	if(defined $email){
	my $cmd_email="qsub -q all.q -hold_jid folder_build -m e -M $email -o ./node_log -e ./node_log -N organize -cwd $cmd_path/organize.sh";
	system($cmd_email);	
	}else{
	my $cmd_email="qsub -q all.q -hold_jid folder_build -o ./node_log -e ./node_log -N organize -cwd $cmd_path/organize.sh";
	system($cmd_email);	
	}	

}elsif($number > 40){
	####Reference genome loading####
	print "ref genome of $organism is loading to $sort[0] and $sort[1] and $sort[2] and $sort[3]\n";
	system("qsub -q all.q -l h=$sort[0] -N $node[0] -cwd $cmd_path/loading.sh $organism");
	system("qsub -q all.q -l h=$sort[1] -hold_jid $node[0] -N $node[1] -cwd $cmd_path/loading.sh $organism");
	system("qsub -q all.q -l h=$sort[2] -hold_jid $node[1] -N $node[2] -cwd $cmd_path/loading.sh $organism");
	system("qsub -q all.q -l h=$sort[3] -hold_jid $node[2] -N $node[3] -cwd $cmd_path/loading.sh $organism");
	####Mapping by STAR####
	my $div = int($number/4);
	open OUT, ">inputfile1_temp.txt" or die $!;
	open OUTT, ">inputfile2_temp.txt" or die $!;
	open OUTTT, ">inputfile3_temp.txt" or die $!;
	open OUTTTT, ">inputfile4_temp.txt" or die $!;
	print OUT @data[0..($div-1)];
	print OUTT @data[$div..($div*2-1)];
	print OUTTT @data[($div*2)..($div*3-1)];
	print OUTTTT @data[($div*3)..$#data];
	my $number1 = $#data-$div*3+1;
	my $mapping_comman1 = "qsub -q all.q -l h=$sort[0] -hold_jid $node[0] -N $mapping_node[0] -t 1\-$div -tc 4 -cwd $cmd_path/mapping_star.sh $organism $seq inputfile1_temp\.txt";
	system($mapping_comman1);
	my $mapping_comman2 = "qsub -q all.q -l h=$sort[1] -hold_jid $node[1] -N $mapping_node[1] -t 1\-$div -tc 4 -cwd $cmd_path/mapping_star.sh $organism $seq inputfile2_temp\.txt";
	system($mapping_comman2);
	my $mapping_comman3 = "qsub -q all.q -l h=$sort[2] -hold_jid $node[2] -N $mapping_node[2] -t 1\-$div -tc 4 -cwd $cmd_path/mapping_star.sh $organism $seq inputfile3_temp\.txt";
	system($mapping_comman3);
	my $mapping_comman4 = "qsub -q all.q -l h=$sort[3] -hold_jid $node[3] -N $mapping_node[3] -t 1\-$number1 -tc 4 -cwd $cmd_path/mapping_star.sh $organism $seq inputfile4_temp\.txt";
	system($mapping_comman4);
	####Remove reference genome on nodes####
	my $remove_comman_1 = "qsub -q all.q -l h=$sort[0] -hold_jid $mapping_node[0] -N $remove_node[0] -cwd $cmd_path/remove_ref.sh $organism";
	system($remove_comman_1);
	my $remove_comman_2 = "qsub -q all.q -l h=$sort[1] -hold_jid $mapping_node[1],$remove_node[0] -N $remove_node[1] -cwd $cmd_path/remove_ref.sh $organism";
	system($remove_comman_2);
	my $remove_comman_3 = "qsub -q all.q -l h=$sort[2] -hold_jid $mapping_node[2],$remove_node[1] -N $remove_node[2] -cwd $cmd_path/remove_ref.sh $organism";
	system($remove_comman_3);
	my $remove_comman_4 = "qsub -q all.q -l h=$sort[3] -hold_jid $mapping_node[3],$remove_node[2] -N $remove_node[3] -cwd $cmd_path/remove_ref.sh $organism";
	system($remove_comman_4);
	####Rsem calculation####
	my $resem_comman1 = "qsub -q all.q -l h=$sort[0] -hold_jid $remove_node[0] -N $rsem_node[0] -t 1\-$div -tc 4 -cwd $cmd_path/rsem_cal.sh $organism $seq inputfile1_temp\.txt";
	system($resem_comman1);
	my $resem_comman2 = "qsub -q all.q -l h=$sort[1] -hold_jid $remove_node[1] -N $rsem_node[1] -t 1\-$div -tc 4 -cwd $cmd_path/rsem_cal.sh $organism $seq inputfile2_temp\.txt";
	system($resem_comman2);
	my $resem_comman3 = "qsub -q all.q -l h=$sort[2] -hold_jid $remove_node[2] -N $rsem_node[2] -t 1\-$div -tc 4 -cwd $cmd_path/rsem_cal.sh $organism $seq inputfile3_temp\.txt";
	system($resem_comman3);
	my $resem_comman4 = "qsub -q all.q -l h=$sort[3] -hold_jid $remove_node[3] -N $rsem_node[3] -t 1\-$number1 -tc 4 -cwd $cmd_path/rsem_cal.sh $organism $seq inputfile4_temp\.txt";
	system($resem_comman4);
	####RseQC####
	my $rseqc_comman1 = "qsub -q all.q -l h=$sort[3] -hold_jid $remove_node[0] -N $rseq[3] -t 1\-$div -tc 2 -cwd $cmd_path/rseqc.sh $organism $seq inputfile1_temp\.txt";
	system($rseqc_comman1);
	my $rseqc_comman2 = "qsub -q all.q -l h=$sort[2] -hold_jid $remove_node[1],$rsem_node[2] -N $rseq[2] -t 1\-$div -tc 2 -cwd $cmd_path/rseqc.sh $organism $seq inputfile2_temp\.txt";
	system($rseqc_comman2);
	my $rseqc_comman3 = "qsub -q all.q -l h=$sort[0] -hold_jid $rsem_node[0] -N $rseq[0] -t 1\-$div -tc 2 -cwd $cmd_path/rseqc.sh $organism $seq inputfile3_temp\.txt";
	system($rseqc_comman3);
	my $rseqc_comman4 = "qsub -q all.q -l h=$sort[1] -hold_jid $rsem_node[1] -N $rseq[1] -t 1\-$number1 -tc 2 -cwd $cmd_path/rseqc.sh $organism $seq inputfile4_temp\.txt";
	system($rseqc_comman4);
	####RSEM_summary####
	my $rsem_summary = "qsub -q all.q -l h=$sort[3] -hold_jid $rseq[3],$rsem_node[0],$rsem_node[1],$rsem_node[2],$rsem_node[3] -N rsem_count_tmp_summary -cwd $cmd_path/gene_results_summary.sh $proj_ID";
	system($rsem_summary);
	####Gene_body_coverage####
	my $gene_body = "qsub -q all.q -l h=$sort[0] -hold_jid $rseq[0],$rsem_node[0],$rsem_node[1],$rsem_node[2],$rsem_node[3] -N geneBody_ -cwd $cmd_path/geneBody.sh $organism $proj_ID";
	system($gene_body);
	my $rseqc_summary = "qsub -q all.q -hold_jid $rseq[3],$rseq[1],$rseq[0],$rseq[2] -N rseqc_summary_ -cwd $cmd_path/rseqc_summary.sh";
	system($rseqc_summary);
	system("qsub -q all.q -hold_jid rseqc_summary_,geneBody_ -N folder_build -cwd $cmd_path/mkdir.sh");
	if(defined $email){
	my $cmd_email="qsub -q all.q -hold_jid folder_build -m e -M $email -o ./node_log -e ./node_log -N organize -cwd $cmd_path/organize.sh";
	system($cmd_email);	
	}else{
	my $cmd_email="qsub -q all.q -hold_jid folder_build -o ./node_log -e ./node_log -N organize -cwd $cmd_path/organize.sh";
	system($cmd_email);	
	}	
}








