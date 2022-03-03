#!/usr/bin/perl -w
use strict;
use Data::Dumper qw(Dumper);
use Getopt::Long;
use Term::ANSIColor   qw(:constants);
$Term::ANSIColor::AUTORESET=1;

=head1 NAME

RNA-seq Mapping(STAR) QC RSEM pipeline v3

=head1 DESCRIPTION

This pipeline integrats the Mapping, gene counts/tpm by RSEM and RseQC

This pipeline is compatible for reads of "single-end" and "paired-end" which is specified by the option "-t".

=head1 USAGE

In the folder with only "fastq.gz" files:

nohup perl Mapping_QC_Auto_v3.pl -t <SE|PE> -o <Human_mRNA|Mouse_mRNA|Human_totalRNA|Mouse_totalRNA|Rat> -p <project_ID> -n <1,2,3,...> -qc -gb > projectid.log.txt 2>&1 &

example:
nohup perl Mapping_QC_Auto_v3.pl -t SE -o Mouse_mRNA -p AA-3370--06--21--2017 -n 23,24,25,26,27 -qc > AA-3370--06--21--2017.log.txt 2>&1 &


=head1 REQUIREMENT

- Perl 5

- perl module: Getopt::Long

=head1 OPTIONS

Running options:

#-e or --email [optional]
#       Provide your email address if you would like to be notified after jobs completed.

-t or --type [required if no samplesheet supplied]
        The sequencing type is "single end (SE)" or "paired end (PE)"

-o or --organism
        the reference genoem: Human or Mouse

-p or --project
        project ID for this run

-qc or --qualitycontrol
                whether run quality contorl (RSeQC) or not after mapping

-gb or --genebody
                whether run genebody test to check 3' or 5' bias

-n or --nodes
                which nodes the jobs will be ran on the all.q.The number of jobs shouldn't smaller than the number of nodes used

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

sub execute {
    my ($cmd) = @_;
    print "$cmd\n";
    system($cmd);
}

$ENV{'SGE_ENABLE_COREDUMP'}='yes';

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
                execute("ls *fastq.gz > input_file_fastq_tmp.txt");
                open IN, "input_file_fastq_tmp.txt" or die $!;
                my @data = <IN>;
                foreach my $line (@data){
                        chomp $line;
                        if ($line =~ /_R1/){
                                my @a = split("_R1",$line);
#                               substr($a[0],-1) = "";
                                my $new_name = join(".",$a[0],"R1","fastq.gz");
                                if($new_name eq $line){
                                        next;
                                }else{
                                        my $cmd = "mv $line $new_name";
                                        execute($cmd);
                                }
                        }elsif($line =~ /_R2/){
                                my @a = split("_R2",$line);
#                               substr($a[0],-1) = "";
                                my $new_name = join(".",$a[0],"R2","fastq.gz");
                                if($new_name eq $line){
                                        next;
                                }else{
                                        my $cmd = "mv $line $new_name";
                                        execute($cmd);
                                }
                                #       }else{
                                #       die "please provide paired-end reads with right format when using the "PE" argument";
                        }else{
                                die "please provide paired-end reads!\n";
                        }
                }

#               execute("rm input_file_fastq_tmp.txt");
		my $cmd_fastq_file = "ls *.fastq.gz|sed 's/\\.R[1\|2].*fastq.gz//g'|sort -u > input_fastq.txt";
		execute($cmd_fastq_file);
#                execute("ls *.fastq.gz|sed 's/\.R[1\|2].*fastq.gz//g'|sort -u > input_fastq.txt");
        }elsif($type eq "SE"){
                #       execute("ls *.fastq|sed 's/.[fastq|fq]//g' > input_fastq.txt");
                my $temp_cmd = "ls *.fastq.gz|sed 's/\\\(.fastq.gz\\|.fq.gz\\)//g' > input_fastq.txt";
                execute($temp_cmd);
                #       execute("ls *.fastq|sed 's/\(.fastq\|.fq\)//g' > input_fastq.txt");
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
                        my $pair1 = $line.".R1".".fastq.gz";
                        my $pair2 = $line.".R2".".fastq.gz";
#                       print $pair1,"\t",$pair2,"\n";

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

#__END__
if (defined $email){
        print "The complete notification will be sent to: ";
        print RED "$email\n";
}else{
        print GREEN "no email notification after job complete.\n";
}

#my $cmd_path="/common/genomics-core/apps/sequencing_data_distri";
my $cmd_path="/common/genomics-core/apps/mapping_qc_auto";
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

my $nodeindex;
my $index;
		

####Reference genome loading####
#### suggest changing this to a loop to allow a different number of hosts
# print "ref genome of $organism is loading to $sort[0] and $sort[1] and $sort[2] and $sort[3]\n";

for ($index=0;$index<$nodes_number;$index++)
{
	execute("qsub -q all.q -l h=$sort[$index] -N genome_load_$sort[$index] -sync y -cwd $cmd_path/loading-setx_new_reference.sh $organism $sort[$index]");
	
}
####Mapping by STAR####
for ($index=0;$index<$number;$index++)
{
	$nodeindex=$index % $nodes_number;
	print $data[$index],"\n";
#	my $mapping_command = "qsub -q all.q -l h=$sort[$nodeindex] -pe smp 3 -N mapping_$sort[$nodeindex] -cwd $cmd_path/mapping_star_notarray_setx.sh $organism $seq $data[$index]";
	my $mapping_command = "qsub -q all.q -l h=$sort[$nodeindex] -pe smp 10 -N mapping -cwd $cmd_path/mapping_star_notarray_setx_new_reference.sh $organism $seq $data[$index]";
	execute($mapping_command);
}

##########Removeing reference genome
for ($index=0;$index<$nodes_number;$index++)
{
	#execute("qsub -q all.q -l h=$sort[$index] -N remove_ref_$sort[$index] -hold_jid mapping_$sort[$nodeindex] -sync y -cwd $cmd_path/remove_ref.sh $organism $sort[$index]");
	execute("qsub -q all.q -l h=$sort[$index] -N remove_ref -hold_jid mapping -sync y -cwd $cmd_path/remove_ref_new_reference.sh $organism $sort[$index]");
	
}

####Rsem calculation####
for ($index=0;$index<$number;$index++)
{
	$nodeindex=$index % $nodes_number;
	print $data[$index],"\n";
	#execute("qsub -q all.q -N rsem_cal -hold_jid remove_ref_$sort[$nodeindex] -pe smp  -cwd $cmd_path/rsem_cal.sh $organism $seq $data[$index]");
	execute("qsub -q all.q -N rsem_cal -hold_jid remove_ref -pe smp 5 -cwd $cmd_path/rsem_cal_new_reference.sh $organism $seq $data[$index]");
	
}
####RSEM_summary####
my $rsem_summary = "qsub -q all.q -hold_jid rsem_cal -N rsem_count_tmp_summary -cwd $cmd_path/gene_results_summary.sh $proj_ID";
execute($rsem_summary);

####FASTQC####
for ($index=0;$index<$number;$index++)
{
	my $fastqc = "qsub -q all.q -hold_jid rsem_count_tmp_summary -pe smp 5 -N fastqc -cwd $cmd_path/fastqc.sh $data[$index]";
	execute($fastqc);
}

if($qc == 1){
####RseQC####
	for ($index=0;$index<$number;$index++){
		$nodeindex=$index % $nodes_number;
		print $data[$index],"\n";
        	my $rseqc_comman1 = "qsub -q all.q -hold_jid fastqc -pe smp 5 -N rseqc -cwd $cmd_path/rseqc_new_reference.sh $organism $seq $data[$index]";
        	execute($rseqc_comman1);
	}


####RseQC summary                
	my $rseqc_summary = "qsub -q all.q -hold_jid rseqc -N rseqc_summary -cwd $cmd_path/rseqc_summary.sh $proj_ID";
	execute($rseqc_summary);

####MultiQC ####
	my $multiqc_summary = "qsub -q all.q -hold_jid rseqc_summary -N multiqc_generate -cwd $cmd_path/multiqc.sh $proj_ID";
	execute($multiqc_summary);
}


if ($genebody == 1){
####Gene_body_coverage####
	my $gene_body = "qsub -q all.q -hold_jid rsem_cal -N geneBody -cwd $cmd_path/geneBody.sh $organism $proj_ID";
	execute($gene_body);
}


####creater folder
if($qc==1 && $genebody ==1){
	execute("qsub -q all.q -hold_jid geneBody,rseqc_summary -N folder_build -cwd $cmd_path/mkdir.sh");
	my $cmd_email="qsub -q all.q -hold_jid folder_build,multiqc_generate -o ./node_log -e ./node_log -N organize -cwd $cmd_path/organize.sh";
	execute($cmd_email);
}elsif($qc==1 && $genebody !=1){
	execute("qsub -q all.q -hold_jid rseqc_summary -N folder_build -cwd $cmd_path/mkdir.sh");
	my $cmd_email="qsub -q all.q -hold_jid folder_build,multiqc_generate -o ./node_log -e ./node_log -N organize -cwd $cmd_path/organize.sh";
	execute($cmd_email);
}elsif($qc !=1 && $genebody ==1){
	execute("qsub -q all.q -hold_jid geneBody -N folder_build -cwd $cmd_path/mkdir.sh");
	my $cmd_email="qsub -q all.q -hold_jid folder_build,multiqc_generate -o ./node_log -e ./node_log -N organize -cwd $cmd_path/organize.sh";
	execute($cmd_email);
}elsif($qc !=1 && $genebody !=1){
	execute("qsub -q all.q -hold_jid fastqc -N folder_build -cwd $cmd_path/mkdir.sh");
	my $cmd_email="qsub -q all.q -hold_jid folder_build,multiqc_generate -o ./node_log -e ./node_log -N organize -cwd $cmd_path/organize.sh";
	execute($cmd_email);
}
	

