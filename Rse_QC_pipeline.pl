#!/usr/bin/perl -w
use strict;
use Getopt::Long;

#########STAR Alignment ######
#STAR --genomeDir /stf/home/cluster/wangyiz/genomics/NGS/refGenomes/STAR/GRCh38 --outSAMmode Full --outSAMunmapped Within --outFilterType Normal --outSAMattributes All --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --sjdbScore 1 --runThreadN 8 --genomeLoad NoSharedMemory --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --outSAMheaderHD @HD --outBAMsortingThreadN 8 --outFileNamePrefix $1 --readFilesIn $1.R1.fastq

##########Other##########

#Rse_QC_pipeline.pl -i <bam_file> -t <SE/PE> -d <OUTPUT_DIR> -p <PREFIX> -r <QC_reference_files.txt> [-AS] [-tin]


=head1 NAME

RNA-seq Quality Control Pipeline v1

=head1 DESCRIPTION

This pipeline is based on "RSeQC v2.6.3" (http://rseqc.sourceforge.net/). The pipeline includs the module of reads mapping statistics calculation (bam_stat.py), 5'/3'bias calculation (geneBody_coverage.py), insert size calculation (inner_distance.py), junction saturation calculation (junction_saturation.py), read distribution (CDS,UTR,intron,ribosomal RNA (split_bam.py) by read_distribution.py ), GC content (read_GC.py), read quality (read_quality.py) and transcript integrity (tin.py). All output and figures would be summarized and combined into final_results_summary.txt and final_results_figures.pdf. Output from the RseQC are saved at prefix_temp folder.

This pipeline is compatible for reads of "single-end" and "paired-end" which is specified by the option "-t".

=head1 USAGE

Rse_QC_pipeline.pl -i BAM_File -t <SE|PE> -d <OUTPUT_DIR> -p <PREFIX> -r <QC_reference_files.txt> [-AS] [-tin]

example: Rse_QC_pipeline.pl -i *.bam -t [PE|SE] -d ./results -p * -r <QC_reference_files.txt> [-AS] [-tin]

=head1 REQUIREMENT

- Perl 5

- perl module: Getopt::Long

- input bam files should be sorted and index by samtools and saved together with .bai index files

- Running on servers, please export such paths in the submission bash script	


export PATH=/home/wangyiz/genomics/anaconda2/bin:$PATH

=head1 OPTIONS

Running options:

-i or --input_single INPUT_FILE  [required]
	Single alignment file in BAM format.

-d or --output_dir OUTPUT_DIR 
	path to the directory where the final results should be stored

-p or --prefix PREFIX [required]
	prefix for the final results 

-AS or --alternative_splicing 
	If the project is relative to alternative splicing analysis, this option is recommended to turn on for checking the junction saturation. It is turned off as default.

-tin or --trans_inte
	Option to calculate the transcript integrety for each gene and overall transcriptome. This is helpful for evaluating the quality of RNA samples - whether RNA degraged or not. Value threshhold is from 0 to 100 - the higher the value, the better quality of RNA sample. NOTE: THIS OPTION WILL TAKE VERY VERY LONG TIME!!!
	
-h or --help
	Help information 


=head1 OUTPUT FILES

[prefix]_final_results_summary.txt
	summarized results of mapping stats, including total number of raw reads, total number of uniquely mapped reads, percetage of uniquely mapped reads, total number of duplicated mapped reads, percetage of multiple mapping reads, total percetage of mapping, mapping percetage on different features (CDS, UTR, Intron, rRNA) and mean of transcript integrity(if -tin was specificed). 

[prefix]_final_results_figures.pdf
	Figures of GC content, 5'/3' mapping bias checking, insert size (PE specificed), junction saturation calculation (-AS specificed),read quality plot and read quality heatmap. 

*_temp 
	Folders saving all files produced during the analysis

=head1 AUTHOR

If you have any question, please email: yizhou.wang@cshs.org

Genomics CORE, 03/05/2016

=cut


use Cwd 'abs_path';

#my ($type,$bam,$bam_dir,$output_dir,$prefix,$bed,$genomebed,$ribobed,$hsbed);
my ($type, $bam, $bam_dir, $output_dir, $prefix, $reference);
my $AS = 0;
my $tin_open = 0;

GetOptions(
	'type|t=s' => \$type,
	'input_single|i=s' => \$bam,
	'input_list|ilist=s' => \$bam_dir,
	'output_dir|d=s' => \$output_dir,
	'prefix|p=s' => \$prefix,
	'reference|r=s' => \$reference,
#	'genome_bed|gb=s' => \$genomebed,
#	'rRNA_bed|rb=s' => \$ribobed,
#	'housekeeping_bed|hb=s' =>\$hsbed,
	'alternative_splicing|AS!' => \$AS,
	'trans_inte|tin!' => \$tin_open,
	'help|h' => sub{exec('perldoc',$0);
									exit(0);},
);


#unless($bam||$bam_dir){
#	die("Require an input file via -i or --input_single for single .bam or -ilist or --input_list for a list of .bam files on the cmdline. See usage with -h\n");
#}else{
#	$file = $bam;
#}

if ($bam){
	my $abs_path_bamfile=abs_path($bam);#obtain the abs path of the file
	my $abs_path_folder = abs_path($output_dir);
#	my $abs_path_genomebed = abs_path($genomebed);
#	my $abs_path_ribobed = abs_path($ribobed);
#	my $abs_path_hsbed = abs_path($hsbed);
#	my @all_files = ($abs_path_bamfile,$abs_path_genomebed,$abs_path_ribobed,$abs_path_hsbed);
	my $reference_bed = abs_path($reference);
	my $bam_name = $bam;
	my $rRNA = 0;
	my $tRNA = 0;
	open REF, $reference_bed or die $!;
	my @ref = <REF>;
	my %ref_bed;
	my @junk_RNA;
	my %junk_RNA_hash;
	foreach my $line (@ref){
		chomp $line;
		my @a = split("\t",$line);
		if ($a[0] =~ /optional/){
			$a[0] =~ s/_optional//g;
			push(@junk_RNA,$a[0]);
		}
		if ($a[0] =~ /required/){
			$a[0] =~ s/_required//g;
		}
		$ref_bed{$a[0]} = $a[1];
	}
	$bam_name =~ s/\.bam//g;
	my $results_folder = $prefix;
	my $inputfile_name = $prefix;
	$prefix = $prefix."_tempfile";
#	File_Exist(@all_files);
	if( -e $abs_path_folder){
		print "output folder has already been built up\n";
	}else{
		mkdir $abs_path_folder;
	}	
	chdir $abs_path_folder;
	system("mkdir $results_folder");
	chdir $results_folder;
	system("samtools view -F 0x0100 $abs_path_bamfile -b > $prefix.unique_primary.bam");
	system("samtools view -s 0.2 -b $abs_path_bamfile -b > $prefix.20_per.bam");
	if ($type eq "PE"){
		system("inner_distance.py -i $prefix.unique_primary.bam -o $prefix.inner_distance -r $ref_bed{\"all\"}");
		system('echo "Done for inner distance calculation!!!"');
		system("samtools view -f 64 -b $prefix.unique_primary.bam > $prefix.first_strand.bam");
		system("samtools view -f 128 -b $prefix.unique_primary.bam > $prefix.second_strand.bam");
		system("bam_stat.py -i $prefix.first_strand.bam > $prefix.first_strand.bam.stat.txt");
		system("bam_stat.py -i $prefix.second_strand.bam > $prefix.second_strand.bam.stat.txt");
		system('echo "Done for bam_stat!!!"');
	}elsif($type eq "SE"){
		system("bam_stat.py -i $prefix.unique_primary.bam > $prefix.bam.stat.txt");
		system('echo "Done for bam_stat!!!"');
		print "Single-end reads provided, so no insert-size calculation.\n";
	}
	if($AS == 1){
		system("junction_saturation.py -i $abs_path_bamfile -r $ref_bed{\"all\"} -o $prefix.junction_saturation");
		system('echo "Done for junction saturation calculation!!!"');
	}
	system("read_distribution.py -i $prefix.unique_primary.bam -r $ref_bed{\"all\"} > $prefix.read_distribution.txt");
	system('echo "Done for read distribution calculation!!!"');
	system("read_GC.py -i $prefix.20_per.bam -o $prefix.GC");
	system('echo "Done for GC content calculation!!!"');
	system("read_quality.py -i $prefix.20_per.bam -o $prefix.readquality");
	system('echo "Done for read quality calculation!!!"');
	foreach my $d (@junk_RNA){
		chomp $d;
		my $command = "coverageBed -abam $prefix.unique_primary.bam -b $ref_bed{$d} -counts -f 0.01 -nonamecheck\|awk \'{print \$NF}\'\|awk \'{s+=\$1} END {print s}\'";
		my $junk_RNA_number = `$command`;
		chomp $junk_RNA_number;
		$junk_RNA_hash{$d} = $junk_RNA_number;
		system("echo \"Done for $d calculation!!!\"")
	}
	if($tin_open == 1){
		system("tin.py -i $abs_path_bamfile -r $ref_bed{\"all\"}");
		system('echo "Done for transcript integrity calculation!!!"');
	}
	#system("geneBody_coverage.py -r $ref_bed{\"housekeeping\"} -i $abs_path_bamfile -o $prefix");
	if(-z "$prefix.read_distribution.txt"){
		print "Read distribution file is NULL probably due to memmory allocation problem!! We need to recalculate the read distribution again!";
		my $count_rseqc=5;
		do{
			system("read_distribution.py -i $prefix.unique_primary.bam -r $ref_bed{\"all\"} > $prefix.read_distribution.txt");
		$count_rseqc--;
		}until(-s "$prefix.read_distribution.txt" || $count_rseqc==0);
	}else{
		print "Read distribution file is no NULL! Continue!\n";
	}
	system('echo "ALL DONE!!! GOOD LUCKY FOR YOUR DATA and HAVE A NICE DAY!!!"');
	
	my $tin = 0;
	if($tin_open == 1){
		my $tin_file = `ls *.summary.*`;
		open INNNN, "$tin_file" or die $!;
		my @data_tin = <INNNN>;
		shift @data_tin;
		foreach my $line (@data_tin){
			chomp $line;
			$line =~ s/\s+/\t/g;
			my @a = split("\t",$line);
			$tin = sprintf("%.2f",$a[2]);
		}
	}

	open INN, "$prefix.read_distribution.txt" or die $!;
#	open INNN, "$prefix.rRNA.txt" or die $!;
	
	my @data_cds = <INN>;
#	my @data_ribo = <INNN>;

############this is to calculate basic stats of reads mapping#########
	my $total = 0;
	my $total_first = 0;
	my $total_second = 0;
	my $unmapped_first = 0;
	my $unmapped_second = 0;
	my $unmapped = 0;
	my $unique = 0;


	my $read1 = 0;
	my $read2 = 0;
	my $duplicate_first = 0;
	my $duplicate_second = 0;
	my $mapped = 0;
	my $mapped_per = 0;
	my $unique_per = 0;
	my $duplicate = 0;
	my $duplicate_per = 0;
	my $proper_pair = 0;

	if ($type eq "PE"){
		open IN_FIRST, "$prefix.first_strand.bam.stat.txt" or die $!;
		open IN_SECOND, "$prefix.second_strand.bam.stat.txt" or die $!;
		my @data_first = <IN_FIRST>;
		my @data_second = <IN_SECOND>;
		
		foreach my $line (@data_first){
			chomp $line;
			$line =~ s/\s+/\t/g;
			if ($line =~ /Total\trecords/){
				my @a = split("\t",$line);
				$total_first = $a[2];
			}
			if ($line =~ /Unmapped\treads/){
				my @a = split("\t",$line);
				$unmapped_first = $a[2];
			}
			if ($line =~ /Read-1/){
				my @a = split("\t",$line);
				$read1 = $a[1];
			}
			if ($line =~ /non-unique/){
				my @a = split("\t",$line);
				$duplicate_first = $a[4];
			}
			if ($line =~ /pairs/){
				my @a = split("\t",$line);
				$proper_pair = $a[5];
			}
		}

		foreach my $line (@data_second){
			chomp $line;
			$line =~ s/\s+/\t/g;
			if ($line =~ /Total\trecords/){
				my @a = split("\t",$line);
				$total_second = $a[2];
			}
			if ($line =~ /Unmapped\treads/){
				my @a = split("\t",$line);
				$unmapped_second = $a[2];
			}
			if ($line =~ /Read-2/){
				my @a = split("\t",$line);
				$read2 = $a[1];
			}
			if ($line =~ /non-unique/){
				my @a = split("\t",$line);
				$duplicate_second = $a[4];
			}
		}
	$mapped = $total_first+$total_second-$unmapped_first-$unmapped_second;
	$mapped_per = sprintf("%.2f",100*$mapped/($total_first+$total_second));
	$unique_per = sprintf("%.2f",100*($read1+$read2)/($total_first+$total_second));
	$duplicate_per = sprintf("%.2f",100*($duplicate_first+$duplicate_second)/($total_first+$total_second)); 

	}elsif($type eq "SE"){	
		open IN, "$prefix.bam.stat.txt" or die $!;
		my @data_all = <IN>;

		foreach my $line (@data_all){
			chomp $line;
			$line =~ s/\s+/\t/g;
			if ($line =~ /Total\trecords/){
				my @a = split("\t",$line);
				$total = $a[2];
			}
			if ($line =~ /Unmapped\treads/){
				my @a = split("\t",$line);
				$unmapped = $a[2];
			}
			if ($line =~ /\(unique\)/){
				my @a = split("\t",$line);
				$unique = $a[4];
			}
		}
	$mapped = $total - $unmapped;
	$mapped_per = sprintf("%.2f",100*$mapped/$total);
	$unique_per = sprintf("%.2f",100*$unique/$total);
	$duplicate = $mapped - $unique;
	$duplicate_per = sprintf("%.2f",100*$duplicate/$total);
	}

#########this is to calculate ribosomal RNA########
#	shift @data_ribo;

#	my $sum = 0;
#	my $ribosomal = 0;

#	foreach my $line (@data_ribo){
#		chomp $line;
#		if ($line =~ /\.in\./){
#			my @a = split(":",$line);
#			$ribosomal = $a[1];
#			$sum = $sum + $ribosomal;
#		}elsif($line =~ /\.ex\./){
#			my @a = split(":",$line);
#			$sum = $sum + $a[1];
#		}
#	}
#	my $result_ribo = sprintf("%.2f",100 * $ribosomal/$sum);
	open JUNK_RNA,">junk_RNA_tempfile.txt" or die $!;
	my @junk_RNA_per;
	foreach my $line (@junk_RNA){
		my $value = sprintf("%.2f",100*$junk_RNA_hash{$line}/$mapped)."%";
		push (@junk_RNA_per, $value);
		print JUNK_RNA join("\t",$line,$junk_RNA_hash{$line},$value),"\n";
	}
	close JUNK_RNA;

system("Rscript /common/genomics-core/apps/mapping_qc_auto/junkRNA_muti.R $mapped $prefix");

###########this is to calculate CDS, UTR and Intron########
	my $total_tag = 0;
	my $cds = 0;
	my $UTR = 0;
	my $intron = 0;

	foreach my $line (@data_cds){
		chomp $line;
		$line =~ s/\s+/\t/g;
		if ($line =~ /Total\tTags/){
			my @a = split("\t",$line);
			$total_tag = $a[2];
		}
		if ($line =~ /CDS_Exons/){
			my @a = split("\t",$line);
			$cds = $a[2];
		}
		if ($line =~ /5\'UTR_Exons/){
			my @a = split("\t",$line);
			$UTR = $UTR+$a[2];
		}
		if ($line =~ /3\'UTR_Exons/){
			my @a = split("\t",$line);
			$UTR = $UTR+$a[2];
		}
		if ($line =~ /Introns/){
			my @a = split("\t",$line);
			$intron = $a[2];
		}
	}

	$cds = sprintf("%.2f",100*$cds/$total_tag);
	$UTR = sprintf("%.2f",100*$UTR/$total_tag);
	$intron = sprintf("%.2f",100*$intron/$total_tag);

###########################################
	my $output_file = $bam_name."_final_results";
	my $output_file_pdf = $bam_name."_final_results_pdf";
	open OUT, ">$output_file" or die $!;

	if($tin_open ==1 ){
		if($type eq "SE"){
			print OUT join("\t","Sample","#_raw_reads","#_unique_reads","%_unique_reads","#_multi_mapping_reads","%_multi_mapping_reads","total%_mapping","CDS","UTR","intron",@junk_RNA,"TIN(median)"),"\n";
			print OUT join("\t",$inputfile_name,$total,$unique,"$unique_per%",$duplicate,"$duplicate_per%","$mapped_per%","$cds%","$UTR%","$intron%",@junk_RNA_per,$tin),"\n";
			close INNNN;
		}elsif($type eq "PE"){
			print OUT join("\t","Sample","#_raw_Read1","#_raw_Read2","#_proper_pairs","#_unique_Reads1","#_unique_Reads2","%_unique_reads_total","#_multi_mapping_Reads1","#_multi_mapping_Read2","%_multi_mapping_reads_total","total%_mapping","CDS","UTR","intron",@junk_RNA,"TIN(median)"),"\n";
			print OUT join("\t",$inputfile_name,$total_first,$total_second,$proper_pair,$read1,$read2,"$unique_per%",$duplicate_first,$duplicate_second,"$duplicate_per%","$mapped_per%","$cds%","$UTR%","$intron%",@junk_RNA_per,$tin),"\n";
			close INNNN;
		}
	}else{
		if($type eq "SE"){
			print OUT join("\t","Sample","#_raw_reads","#_unique_reads","%_unique_reads","#_multi_mapping_reads","%_multi_mapping_reads","total%_mapping","CDS","UTR","intron",@junk_RNA),"\n";
			print OUT join("\t",$inputfile_name,$total,$unique,"$unique_per%",$duplicate,"$duplicate_per%","$mapped_per%","$cds%","$UTR%","$intron%",@junk_RNA_per),"\n";
			close INNNN;
		}elsif($type eq "PE"){
			print OUT join("\t","Sample","#_raw_Read1","#_raw_Read2","#_proper_pairs","#_unique_Reads1","#_unique_Reads2","%_unique_reads_total","#_multi_mapping_Reads1","#_multi_mapping_Read2","%_multi_mapping_reads_total","total%_mapping","CDS","UTR","intron",@junk_RNA),"\n";
			print OUT join("\t",$inputfile_name,$total_first,$total_second,$proper_pair,$read1,$read2,"$unique_per%",$duplicate_first,$duplicate_second,"$duplicate_per%","$mapped_per%","$cds%","$UTR%","$intron%",@junk_RNA_per),"\n";
			close INNNN;
		}
	}

	close OUT;
	close IN;
	close INN;
#	close INNN;

	my $read_dis_figure = "\/common\/genomics-core\/bin\/read_distribution\.R";
	system("cat $read_dis_figure $prefix.*.r > $prefix.all_temp.r");

	open IN, "$prefix.all_temp.r" or die $!;

	my @script = <IN>;

	foreach my $line (@script){
	
		$line =~ s/pdf\(.*//g;
		$line =~ s/dev\.off\(\)//g;
	}

	push @script,"dev.off()";
	unshift @script,"pdf(\"$output_file_pdf\",10,8)\n";


	open OUT, ">$prefix.all_new.r" or die $!;
	print OUT "@script\n";

	close OUT;
	close IN;

#	my $cmd_pdf_combine='convert -density 250 $(ls -rt *pdf) final_results_pdf';	
	my $cmd_temp_folder = 'mkdir '.$results_folder.'_temp';

#	system($cmd_pdf_combine);
	system("Rscript $prefix.all_new.r $output_file");
	system($cmd_temp_folder);

#	my $cmd_txt = "mv $prefix*.txt ".$prefix."_temp";
#	my $cmd_xls = "mv $prefix*.xls ".$prefix."_temp";
	my $cmd_all = "mv *_tempfile* ".$results_folder."_temp";
	if ($tin_open ==1 ){
		my $cmd_sample_xls = "mv *.xls ".$results_folder."_temp";
		my $cmd_sample_summary = "mv *.summary.* ".$results_folder."_temp";
		system($cmd_sample_xls);
		system($cmd_sample_summary);
	}
	my $cmd_rename1 = "mv $output_file ".$results_folder.".final_results_summary.txt";
	my $cmd_rename2 = "mv $output_file_pdf ".$results_folder.".final_results_figures.pdf";

#	system($cmd_txt);
#	system($cmd_xls);
	system($cmd_all);
	system($cmd_rename1);
	system($cmd_rename2);
#	system("rm log.txt");
#}elsif($bam_dir){
#	open IN, $bam_dir or die "No list of bam file was detected!";
#	my @bam_all = <IN>;
#	my $abs_path_folder = abs_path($output_dir);
#	my $abs_path_genomebed = abs_path($genomebed);
#	my $abs_path_ribobed = abs_path($ribobed);
#	my $abs_path_hsbed = abs_path($hsbed);
#	my $abs_path_allbam = abs_path($bam_dir);
#	my @all_rRNA;
#	my @all_files = ($abs_path_genomebed,$abs_path_ribobed,$abs_path_hsbed,$abs_path_allbam);
#	File_Exist(@all_files);
#	File_Exist(@bam_all);
#	if(-e $abs_path_folder){
#		print "output folder has been already built up\n";
#	}else{
#		mkdir $abs_path_folder;
#	}	
#	my @all_distribution;
#	my @all_tin;
#	my @all_basic_stat;
#	my $sample_name;
##	my $flag = 0;
##	my $path = abs_path("\.\/haha");
##	print $path,"\n";
##	my $prefix = "lala";
#	my $output_temp_path = $abs_path_folder.'/'.$prefix.'_temp';
#	if(-e $output_temp_path){
#		print "temp folder has been already built up\n";
#	}else{
#	my $cmd_temp_folder = 'mkdir '.$abs_path_folder.'/'.$prefix.'_temp';
##	print $cmd_temp_folder, "\n";
#	system($cmd_temp_folder);
#	}
##	my $cmd_temp_folder = 'mkdir '.$prefix.'_temp';
##	system($cmd_temp_folder);
#
##	foreach my $single_bam(@bam_all){
##		chomp $single_bam;i
##		my @a = split(/\//,$single_bam);
##		$sample_name = $a[$#a];
##		$hash{$sample_name}=$single_bam;
##	}
#
#
#	foreach my $single_bam(@bam_all){
#		chomp $single_bam;
#		my @a = split(/\//,$single_bam);
#		$sample_name = $a[$#a];
#		$sample_name =~ s/\.bam//g;
#		chdir $abs_path_folder;
#		system("samtools view -s 0.2 -b $single_bam -b > $sample_name.20_per.bam");
#		system("samtools view -F 0x0100 $single_bam -b > $sample_name.unique_primary.bam");
#		system("bam_stat.py -i $sample_name.unique_primary.bam > $sample_name.bam.stat.txt");
#		system('echo "Done for bam_stat!!!"');
#		if ($type eq "PE"){
#			system("inner_distance.py -i $sample_name.unique_primary.bam -o $sample_name -r $abs_path_genomebed");
#			system('echo "Done for inner distance calculation!!!"');
#		}elsif($type eq "SE"){
#			print "Single-end reads provided, so no insert-size calculation.\n";
#		}
#		if($AS == 1){
#				system("junction_saturation.py -i $single_bam -r $abs_path_genomebed -o $sample_name.junction_saturation");
#				system('echo "Done for junction saturation calculation!!!"');
#		}
#		system("read_distribution.py -i $sample_name.unique_primary.bam -r $abs_path_genomebed > $sample_name.read_distribution.txt");
#		system('echo "Done for read distribution calculation!!!"');
##			system("read_GC.py -i $sample_name.unique_primary.bam -o $sample_name.GC");
#		system("read_GC.py -i $sample_name.20_per.bam -o $sample_name.GC");
#		system('echo "Done for GC content calculation!!!"');
##			system("read_quality.py -i $sample_name.unique_primary.bam -o $sample_name.readquality");
#		system("read_quality.py -i $sample_name.20_per.bam -o $sample_name.readquality");
#		system('echo "Done for read quality calculation!!!"');
#		system("split_bam.py -i $sample_name.unique_primary.bam -r $abs_path_ribobed -o $sample_name.rRNA > $sample_name.rRNA.txt");
#		system('echo "Done for ribosomal RNA calculation!!!"');
#		system("tin.py -i $single_bam -r $abs_path_genomebed");
#		system('echo "Done for transcript integrity calculation!!!"');
#		my $tin_file = `ls $sample_name.summary.*`;
#		chomp $tin_file;
#		system("mv $tin_file $sample_name.tin.txt");
#		system('echo "ALL DONE!!! GOOD LUCKY FOR YOUR DATA and HAVE A NICE DAY!!!"');
#
#		open IN, "$sample_name.bam.stat.txt" or die $!;
#		open INN, "$sample_name.read_distribution.txt" or die $!;
#		open INNN, "$sample_name.rRNA.txt" or die $!;
#		open INNNN, "$sample_name.tin.txt" or die $!;
#
#		my @data_all = <IN>;
#		my @data_cds = <INN>;
#		my @data_ribo = <INNN>;
#		my @data_tin = <INNNN>;
#
#############this is to calculate basic stats of reads mapping#########
#
#		my $total = 0;
#		my $unmapped = 0;
#		my $unique = 0;
#
#
#		foreach my $line (@data_all){
#			chomp $line;
#			$line =~ s/\s+/\t/g;
#			if ($line =~ /Total\trecords/){
#				my @a = split("\t",$line);
#				$total = $a[2];
#			}
#			if ($line =~ /Unmapped\treads/){
#				my @a = split("\t",$line);
#				$unmapped = $a[2];
#			}
#			if ($line =~ /\(unique\)/){
#				my @a = split("\t",$line);
#				$unique = $a[4];
#			}
#		}
#
#		my $mapped = $total - $unmapped;
#		my $mapped_per = sprintf("%.2f",100*$mapped/$total);
#		my $unique_per = sprintf("%.2f",100*$unique/$mapped);
#		my $duplicate = $mapped - $unique;
#		my $duplicate_per = sprintf("%.2f",100*$duplicate/$mapped);
#
##########this is to calculate ribosomal RNA########
#		shift @data_ribo;
#
#		my $sum = 0;
#		my $ribosomal = 0;
#
#		foreach my $line (@data_ribo){
#			chomp $line;
#			if ($line =~ /\.in\./){
#				my @a = split(":",$line);
#				$ribosomal = $a[1];
#				$sum = $sum + $ribosomal;
#			}elsif($line =~ /\.ex\./){
#				my @a = split(":",$line);
#				$sum = $sum + $a[1];
#			}
#		}
#		my $result_ribo = sprintf("%.2f",100 * $ribosomal/$sum);
#
############this is to calculate CDS, UTR and Intron########
#		my $total_tag = 0;
#		my $cds = 0;
#		my $UTR = 0;
#		my $intron = 0;
#
#		foreach my $line (@data_cds){
#			chomp $line;
#			$line =~ s/\s+/\t/g;
#			if ($line =~ /Total\tTags/){
#				my @a = split("\t",$line);
#				$total_tag = $a[2];
#			}
#			if ($line =~ /CDS_Exons/){
#				my @a = split("\t",$line);
#				$cds = $a[2];
#			}
#			if ($line =~ /5\'UTR_Exons/){
#				my @a = split("\t",$line);
#				$UTR = $UTR+$a[2];
#			}
#			if ($line =~ /3\'UTR_Exons/){
#				my @a = split("\t",$line);
#				$UTR = $UTR+$a[2];
#			}
#			if ($line =~ /Introns/){
#				my @a = split("\t",$line);
#				$intron = $a[2];
#			}
#		}
#
#		$cds = sprintf("%.2f",100*$cds/$total_tag);
#		$UTR = sprintf("%.2f",100*$UTR/$total_tag);
#		$intron = sprintf("%.2f",100*$intron/$total_tag);
#
##################this is to calculate TIN#############
#		shift @data_tin;
#
#		my $tin = 0;
#
#		foreach my $line (@data_tin){
#			chomp $line;
#			$line =~ s/\s+/\t/g;
#			my @a = split("\t",$line);
#			$tin = sprintf("%.2f",$a[1]);
#		}
#
############################################
##		my $output = $sample_name.".final_results";
##		open OUT, ">$output" or die $!;
#
#		open OUT, ">final_results_summary" or die $!;
#		print OUT join("\t","Sample","#_raw_reads","#_unique_reads","%_unique_reads","#_duplicated_reads","%_duplicated_reads","total%_mapping","CDS","UTR","intron","rRNA","TIN(mean)"),"\n";
#		print OUT join("\t",$sample_name,$total,$unique,"$unique_per%",$duplicate,"$duplicate_per%","$mapped_per%","$cds%","$UTR%","$intron%","$result_ribo%",$tin),"\n";
#
#		close OUT;
#		close IN;
#		close INN;
#		close INNN;
#		close INNNN;
#
#		system("cat $sample_name.*.r > $sample_name.all_temp.r");
#
#		open IN, "$sample_name.all_temp.r" or die $!;
#
#		my @script = <IN>;
#
#		foreach my $line (@script){
#	
#			$line =~ s/pdf\(.*//g;
#			$line =~ s/dev\.off\(\)//g;
#		}
#
#		push @script,"dev.off()";
#		unshift @script,"pdf(\"final_results_pdf\")";
#
#
#		open OUT, ">$sample_name.all_new.r" or die $!;
#		print OUT "@script\n";
#
#		close OUT;
#		close IN;
#		
##		my $cmd_pdf_combine='convert -density 250 $(ls -rt '.$sample_name.'*.pdf) final_results_pdf';
#		system("Rscript $sample_name.all_new.r");
#	
##		system($cmd_pdf_combine);
##		my $cmd_txt = "mv $sample_name*.txt ".$prefix."_temp";
##		my $cmd_xls = "mv $sample_name*.xls ".$prefix."_temp";
#		my $cmd_all = "mv $sample_name\.* ".$prefix."_temp";
#		my $cmd_rename1 = "mv final_results_summary ".$sample_name.".final_results_summary.txt";
#		my $cmd_rename2 = "mv final_results_pdf ".$sample_name.".final_results_figures.pdf";
#
##		system($cmd_txt);
##		system($cmd_xls);
#		system($cmd_all);
#		system($cmd_rename1);
#		system($cmd_rename2);
##	$flag = 1;
#	}
#
#	my @final_results = `ls *.final_results_summary.txt`;
#	open OUT, ">$prefix.results_all.txt" or die $!;
#	
#	print OUT join("\t","Sample","#_raw_reads","#_unique_reads","%_unique_reads","#_duplicated_reads","%_duplicated_reads","total%_mapping","CDS","UTR","intron","rRNA","TIN(mean)"),"\n";
#	
#	foreach my $file (@final_results){
#		open IN, $file or die $!;
#		my @data = <IN>;
#		print OUT $data[1];
#		close IN;	
#	}
#				
#	close OUT;
#	system("geneBody_coverage.py -r $hsbed -i $abs_path_allbam -o $prefix");
#	my $cmd_txt = "mv *.geneBodyCoverage.txt ".$prefix."_temp";
#	my $cmd_r = "mv *.geneBodyCoverage.r ".$prefix."_temp";
#	my $cmd_final_results = 'mv *.final_results_summary.txt '.$prefix."_temp";
#	my $cmd_log = 'mv log.txt '.$prefix."_temp";
#	system($cmd_txt);
#	system($cmd_r);
#	system($cmd_final_results);
#	system($cmd_log);	
#}elsif($bam && $bam_dir){
#	die "-i and -ilist can't be used at the same time, please re-specify the input bam files!";
}else{
	die "please specific input bam files!";
}

sub File_Exist {
	my @input = @_;
	foreach my $file (@input){
		chomp $file;
		if(-e $file){
			next;
		}else{
			die "FATAL ERROR: $file doesn't exist!";
			
		}
	}
}
__END__
	








