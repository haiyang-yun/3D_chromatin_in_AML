#! /usr/bin/perl -w

##############################################################################
# Script to QC the raw data, aligned the paired reads, convert the aligned 
# reads in SAM format to bedgraph and bigwig and finally output HT-seq counts
# Step 1: Process each fastq file (fastqc) and collect joint list of adapters
# Step 2: Run STAR on trimmed files
# Step 3: Convert sam to bigwig 
# Step  4: Sort uniquely mappable sam file and generate raw counts
##############################################################################

# Example usage of running this script: runRNA_STAR_paired.pl [INPUT_r_1_FASTQ] [INPUT_r_2_FASTQ] [GENOTYPE] mm10 STAR-GENOMES-mm10.gencode.vM7.comprehensive gencode.vM7.comprehensive.annotation.gtf [exons y/n]
# [INPUT_r_1_FASTQ] [INPUT_r_2_FASTQ]: input paired reads
# [GENOTYPE]: folder named as genotype (e.g. WT, FLT3, NPM1, DM) where fastq files are stored
# STAR-GENOMES-mm10.gencode.vM7.comprehensive: path to STAR genome folder
# gencode.vM7.comprehensive.annotation.gtf: annoataion file

use strict;
use File::Basename;

if (@ARGV != 7){ print("\n
Please provide: 
<fastq file 1> 
<fastq file 2>
<sample name> 
<genome mm10/hg19> 
<path to STAR genome folder >
<annoataion File>
<exons y/n>
\n\n"); exit;}

my @paired_files=();

my $file1 = $ARGV[0];
chomp $file1;

push(@paired_files, $file1);

my $file2 = $ARGV[1];
chomp $file2;

push(@paired_files, $file2);

my $file1_basename = basename($file1,".fastq");
my $file2_basename = basename($file2,".fastq");

my $GSM = $ARGV[2];
chomp $GSM;

my $genome = $ARGV[3];
chomp $genome;

my $genomeDir="/serenity/data/reference-genomes/".$ARGV[4];
chomp $genomeDir;

my $annotationFile = $genomeDir."/".$ARGV[5];
chomp $annotationFile;

my $exons = $ARGV[6];
chomp $exons;

#unless ($genome eq "hg19"){print "\nCurrently only set up to run hg19... Exiting.\n"; exit;}
#unless ($genome eq "mm10" || $genome eq "hg19"){print "\n\n'Genome not recognised. Please enter 'mm10' or 'hg19' \n"; exit;}
unless ($exons eq "y" || $exons eq "n"){print "\n\n'exons_only' choice not recognised. Please enter only either 'y' or 'n' (lowercase).\n"; exit;}


my $total_overrep = 0;

## Step 1: Process each fastq file (fastqc) and collect joint list of adapters
foreach my $file (@paired_files)
{
### run fastqc on file1 and file2
print "\nProcessing $file...\n";
print "\nRunning fastqc...\n\n";
my $cmd = "fastqc --quiet -f fastq $file";
system($cmd);
### check if overrepresented sequences found
my $basename = basename($file,".fastq");
my $fastqc_data = "./".$basename."_fastqc/fastqc_data.txt";
my $overrep = `grep \">>Overrepresented sequences\" $fastqc_data`;
$cmd= "grep \">>Overrepresented sequences\" $fastqc_data > grep_test.txt";
system($cmd);
if ($overrep =~ m/fail/)
{
	$total_overrep++;
	### retrieve overrep seqs from report
	my $adapters_out = $basename."_adapters.txt";
	my $cmd= "grep -i \"adapter\" $fastqc_data > $adapters_out";
	system($cmd);
	
	##retrieve illumina quality encoding info from report
	my $illumina_out = $basename."_illumina.txt";
	$cmd= "grep -i \"Encoding\" $fastqc_data > $illumina_out";
	system($cmd);
} # close for if 
}#close foreach

if ($total_overrep > 0)
{
### Combine adapters from both files, sort and remove duplicates
my $cmd = "cat *_adapters.txt | sort | uniq > all_adapters.txt";
system($cmd);
### Combine illumina encoding, sort, remove duplicates, extract score and set Phred
my $illumina = `cat *_illumina.txt | sort | uniq`;
chomp $illumina;
my $phred = 0;
my $score = 0;
if ($illumina =~ m/\s+(\d+\.\d+)\s*$/)
{
	$score = $1;
	chomp $score;
} # close for if
if ($score == 1.5)
{
	$phred = "phred64";
}
elsif ($score == 1.9)
{
	$phred = "phred33";
}
else
{
print "\n\nIllumina quality encoding not recognised!\n\n"; 
exit;
} # close for if

### check if overrep seqs contains adapters specifically
open (TEMP, "all_adapters.txt");
my @adapters = <TEMP>;
close(TEMP);
my $adaps_found = 0;
foreach my $adapter_line (@adapters)
{
	if ($adapter_line =~ m/adapter/i) { $adaps_found = 1; last;}
}
### if adapters found, iterate trimming for all adapters found
if ($adaps_found == 1)
{
print "\nAdapters found...\n";
### make copies of the original untrimmed fastq files before trimming
my $cmd= "mkdir originals";
system($cmd);
$cmd= "cp $file1 ./originals/";
system($cmd);
$cmd= "cp $file2 ./originals/";
system($cmd);
	my $j = 0;	
	foreach my $adapter_line (@adapters)
	{
		$j++;
		print "\nTrimming: $adapter_line\n\n";
		my $length = 0;
		my $adapter = "";
		if ($adapter_line =~ m/^(\w+)\t\d+/)
		{
			$adapter = $1;
			chomp $adapter;
			$length = length($adapter);
		}

	my $cmd= "trim_galore --paired --$phred --quality 20 -a $adapter --stringency $length $file1 $file2";
	print "\n$cmd\n\n";
	system($cmd);
	
	my $trimmed1 = $file1_basename."_val_1.fastq";
	$cmd= "mv $trimmed1 $file1";
	system($cmd);
	
	my $trimmed2 = $file2_basename."_val_2.fastq";
	$cmd= "mv $trimmed2 $file2";
	system($cmd);
	
	my $trimmingLog = $file1."_trimming_report.txt";
	my $newTrimmingLog = $file1."_trimming_report".$j.".txt";
	$cmd= "mv $trimmingLog $newTrimmingLog";
	system($cmd);
	
	$trimmingLog = $file2."_trimming_report.txt";
	$newTrimmingLog = $file2."_trimming_report".$j.".txt";
	$cmd= "mv $trimmingLog $newTrimmingLog";
	system($cmd);
	
	}#close foreach adapater
}#close if adaps found
}#close if overrep seqs found ($total_overrep)


## Step 2: Run STAR and extract uniquely mappable reads 
my $cmd= "mkdir STAR_out";
system($cmd);
chdir("STAR_out");

### Run STAR
$cmd= "/serenity/software/STAR_2.4.0k/bin/Linux_x86_64/STAR --genomeDir $genomeDir --readFilesIn ../$file1 ../$file2 --runThreadN 5 --outFilterIntronMotifs None --outFilterMismatchNmax 2 --outFileNamePrefix $GSM.";
system($cmd);

### extract only uniquely mappable reads (NH:i:1 and MAPQ=255) and sam header, remove flag 'nM:i:n' which is unrecognised by sam2bed and remove white space from the ends of each line:
print "\n\nExtracting uniquely mappable reads from sam...\n";
$cmd= "grep \'NH:i:1\\s*\\|NH:i:1\$\\|^\@\' $GSM.Aligned.out.sam | awk \'{if(\$5 == 255) print \$0; else if (\$1 ~ /^\@/) print \$0;}\' | sed \'s/nM:i:[0-9]*//g\' | sed \'s/\\s*\$//\' > $GSM.uniquely_mappable.sam";
system($cmd);

## Step  3: Generate raw counts using HTSeq
### run HTseq
print "\n\nRunning HTSeq...\n";
$cmd= "python -m HTSeq.scripts.count --stranded no --type exon -m intersection-nonempty $GSM.uniquely_mappable.sam $annotationFile --samout $GSM.htseq.samout > $GSM.htseq.counts";
#$cmd= "/usr/bin/python -m HTSeq.scripts.count --stranded no --type exon -m intersection-nonempty $GSM.uniquely_mappable.sam $annotationFile > $GSM.htseq.counts";
system($cmd);
print "\n\nhtseq.counts ready for running DE-Seq! :-)\n\n";

exit;
