#! /usr/bin/perl -w

use strict;


my $input = $ARGV[0]; # .BED file (sample file)
chomp $input;
my $val = $ARGV[1]; # p value
chomp $val;

#my $control = $ARGV[1]; # control file in .BED format if present
#chomp $control;



unless ($input =~ m/\.BED$/)
	{print "\n\nInput must be .BED\n\n"; exit;}

my $inputFile="../../".$input;
my @filename1 = split(/\.BED/, $input);


my $name = join ("", @filename1);
print "\n Processing $name\n\n";


my $cmd = "";
my $genome ="mm10";
my $gSize="mm";

my $name2 = $name.".NoModel.NoLambda_p".$val;
$cmd = "macs2 callpeak -t $inputFile  -g $gSize -n $name2 -f BED -p $val --nomodel --nolambda  --bdg ";
print "\n\nRunning command: $cmd\n\n";
system($cmd);

my $new_input = $name2."_peaks.xls";
my $new_output = $name2."_400bp.bed";
$cmd = "peak_finding_result2bed_sv.pl $new_input $new_output 4 $genome " ;
#system($cmd);


#print "\n $name FINISHED. Results in ./MACS2/$folder_tmp3/\n\n";
print "\n $name FINISHED.\n\n";

exit;
