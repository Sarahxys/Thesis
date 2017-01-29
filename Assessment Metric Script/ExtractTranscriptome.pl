#!/usr/bin/perl

use strict;
use warnings;

my $inputfile1 =$ARGV[0];
my $inputfile2 =$ARGV[1];
my $outputfile = "Extracted$inputfile2";
my %name_hash;
my $pass = 0;

open(INPUT,"<", "$inputfile1") or die "could not open the input file1";

while (my $line = <INPUT>){
	chomp($line);
	if ($line =~ "Summary"){
		last;
	}
	else{ $name_hash{$line} = 1;}
}
close INPUT;

open(INPUT,"<", "$inputfile2") or die "could not open the input file2";
open (OUTPUT, ">", $outputfile) or die "the output file can not be created"; #opent the first input file
my @line;
while (my $line = <INPUT>){
	chomp($line);
	if ($line !~ /^>/ && $pass == 1){
		print OUTPUT $line,"\n";	
	}
	elsif($line =~ /^>/){
		@line = split (/\s+/,$line); #split the line with spaces (\s = space; and \s+ = more than one spece);
		my $seqid = $line[0]; 
		$pass = 0;
		$seqid =~ s/^>//g;#extract the seqid from the line;
		if (exists $name_hash{$seqid} && $name_hash{$seqid} == 1){
			print OUTPUT $line,"\n";
			$pass = 1;
		}

	}	
}
close INPUT;
close OUTPUT;
exit;



