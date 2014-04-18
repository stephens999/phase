#!/usr/bin/perl

# converts output from PHASE applied to infinite "hotspot"
#data,  created by merging consecutive data sets in the truth
# file, into format suitable for comparison with the truth file
# usage is
# perl unmergehotspotouts.pl phase.out truth.in phase.unmerged.out numberofchromoosmes

open(MERGED,"$ARGV[0]");
open(TRUTH,"$ARGV[1]");
open(UNMERGED,">$ARGV[2]");

$nchr = $ARGV[3];

while(<TRUTH>){
    $line1[0] =$_;

    for($i=1; $i<$nchr; $i=$i+1){
	$line1[$i]=<TRUTH>;
    }
    for($i=0; $i<$nchr; $i=$i+1){
	$line2[$i]=<TRUTH>;
    }

    $nloci1 = length($line1[0])-1;
    $nloci2 = length($line2[0])-1;
    
    for($i=0; $i<$nchr; $i=$i+1){
	$line[$i] = <MERGED>;
	print $line[$i];
	@chrom = split(/ /,$line[$i]);
	print join(' ',@chrom); 
	$temp = join(' ',splice(@chrom, 0,$nloci1, ()));
	print UNMERGED "$temp\n";
	$line[$i] = join(' ',@chrom);
    }
    
    for($i=0; $i<$nchr; $i=$i+1){
	print UNMERGED join(' ',$line[$i]);
    }
    print UNMERGED "\n";
    $dummy = <MERGED>; # read blank line


}
    
