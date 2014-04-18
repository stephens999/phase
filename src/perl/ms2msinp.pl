#!/usr/bin/perl

# converts output from microsat program into format
# suitable for input into PHASE
# usage is
# perl ms2msinp.pl ms.in ms.inp

open(IN,"$ARGV[0]");
open(OUT,">$ARGV[1]");

$firstline=<IN>; # first line is 2*$N
$N=$firstline;

while(<IN>){
    $colonR=$_; # next line is :R
    $colonR =~ m/:(\d+)/;
    $R=$1;
    print OUT "$N";
    
    $dummy = <IN>; # dummy line of positions
   
    #and add further haplotypes
    for($i=0;$i<(2*$N);$i+=1){
	$temp=<IN>;
	push(@hap,$temp);
    }
    
    print OUT $R."\n";
    print OUT 'M' x $R."\n";

    # now output haplotypes in a random order
    # (as per Perl for Dummies, p 304)

    for($i=0;$i<(2*$N);$i+=1){
	print OUT splice(@hap,rand(@hap),1);
    }
    print OUT "\n";
}

