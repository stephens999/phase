#!/usr/bin/perl

# converts "truth" files from hudson program into format
# suitable for input into phase
# usage is
# perl hudsontruth2phase.pl truth.in hudson.inp nchr
# where nchr is the number of chromosomes in each data set

open(IN,"$ARGV[0]");
open(OUT,">$ARGV[1]");
$nchr=$ARGV[2];
$N=$nchr/2;

while(<IN>){
    $firstline=$_;
    $nloci=length($firstline)-1;    
    print OUT "$N\n";
    print OUT "$nloci\n";
    print OUT 'S' x $nloci."\n";
    print OUT "$firstline";
    
    for($i=1; $i<$nchr; $i=$i+1){
	$line=<IN>;
	print OUT $line;
    }
}
    

