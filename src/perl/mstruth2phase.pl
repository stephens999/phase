#!/usr/bin/perl

# converts "truth" files from ms program into format
# suitable for input into phase
# usage is
# perl mstruth2phase.pl truth.in hudson.inp nchr
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
    print OUT "P 0.0909090909090909 0.181818181818182 0.272727272727273 0.363636363636364 0.454545454545455 0.545454545454545 0.636363636363636 0.727272727272727 0.818181818181818 0.909090909090909\n";
    print OUT 'M' x $nloci."\n";
    
    @splitfirstline = split(/\s*/,$firstline);
    for($i=0;$i<$nloci;$i=$i+1){
	print OUT ord("$splitfirstline[$i]").' ';
    }
    print OUT "\n";

    for($n=1; $n<$nchr; $n=$n+1){
	$line=<IN>;
	@splitline = split(/\s*/,$line);
	for($i=0;$i<$nloci;$i=$i+1){
	    print OUT ord("$splitline[$i]").' ';
	}
	print OUT "\n";
    }
}
    

