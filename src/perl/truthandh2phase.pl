#!/usr/bin/perl

# converts output from hudson program, plus truth file,
# into format suitable for input into PHASE-R
# usage is
# perl truthandhh2phase.pl truth.in hudson.in phase-R.inp

open(TRUTH,"$ARGV[0]");
open(HUDSON,"$ARGV[1]");
open(OUT,">$ARGV[2]");

$nchr = $ARGV[3];

while(<TRUTH>){
    $firstline=$_;
    #$firstline=<HUDSON>;
    #$secondline=<HUDSON>;
    #@params = split(/\s+/,$secondline);
    #$N=$params[1]/2;

    $nloci=length($firstline)-1;    
    $N=$nchr/2;
    print OUT "$N\n";
    print OUT "$nloci\n";
    print OUT "P ";
    
    #print out positions!
    $dummy = <HUDSON>;
    until($dummy =~ m/segsites/){
	$dummy=<HUDSON>;
    }
    $dummy = <HUDSON>;
    while($dummy =~ m/\./){
	print OUT "$dummy";
	$dummy=<HUDSON>;
    }

    print OUT 'S' x $nloci."\n";
    print OUT "$firstline";
    
    for($i=1; $i<$nchr; $i=$i+1){
	$line=<TRUTH>;
	print OUT $line;
    }
}
    
