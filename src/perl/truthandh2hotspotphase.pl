#!/usr/bin/perl

# converts output from hudson program, plus truth file,
# into format suitable for input into PHASE-R
# Merges pairs of consecutive data sets to make an
# "infinite hotspot" in the middle. Flags that
# hotspot by a large change in positions - not necessary
# usage is
# perl truthandh2hotspotphase.pl truth.in hudson.in phase-R.inp numberofchromoosmes

open(TRUTH,"$ARGV[0]");
open(HUDSON,"$ARGV[1]");
open(OUT,">$ARGV[2]");

$nchr = $ARGV[3];

while(<TRUTH>){
    $line1[0] =$_;
    chomp($line1[0]);

    for($i=1; $i<$nchr; $i=$i+1){
	$line1[$i]=<TRUTH>;
	chomp($line1[$i]);
    }
    for($i=0; $i<$nchr; $i=$i+1){
	$line2[$i]=<TRUTH>;
    }
    #the -1 in the following is because line2 has an endl
    $nloci = length($line1[0])+length($line2[0])-1;
    
    $N=$nchr/2;
    print OUT "$N\n";
    print OUT "$nloci\n";
    print OUT "P ";
    
    #print out positions of first data set
    $dummy = <HUDSON>;
    until($dummy =~ m/segsites/){
	$dummy=<HUDSON>;
    }
    $dummy = <HUDSON>;
    while($dummy =~ m/\./){
	print OUT "$dummy";
	$dummy=<HUDSON>;
    }
    # now positions of second data set
    until($dummy =~ m/segsites/){
	$dummy=<HUDSON>;
    }
    $dummy = <HUDSON>;
    while($dummy =~ m/\./){
	
	@pos = split(/ /,$dummy);
	
	for($i=0; $i<($#pos); $i=$i+1){
	    $pos[$i] = $pos[$i]+1000;
	    print OUT "$pos[$i] ";
	}
#	print OUT "$dummy";
	$dummy=<HUDSON>;
    }
    
    print OUT "\n";
    print OUT 'S' x $nloci."\n";
    
    for($i=0; $i<$nchr; $i=$i+1){
	print OUT "$line1[$i]$line2[$i]";
    }    

}
    
