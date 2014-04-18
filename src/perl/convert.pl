#!/usr/bin/perl

# converts output from hudson program into format
# suitable for input into PHASE
# usage is
# perl convert.pl hudson.in hudson.inp

open(IN,"$ARGV[0]");
open(OUT,">$ARGV[1]");

$firstline=<IN>;
$secondline=<IN>;
@params = split(/\s+/,$secondline);
$N=$params[1]/2;

while(<IN>){
    print OUT "$N\n";
    $dummy=<IN>;

    # find first line of positions
    until($dummy =~ m/\./){
	$dummy=<IN>;
    }
    #cycle through lines of positions
    while($dummy =~ m/\./){
	$dummy=<IN>;
    }
    #strip out white space
    chomp($dummy);
    $dummy =~ s/\s+//g;
    $R = length($dummy);

    #now put first haplotype in list
    push(@hap,"$dummy\n");
    
    #and add further haplotypes
    for($i=1;$i<(2*$N);$i+=1){
	$temp=<IN>;
	push(@hap,$temp);
    }
    
    print OUT $R."\n";
    print OUT 'S' x $R."\n";

    # now output haplotypes in a random order
    # (as per Perl for Dummies, p 304)

    for($i=0;$i<(2*$N);$i+=1){
	print OUT splice(@hap,rand(@hap),1);
    }
    print OUT "\n";
}

