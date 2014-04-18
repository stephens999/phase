#!/usr/bin/perl

# checks output from postprocess with truth
# usage is
# perl check.pl truth.in guess.in results.out

open(TRUTH,"$ARGV[0]");
open(GUESS,"$ARGV[1]");
open(RESULTS,">$ARGV[2]");


sub remove_whitespace{
    my($string)= @_;
    chomp($string);
    $string =~ s/\s+//g;
    return($string);
}

sub min{
    my($first)= shift(@_);
    my($second)= shift(@_);
    if($first>$second){
	return $second;
    } else {
	return $first;
    }
}


while(<GUESS>){
    print "First char is $_\n";
    $N=<TRUTH>;
    while(!($N=~/\S/)){ #while N is just blanks
	$N=<TRUTH>;}
    $R=<TRUTH>;
    $loci=<TRUTH>;
    $R=<GUESS>;
    print "Second line is $R\n";
    $loci=<GUESS>;
    print "Third line is $loci\n";
    #remove whitespace
    $loci=remove_whitespace($loci);
    @loci=split(//,$loci);
    $inderrors=0;
    $poserrors=0;
    $ambigind=0; # number of ambiguous individuals
    $ambigpos=0; # number of ambiguous positions
    for($n=0;$n<$N;$n=$n+1){
	$errors0=0; # number of errors counting each of the two possible ways
	$errors1=0;
	$hetpos=0; #number of heterozygous positions
	$firsthapt=<TRUTH>;
	$secondhapt=<TRUTH>;
	$firsthapg=<GUESS>;
	$secondhapg=<GUESS>;
	
	foreach $i (@loci){
	    if($i eq 'S'){	
		$firsthapt =~ s/\s*(\S?)//; # first non whitespace character	
		$firstallelet = $1;
		
		$firsthapg =~ s/\s*(\S?)//; # first non whitespace character
		$firstalleleg = $1;
		
		$secondhapt =~ s/\s*(\S?)//; # first non whitespace character
		$secondallelet = $1;
		
		$secondhapg =~ s/\s*(\S?)//; # first non whitespace character
		$secondalleleg = $1;
		
	    } else {
		$firsthapt =~ s/\s*(\S+)//; # first non whitespace characters	
		$firstallelet = $1;

		$firsthapg =~ s/\s*(\S+)//; # first non whitespace characters
		$firstalleleg = $1;		

		$secondhapt =~ s/\s*(\S+)//; # first non whitespace characters
		$secondallelet = $1;
		
		$secondhapg =~ s/\s*(\S+)//; # first non whitespace characters
		$secondalleleg = $1;
		
	    }
	    $errors0 += ($firstalleleg ne $firstallelet);
	    $errors0 += ($secondalleleg ne $secondallelet);
	    $errors1 += ($firstalleleg ne $secondallelet);
	    $errors1 += ($secondalleleg ne $firstallelet);
	    $hetpos += ($firstallelet ne $secondallelet);
	}
	$errors= min($errors0,$errors1);
	$inderrors += ($errors>0);
	$poserrors += $errors;
	$ambigpos += $hetpos/2;
	$ambigind += ($hetpos > 1);
    }
    print RESULTS "$ambigind $inderrors  $ambigpos $poserrors\n";
}

	    
