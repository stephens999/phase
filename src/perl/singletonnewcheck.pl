#!/usr/bin/perl

# checks output with truth
# usage is
# perl newcheck.pl truth.in guess.in results.out
# -truth file must be of form N R SSSS etc
# -guess file must not have N,R SSS etc.

open(TRUTH,"$ARGV[0]");
open(GUESS,"$ARGV[1]");
open(RESULTS,">$ARGV[2]");


sub remove_whitespace{
    my($string)= @_;
    chomp($string);
    $string =~ s/\s+//g;
    return($string);
}

sub remove_brackets{
    my($string)= @_;
    $string =~ s/\(//g;
    $string =~ s/\)//g;	 
    $string =~ s/\[//g;
    $string =~ s/\]//g;	 
 
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
    $firsthapg=$_;
    if(($firsthapg=~/\S/)){ #if not just blanks
	print "$firsthapg\n";
	$N=<TRUTH>;
	print "First char is $N\n";
	while(!($N=~/\S/)){ #while N is just blanks
	    $N=<TRUTH>;}
	$R=<TRUTH>;
       
	
	$loci=<TRUTH>;
	until($loci=~/[SM]/){
	    $loci=<TRUTH>; 
	}
#$R=<GUESS>;
	#print "Second line is $R\n";
	#$loci=<GUESS>;
	#print "Third line is $loci\n";
	#remove whitespace
	$loci=remove_whitespace($loci);
	@loci=split(//,$loci);
	$inderrors=0;
	$poserrors=0;
	$indphaseerrors = 0;
	$posphaseerrors = 0;
	$ambigind=0; # number of ambiguous individuals
	$ambigpos=0; # number of ambiguous positions

	$locus = 0;
	foreach $i (@loci){
	    $freq[$locus] =0;
	    $locus = $locus +1;
	}
	
	for($n=0;$n<$N;$n=$n+1){
	    $locus = 0;
	    foreach $i (@loci){
		$phaseerrors0[$locus]=0;
		$phaseerrors1[$locus]=0;	  
		$errors0[$locus]=0; # number of errors counting each of the two possible ways
		$errors1[$locus]=0;
		
		$phaseerrorsarray[$locus] =0;
		$locus = $locus +1;
	    }
	    $sumphaseerrors0 = 0;
	    $sumphaseerrors1 =0;

	    $hetpos=0; #number of heterozygous positions
	    $firsthapt=<TRUTH>;
	    $secondhapt=<TRUTH>;
	    if($n>0){
		$firsthapg=<GUESS>;	    
	    }
	    $secondhapg=<GUESS>;

	    $firsthapg=remove_brackets($firsthapg);
	    $secondhapg=remove_brackets($secondhapg);

	    print "n= $n\n";
	    print "First hap is $firsthapt\n";
	    print "First hap is $firsthapg\n";
	    print "2nd hap is $secondhapt\n";
	    print "2nd hap is $secondhapg\n";
	    
	    $locus = 0;
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

		    $freq[$locus] += $firstallelet + $secondallelet;

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
		$errors0[$locus] += ($firstalleleg ne $firstallelet);
		$errors0[$locus] += ($secondalleleg ne $secondallelet);
		$errors1[$locus] += ($firstalleleg ne $secondallelet);
		$errors1[$locus] += ($secondalleleg ne $firstallelet);
		$hetpos += ($firstallelet ne $secondallelet);
		$phaseerrors0[$i] += ($firstallelet ne $secondallelet) & ($firstalleleg eq $secondallelet) & ($firstallelet eq $secondalleleg);
		$phaseerrors1[$i] += ($firstallelet ne $secondallelet) & ($firstalleleg eq $firstallelet) & ($secondallelet eq $secondalleleg);
		$sumphaseerrors0 += ($firstallelet ne $secondallelet) & ($firstalleleg eq $secondallelet) & ($firstallelet eq $secondalleleg);
		$sumphaseerrors1 += ($firstallelet ne $secondallelet) & ($firstalleleg eq $firstallelet) & ($secondallelet eq $secondalleleg);
		$locus = $locus +1;

	    }
	    
	    print "sumphaseerrors0 = $sumphaseerrors0; sumphaseerrors1 = $sumphaseerrors1\n";
	    $locus =0;
	    foreach $i (@loci){
		if($sumphaseerrors0<$sumphaseerrors1){
		    $phaseerrorsarray[$locus] += $phaseerrors0[$locus];
		} else {
		    $phaseerrorsarray[$locus] += $phaseerrors1[$locus];
		}
		$locus = $locus +1;
	    }

	    $phaseerrors= min($sumphaseerrors0,$sumphaseerrors1);
	    $indphaseerrors += ($phaseerrors>0);
	    $posphaseerrors += $phaseerrors;
	    print "posphaseerrors = $posphaseerrors\n";
	    $ambigpos += $hetpos/2;
	    $ambigind += ($hetpos > 1);
	}
	
	$singletonerrors = 0;
	$locus = 0;
	foreach $i (@loci){
	    if(($freq[$locus]==1) or ($freq[$locus]==(2*$N-1))){
		$singletonerrors += $phaseerrorsarray[$locus];
	    }
	    $locus = $locus +1;
	}

	print RESULTS "$ambigind $indphaseerrors  $ambigpos $posphaseerrors $singletonerrors \n";
    }
}

	    
