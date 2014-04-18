#!/usr/bin/perl

# checks output with truth and predicted overall error probs
# usage is
# perl newcheck.pl truth.in guess.in guess.errprobs results.out results.err results.corr
# -truth file must be of form N R SSSS etc
# -guess file must not have N,R SSS etc.

open(TRUTH,"$ARGV[0]");
open(GUESS,"$ARGV[1]");
open(ERRPROBS,"$ARGV[2]");
open(RESULTS,">$ARGV[3]");
open(ERRORS,">$ARGV[4]");
open(CORRECT,">$ARGV[5]");


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
	for($n=0;$n<$N;$n=$n+1){
	    $phaseprob = <ERRPROBS>;

	    $phaseerrors0=0;
	    $phaseerrors1=0;
	  
	    $errors0=0; # number of errors counting each of the two possible ways
	    $errors1=0;
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
		$phaseerrors0 += ($firstallelet ne $secondallelet) & ($firstalleleg eq $secondallelet) & ($firstallelet eq $secondalleleg);
		$phaseerrors1 += ($firstallelet ne $secondallelet) & ($firstalleleg eq $firstallelet) & ($secondallelet eq $secondalleleg);

	    }
	    print "phaseerrors0 = $phaseerrors0; phaseerrors1 = $phaseerrors1\n";
	    $phaseerrors= min($phaseerrors0,$phaseerrors1);
	    $indphaseerrors += ($phaseerrors>0);
	    if($phaseerrors>0){
		print ERRORS "$phaseprob";
	    }
	    else{
		print CORRECT "$phaseprob";
	    }
	    $posphaseerrors += $phaseerrors;
	    print "posphaseerrors = $posphaseerrors\n";
	    use integer;
	    $ambigpos += $hetpos/2;
	    $ambigind += ($hetpos > 1);
	}
	print RESULTS "$ambigind $indphaseerrors  $ambigpos $posphaseerrors\n";
    }
}

	    
