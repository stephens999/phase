#!/usr/bin/perl

# counts number of genotype errors compared with truth
# usage is
# perl countgenotypeerrors.pl phase.inp guess.inp results.out
# -truth file must 
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


while(<TRUTH>){
    $N=$_;
    if(($N=~/\S/)){ #if not just blanks	
	$R=<TRUTH>;
	$P=<TRUTH>;
	while(!($P=~/[MS]/)){
	    $P = <TRUTH>;
	}
	$loci=$P;
	$loci=remove_whitespace($loci);
	@loci=split(//,$loci);

	$genotypeerrors = 0;
	$inderrors=0;
	$poserrors=0;
	$indphaseerrors = 0;
	$posphaseerrors = 0;
	$ambigind=0; # number of ambiguous individuals
	$ambigpos=0; # number of ambiguous positions
	$totalhetpos = 0;

	for($n=0;$n<$N;$n=$n+1){
	    $phaseerrors0=0;
	    $phaseerrors1=0;
	  
	    $errors0=0; # number of errors counting each of the two possible ways
	    $errors1=0;
	    $hetpos=0; #number of heterozygous positions

	    $firsthapt=<TRUTH>;
	    $secondhapt=<TRUTH>;
	    
	    $firsthapg=<GUESS>;
	    $isblank = ($firsthapg=~/\S/);
	   
	    while(!($firsthapg=~/\S/)){ #while just blanks
		$firsthapg = <GUESS>;
	    }

	    $secondhapg = <GUESS>;

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
		
		$e0 = ($firstalleleg ne $firstallelet) + ($secondalleleg ne $secondallelet);
		$e1 = ($firstalleleg ne $secondallelet) + ($secondalleleg ne $firstallelet);
		$errors0 += $e0;	
		$errors1 += $e1;
		$genotypeerrors += min($e0,$e1);
		$hetpos += ($firstallelet ne $secondallelet);
		$phaseerrors0 += ($firstallelet ne $secondallelet) & ($firstalleleg eq $secondallelet) & ($firstallelet eq $secondalleleg);
		$phaseerrors1 += ($firstallelet ne $secondallelet) & ($firstalleleg eq $firstallelet) & ($secondallelet eq $secondalleleg);
	    }
	    print "phaseerrors0 = $phaseerrors0; phaseerrors1 = $phaseerrors1\n";
	    $phaseerrors= min($phaseerrors0,$phaseerrors1);
	    $indphaseerrors += ($phaseerrors>0);
	    $posphaseerrors += $phaseerrors;
	    print "posphaseerrors = $posphaseerrors\n";
	    use integer;
	    $totalhetpos += $hetpos;
	    $ambigpos += $hetpos/2;
	    $ambigind += ($hetpos > 1);
	}
	print RESULTS "$ambigind $indphaseerrors  $ambigpos $posphaseerrors $totalhetpos $genotypeerrors\n";
    }
}

	    
