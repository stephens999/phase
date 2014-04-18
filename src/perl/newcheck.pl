#!/usr/bin/perl

# checks output with truth
# usage is
# perl newcheck.pl truth.in guess.in results.out
# -truth file must be of form N R SSSS etc (may have P line)
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
        $indphaseerrorsnonmissing = 0;

	$posphaseerrors = 0;
	$genotypeerrors = 0;
	$ambigind=0; # number of ambiguous individuals
	$ambigpos=0; # number of ambiguous positions
	$totalhetpos = 0;
        $switchtot = 0;
	$sumswitchaccuracy = 0;

	$ambigindnonmissing =0; # number of ambiguous inds, looking only at pos with no missing alleles

	for($n=0;$n<$N;$n=$n+1){
	    $phaseerrors0=0;
	    $phaseerrors1=0;
	    $phaseerrors0nonmissing = 0;
	    $phaseerrors1nonmissing = 0;
	    $lastphaseerr0 = -1;
	    $ns =0;
	    $errors0=0; # number of errors counting each of the two possible ways
	    $errors1=0;
	    $hetpos=0; #number of heterozygous positions
	    $nonmissinghets  =0; #number of non-missing het positions
	    
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

		$nmissing = ($firstallelet eq '?') + ($secondallelet eq '?');
		
		$e0 = (($firstalleleg ne $firstallelet) & ($firstallelet ne '?')) + (($secondalleleg ne $secondallelet) & ($secondallelet ne '?')) ;
		$e1 = (($firstalleleg ne $secondallelet) & ($secondallelet ne '?')) + (($secondalleleg ne $firstallelet) & ($firstallelet ne '?'));
		$errors0 += $e0;
		$errors1 += $e1;

		if(($e0>0) & ($e1>0)){
		    print "ERROR! guess is inconsistent with input file!\n";
		    exit;
		}

		$genotypeerrors += min($e0,$e1);
		$hetpos += ($firstallelet ne $secondallelet);
		$nonmissinghets += ($nmissing == 0) * ($firstallelet ne $secondallelet);
		#print "firstallelet = $firstallelet\n";
		#print "firstalleleg = $firstalleleg\n";
		#print "secondallelet = $secondallelet\n";
		#print "secondalleleg = $secondalleleg\n";
		

		$phaseerr1 = ($firstallelet ne $secondallelet) & ((($firstalleleg ne $secondallelet) & ($secondallelet ne '?')) | (($secondalleleg ne $firstallelet) & ($firstallelet ne '?')));
		$phaseerr0 = ($firstallelet ne $secondallelet) & ((($firstalleleg ne $firstallelet) & ($firstallelet ne '?')) | (($secondalleleg ne $secondallelet) & ($secondallelet ne '?')));

		if(($firstallelet ne $secondallelet) & ($nmissing == 0)){
			if(($phaseerr0 ne $lastphaseerr0) &($lastphaseerr0 ne (-1))){
			$ns = $ns + 1; 
			}
		    $lastphaseerr0 = $phaseerr0;
		}

		$phaseerrors0 += $phaseerr0;	
          	$phaseerrors1 += $phaseerr1;
		$phaseerrors0nonmissing += $phaseerr0 * ($nmissing == 0);
		$phaseerrors1nonmissing += $phaseerr1 * ($nmissing == 0);

#		print "$phaseerrors0\n";
#		print "$phaseerrors1\n";
#		print "$ns\n";

	    }
	    print "phaseerrors0 = $phaseerrors0; phaseerrors1 = $phaseerrors1\n";
		$phaseerrors= min($phaseerrors0,$phaseerrors1);
		$phaseerrorsnonmissing = min($phaseerrors0nonmissing, $phaseerrors1nonmissing);  
        	$indphaseerrors += ($phaseerrors>0);
                $indphaseerrorsnonmissing += ($phaseerrorsnonmissing>0);

	   	$posphaseerrors += $phaseerrors;
	    	print "posphaseerrors = $posphaseerrors\n";
	    	if($nonmissinghets >1){
			$sumswitchaccuracy = $sumswitchaccuracy + ($nonmissinghets-1-$ns);
	    	        $switchtot = $switchtot + ($nonmissinghets-1);
                }
	    	use integer;
	    	$totalhetpos += $hetpos;
	    	$ambigpos += $hetpos/2;
	    	$ambigind += ($hetpos > 1);
	    	$ambigindnonmissing += ($nonmissinghets > 1);
	
}
	$averageswitch = $sumswitchaccuracy / $switchtot;
	print RESULTS "$ambigind $indphaseerrors  $ambigpos $posphaseerrors $totalhetpos $genotypeerrors $sumswitchaccuracy $switchtot $averageswitch $ambigindnonmissing $indphaseerrorsnonmissing $N";
    }
}

	    
