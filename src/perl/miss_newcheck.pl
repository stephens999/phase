#!/usr/bin/perl

# checks output from postprocess with truth
# usage is
# perl check.pl truth.in guess.in inputfile.in results.out
# -truth file must be of form N R SSSS etc
# -guess file must not have N, R,SSS etc
# inputfile.in includes info on which positions
# were missing

open(TRUTH,"$ARGV[0]");
open(GUESS,"$ARGV[1]");
open(INPUTFILE,"$ARGV[2]");
open(RESULTS,">$ARGV[3]");


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


sub first_nonwhites{
    my($first)=shift(@_);
    $first =~ s/\s*(\S+)//; # first non whitespace characters	
    return($1);
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
	
	$dummy=<INPUTFILE>;
	$dummy=<INPUTFILE>;
	$dummy=<INPUTFILE>;

	#$R=<GUESS>;
	#print "Second line is $R\n";
	#$loci=<GUESS>;
	#print "Third line is $loci\n";
	#remove whitespace
	$loci=remove_whitespace($loci);
	@loci=split(//,$loci);
	$inderrors=0;
	$poserrors=0;
	$ambigind=0; # number of ambiguous individuals
	$ambigpos=0; # number of ambiguous positions
	$unknownpos=0;
	for($n=0;$n<$N;$n=$n+1){
	    $errors0=0; # number of errors counting each of the two possible ways
	    $errors1=0;
	    $hetpos=0; #number of heterozygous positions
	    $missing = 0; # number of missing positions
	    $misshomopos =0;

	    $firsthapt=<TRUTH>;
	    $secondhapt=<TRUTH>;

	    $firsthapi=<INPUTFILE>;
	    $secondhapi=<INPUTFILE>;
	    	    
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
		    $firsthapg=~ s/\s*(\S?)//;
		    $firstalleleg = $1;
		    $secondhapt=~ s/\s*(\S?)//;
		    $secondallelet = $1;	       
		    $secondhapg=~ s/\s*(\S?)//;
		    $secondalleleg = $1;
		    $firsthapi=~ s/\s*(\S?)//;
		    $firstallelei = $1;
		    $secondhapi=~ s/\s*(\S?)//;
		    $secondallelei = $1;  
		} else {
		   $firsthapt =~ s/\s*(\S+)//; # first non whitespace characters	
		    $firstallelet = $1;
		    $firsthapg=~ s/\s*(\S+)//;
		    $firstalleleg = $1;
		    $secondhapt=~ s/\s*(\S+)//;
		    $secondallelet = $1;	       
		    $secondhapg=~ s/\s*(\S+)//;
		    $secondalleleg = $1;
		    $firsthapi=~ s/\s*(\S+)//;
		    $firstallelei = $1;
		    $secondhapi=~ s/\s*(\S+)//;
		    $secondallelei = $1;   
		}
		
		$errors0 += (($firstalleleg ne $firstallelet)  || ($secondalleleg ne $secondallelet));		
		$errors1 += (($firstalleleg ne $secondallelet) || ($secondalleleg ne $firstallelet));		
		$hetpos += ($firstallelet ne $secondallelet);
		if($firstallelei == '-' || $firstallelei == -1){
		    $missing += 1;
		    if($firstallelet == $secondallelet){
			$misshomopos += 1;
		    }
		}
	    }
	    $errors= min($errors0,$errors1);
	    print "Errors: $errors0 $errors1\n";
	    $inderrors += ($errors>0);
	    $poserrors += $errors;
	    use integer;
	    $ambigpos += ($hetpos/2) + $misshomopos; 
	    no integer;
	    $ambigind += (($hetpos+$missing) > 1);
	}
	print RESULTS "$ambigind $inderrors  $ambigpos $poserrors\n";
    }
}

	    
