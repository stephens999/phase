#!/usr/bin/perl

# adds missing alleles to PHASE input file
# note that input file must have a position line
# with all positions on the same line
#
# usage is
# perl addmiss.pl hudson.inp hudson.miss.inp missing_proportion

sub remove_whitespace{
    my($string)= @_;
    chomp($string);
    $string =~ s/\s+//g;
    return($string);
}


open(IN,"$ARGV[0]");
open(OUT,">$ARGV[1]");
$proportion = $ARGV[2];


while(<IN>){
    $N=$_;
    if(($N=~/\S/)){ #if not just blanks	
	print OUT $N;

	$R=<IN>;
	print OUT $R;

	$P=<IN>;
	while(!($P=~/[MS]/)){
	    print OUT $P;
	    $P = <IN>;
	}

	$loci=$P;
	print OUT $loci;
    
	$loci=remove_whitespace($loci);
	@loci=split(//,$loci);
	
	for($n=0;$n<$N;$n=$n+1){
	    
	    $firsthap=<IN>;
	    $secondhap=<IN>;
	    
	    $newfirsthap = '';
	    $newsecondhap = '';

	    foreach $i (@loci){
		if($i eq 'S'){	
		    		    
		    $firsthap =~ s/\s*(\S?)//; # first non whitespace character	
		    $firstallele = $1;
		    		    		    
		    $secondhap =~ s/\s*(\S?)//; # first non whitespace character
		    $secondallele = $1;
		    
		    if(rand()<$proportion){
			$firstallele = '?';
			$secondallele = '?';
		    }

		    
		} else {
		    $firsthap =~ s/\s*(\S+)//; # first non whitespace character	
		    $firstallele = $1;
		    		    		    
		    $secondhap =~ s/\s*(\S+)//; # first non whitespace character
		    $secondallele = $1;
		    
		    if(rand()<$proportion){
			$firstallele = -1;
			$secondallele = -1;
		    }
		}
	       
		$newfirsthap = $newfirsthap.$firstallele;
		$newsecondhap = $newsecondhap.$secondallele;
	    }
	    print OUT "$newfirsthap\n";
	    print OUT "$newsecondhap\n";
		
	}	
    }
}

	    
