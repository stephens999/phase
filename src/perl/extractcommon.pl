#!/usr/bin/perl

# adds missing alleles to PHASE input file
# note that input file must have a position line
# with all positions on the same line
#
# usage is
# perl extractcommon.pl hudson.inp hudson.common.inp cutofffreq 

sub remove_whitespace{
    my($string)= @_;
    chomp($string);
    $string =~ s/\s+//g;
    return($string);
}


open(IN,"$ARGV[0]");
open(OUT,">$ARGV[1]");
$cutoff = $ARGV[2];


while(<IN>){
    $N=$_;
    if(($N=~/\S/)){ #if not just blanks	
	$R=<IN>;
	$loci=<IN>;
	$P = '';
	while(!($loci=~/[MS]/)){
	    $P = $P.$loci; 
		$loci = <IN>;
	}
	    
	$loci=remove_whitespace($loci);
	@loci=split(//,$loci);
	@firsthap = ();
	@secondhap = ();	
	for($n=0;$n<$N;$n=$n+1){
	    
	    $firsthap[$n]=<IN>;
	    $secondhap[$n]=<IN>;
	}
 	$j=0;
	@uselocus = ();
	foreach $i (@loci){
		$uselocus[$j] =1;
		$totalallele =0;
		%count = ();
		if($i eq 'S'){	
			for($n=0;$n<$N;$n=$n+1){

				$fh = remove_whitespace($firsthap[$n]);
				$sh = remove_whitespace($secondhap[$n]);			    		    
		    		@fha = split(//,$fh);
				@sha = split(//,$sh);
				$firstallele = $fha[$j];
		    		$secondallele = $sha[$j];    		    
		    		if($firstallele ne '?'){
					$count{"$firstallele"} += 1;
					$totalallele+=1;			
				}
				if($secondallele ne '?'){
					$count{"$secondallele"} +=1;
					$totalallele+=1;
				}
			}
 		}
		foreach $num (values(%count)){
			if($num < ($totalallele*$cutoff)) {
				$uselocus[$j] =0;
			}
 		}	
		$j+=1;
	}		
	$nloci = $j;
	$nuse = 0;
	for($j=0;$j<$nloci; $j+=1){
		if($uselocus[$j] ==1){
			$nuse+=1;
		}
	}

	if($nuse>0){

	print OUT $N;
	print OUT "$nuse\n";
	print OUT "P ";
	@pos = split(/\s+/,$P);
	for($j=0;$j<$nloci; $j+=1){
                if($uselocus[$j] ==1){
                       print OUT "$pos[$j+1] ";
                }
	}
	print OUT "\n";
  for($j=0;$j<$nloci; $j+=1){
                if($uselocus[$j] ==1){
                       print OUT $loci[$j];
                }
        }

	print OUT "\n";

	for($n=0; $n<$N; $n+=1){
		$fh = remove_whitespace($firsthap[$n]);
		@fha = split(//,$fh);
		for($j=0; $j<$nloci; $j+=1){
			if($uselocus[$j]==1){
				print OUT $fha[$j];
			}
		}
		print OUT "\n";
		$fh = remove_whitespace($secondhap[$n]);
                @fha = split(//,$fh);
                for($j=0; $j<$nloci; $j+=1){
                        if($uselocus[$j]==1){
                                print OUT $fha[$j];
                        }
                }				
		print OUT "\n";		
	}
	}
    }
}

	    
