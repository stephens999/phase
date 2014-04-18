#!/usr/bin/perl

# usage is 
# perl summarise.pl file.res nind

open(SUMMARY,"$ARGV[0]");
$nind = $ARGV[1];

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


$ourerror = 0;
$nqxlerror = 0;
$switcherror = 0; # holds numerator of switch errors
$switchtot = 0; # holds denominator of switch errors
$ndatasets = 0;
$totalwrong = 0;
$ambigind = 0;
$totalwrongnonmissing  = 0;
$ambigindnonmissing = 0;

while(<SUMMARY>){
    $line=$_;
    $ndatasets += 1;
    @linearray=split(/\s+/,$line);
    if($nind == 0){
	$nind = $linearray[11];
	}
    if($linearray[0]>0){
	$ourerror += $linearray[1]/$linearray[0];
    }
    if($linearray[2]>0){
	$poserror += $linearray[3]/$linearray[2];
    }
    $nqxlerror += $linearray[1]/$nind;
    $switcherror += $linearray[6];
    $switchtot += $linearray[7];

    $ambigind += $linearray[0];
    $totalwrong += $linearray[1];
    $ambigindnonmissing += $linearray[9];
    $totalwrongnonmissing += $linearray[10]; 		
}

$ourerror = $ourerror/$ndatasets;
$poserror = $poserror/$ndatasets;
$nqxlerror = $nqxlerror / $ndatasets;
$switcherror = $switcherror / $switchtot; 
$shinerror = $totalwrong / $ambigind;
$shinerrornonmissing = $totalwrongnonmissing / $ambigindnonmissing;

print "File: $ARGV[0]\n";
print "Number of datasets = $ndatasets\n";

printf "%10s %10s %10s %10s %10s %10s\n",'SSD','NQXL','Pos','Switch Err','Shin SE','Shin,nomiss';
printf "%10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n", $ourerror,$nqxlerror,$poserror,1-$switcherror,$shinerror,$shinerrornonmissing;



	    
