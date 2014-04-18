#!/usr/bin/perl

# converts michael's cdat files into PHASE format

# perl cdat2phase.pl file.cdat file.inp

open(IN,"$ARGV[0]");
open(OUT,">$ARGV[1]");

$firstline=<IN>;
@params = split(/\s+/,$firstline);

$N=$params[0]/2;
$R= $params[1];
$pos = <IN>;

print OUT "$N\n";
print OUT "$R\n";
print OUT "P $pos";
print OUT 'S' x $R."\n";

while(<IN>){
    
    $hap1=$_; $hap1 =~ s/\s+//g;
    $hap2=<IN>; $hap2 =~s/\s+//g;
    if($hap1 =~ s/N//){
	print OUT "0\n";
    } else {
	$hap1 =~ s/D//;
	print OUT "1\n";
    }
    $hap2 =~ s/[ND]//;
    print OUT "$hap1\n";
    print OUT "$hap2\n";


}

