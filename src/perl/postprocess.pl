open(INPUT,"testpost.out");

%list=('A',0,'G',0);

print %list,"\n";

$list{'C'}++;
print %list,"\n";


@Pop=();

$Nind=10;

for($i=0;$i<$Nind;$i++){    
    do{
	$Line1=<INPUT>;
	@Hap1=split(/\s+/,$Line1);
    }until($#Hap1+1 > 0);
    $Line2=<INPUT>;	
    @Hap2=split(/\s+/,$Line2);
    @Ind=();
    push(@Ind,[@Hap1]);
    push(@Ind,[@Hap2]);   
    push(@Pop,[@Ind]);	
}

$Nind=$#Pop+1;
$Nloci=$#Hap1+1;
print $Nind;
print $Nloci;
print $i;

for($ind=0; $ind<$Nind; $ind++){
    for($c=0;$c<2;$c++){
	print("$ind,$c");
	for($r=0;$r<$Nloci;$r++){
	    print($Pop[$ind][$c][$r]);
	}
	print"\n";
    }
    print"\n";
}

#%list=(0,'first item',1,'second item');
#$temp=$Hap1[1];
#print $temp;
#print $list{$temp};
#print(@Hap1,@Hap2);


