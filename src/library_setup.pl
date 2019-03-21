use strict;
use warnings;

## get the singleton rate of the coding region of gene, k bases where k decided by number of sites occured in cadd.
## where Func is in [exonic, splicing], which sites do not have exac annotation ==> get the singleton rate
## the input is ExAC dataset.
our $step=0.05;
our $maxC=70;
our $N_bin=int($maxC/$step)+1;
our $EXAC_SIZE=60706;

sub singleton_summary {
	my ($fin,$fgene,$fout)=@_;
	my $globle_singleton=0;
	my $gene_singleton=0;
	my $lastgene="";
	my %geneMap=();
	open GIN,"$fgene";
	while(my $l=<GIN>){
		chomp($l);
		my @ss=split(/\t+/,$l);
		$geneMap{$ss[0]}=$ss[1];
	}
	close GIN;
	open IN, "zcat $fin|";
	open OUT, ">$fout";
	my $line=<IN>;
	my %map=();
	my @sets=split(/\t+/,$line);
	for(my $i=0;$i<@sets;$i++){
		$map{$sets[$i]}=$i;
	}
	if(!exists($map{"Func.wgEncodeGencodeBasicV19"})||!exists($map{"AC_Adj"})){ 
		print ">=1 keyColumn is not found, Func.refGene or ExAC_ALL";
		exit;
	}
	while($line=<IN>){
		 chomp($line);
		 @sets=split(/\t+/,$line);
		 my $gene=$sets[$map{"Gene.wgEncodeGencodeBasicV19"}];
		 if($lastgene eq ""){$lastgene=$gene}
		 if($gene !~ /$lastgene/) {
			 print OUT $lastgene."\t".$gene_singleton/$geneMap{$lastgene}."\n";
		 	 $gene_singleton=0;
			 $lastgene=$gene;
		 }else{
			 if($sets[$map{"Func.wgEncodeGencodeBasicV19"}] !~ /^exonic|^splic/){next}
			 if(length($sets[$map{"Ref"}]) != length($sets[$map{"Alt"}])){next;}
			 if($sets[$map{"AC_Adj"}]==0){next}
			 if($sets[$map{"AC_Adj"}]==1){
				 $gene_singleton=$gene_singleton+1;
				 $globle_singleton+=1;
			 }
		}
	}
	print OUT $lastgene."\t".$gene_singleton/($EXAC_SIZE*3*$geneMap{$lastgene})."\n";
	print OUT "Overall\t".$globle_singleton/($EXAC_SIZE*3*$geneMap{"ALL"})."\n";
	close IN;
	close OUT;
	
}

sub sample_generate {
	my ($freq_arr,$cadd_arr)=@_;
	my $het_cadd=0;
	my $hom_cadd=0;
	my $chet_cadd=0;
	my $i=0;
	
	for(my $i=0;$i< @$freq_arr;$i++){
		#print "CADD $i\t"."\n";
		my $rf=rand(1);
		my $cadd=@$cadd_arr[$i];
		my $freq=@$freq_arr[$i];
	#	print $rf.";".$freq.":\n";
		if($rf <(1-$freq)*(1-$freq)){
			next;
		}elsif($rf >1-$freq*$freq){
		#	push(@gts,2);
			$hom_cadd=$hom_cadd < $cadd? $cadd :$hom_cadd;
		}else{
			#push(@gts,1);
			if($het_cadd == $cadd){
				$chet_cadd=$cadd;
			}elsif($het_cadd < $cadd){
				$chet_cadd=$het_cadd;
				$het_cadd=$cadd;
			}elsif($chet_cadd < $cadd){
				$chet_cadd=$cadd
			}
		#$i=$i+1;
		}
	}
	#print "33CADD\t".$het_cadd."\n";
	my @rs=(int($het_cadd/$step),int($chet_cadd/$step),int($hom_cadd/$step));
#	print "CADD\t".$het_cadd."\t$chet_cadd\t$hom_cadd\n";
#	print "CADD\t".join("\t",@rs)."\n";
	return (\@rs); 
}

sub cohort_generate {
	my ($freq_arr,$cadd_arr,$N_samples)=@_;
	my @fhet=(0) x $N_bin;
	my @fchet=(0) x $N_bin;
	my @fhom=(0) x $N_bin;
	
	for(my $i=0;$i<$N_samples;$i++){
		#print "DDD ".$freq_arr."\t".@$freq_arr."\t".$cadd_arr."\t".@$cadd_arr."\n";
		my ($id_het,$id_chet,$id_hom)=@{sample_generate($freq_arr,$cadd_arr)};
	#	print "$i\t$id_het\t$id_chet\t$id_hom\n";
		for(my $j=0;$j<$id_het;$j++){$fhet[$j]=$fhet[$j]+1;}
		for(my $j=0;$j<$id_chet;$j++){$fchet[$j]=$fchet[$j]+1;}
		for(my $j=0;$j<$id_hom;$j++){$fhom[$j]=$fhom[$j]+1;}
	}
	
	for(my $i=0;$i<$N_bin;$i++){
		$fhet[$i]=$fhet[$i]/$N_samples;
		$fchet[$i]=$fchet[$i]/$N_samples;
		$fhom[$i]=$fhom[$i]/$N_samples;
	}
	#print "length: ".@fhet."\t".@fchet."\t".@fhom."\n";
#	print "het pop: ".@fhet."\n";
	my @rs=(\@fhet,\@fchet,\@fhom);
	return \@rs;
}

sub simulation_sample {
	my ($fin,$fout,$fsingle,$Ntimes,$chr)=@_;
	## per line per variant, per cadd, annotation, per frequency
	#my %map=(); # gene=> singleton
	my %geneMap=();
	open SIN,$fsingle;
	while(my $l=<SIN>){
		chomp($l);
		my @ss=split(/\t/,$l);
		$geneMap{$ss[0]}=$ss[1];
	}
	close SIN;
	if($chr ==0){
		open FIN,"zcat $fin|";
	}else{
		open FIN,"tabix $fin $chr -h|";
	}
	
	open OUT,"> $fout";
	print OUT "### GeneName\t, CADD from 0 to $maxC by 0.05\n";
	#print OUT "GeneName\t".join("\t",)
	my @cadd_arr=();
	my @freq_arr=();
	
	my $line=<FIN>;
	my %map=();
	my @sets=split(/\t+/,$line);
	for(my $i=0;$i<@sets;$i++){
		$map{$sets[$i]}=$i;
	}
	
	if(!exists($map{"Func.wgEncodeGencodeBasicV19"})||!exists($map{"ExAC_ALL"})){ print ">=1 keyColumn is not found, Func.refGene or ExAC_ALL";exit;}
	my $lastgene="";
	while(my $line=<FIN>){
		## got genename, .....
		chomp($line);
		@sets=split(/\t+/,$line);
		my $gene=$sets[$map{"Gene.wgEncodeGencodeBasicV19"}];
		my $freq=$sets[$map{"ExAC_ALL"}];
		my $cadd=$sets[$map{"CADD_Phred"}];
		
		if($lastgene ne "" && $gene !~ /$lastgene/){
			my ($fhet,$fchet,$fhom)=@{cohort_generate(\@freq_arr,\@cadd_arr,$Ntimes)};
			print OUT "Dominant\t".$lastgene."\t".join("\t",@$fhet)."\n";
			print OUT "Compound-Het-Recessive\t".$lastgene."\t".join("\t",@$fchet)."\n";
			print OUT "Hom-Recessive\t".$lastgene."\t".join("\t",@$fhom)."\n";
			@freq_arr=();
			@cadd_arr=();
			$lastgene=$gene;
		}
		if($lastgene eq ""){
			$lastgene=$gene;
		}	
	 	if($sets[$map{"Func.wgEncodeGencodeBasicV19"}] !~ /^exonic|^splic/){next}
	 	if($freq eq "." || $freq==0 ){
			if(exists($geneMap{$gene})){
				$freq=$geneMap{$gene}
			}else{
				$freq=$geneMap{"Overall"}
			}
		}
		if($cadd eq "."){next}
	
		push @freq_arr, $freq;
		push @cadd_arr, $cadd;
	}
	#print "start2 ".$lastgene."\n";
	my ($fhet,$fchet,$fhom)=@{cohort_generate(\@freq_arr,\@cadd_arr,$Ntimes)};
	print OUT "Dominant\t".$lastgene."\t".join("\t",@$fhet)."\n";
	print OUT "Compound-Het-Recessive\t".$lastgene."\t".join("\t",@$fchet)."\n";
	print OUT "Hom-Recessive\t".$lastgene."\t".join("\t",@$fhom)."\n";
	close FIN;	
	close OUT;
 }




sub generate_gene {
	my ($fcadd,$fgene)=@_;
	open FIN," zcat $fcadd |";
	open FOUT," > $fgene";
	
	my $line=<FIN>;
	chomp($line);
	my %map=();
	my @sets=split(/\t+/,$line);
	for(my $i=0;$i<@sets;$i++){
		$map{$sets[$i]}=$i;
	}
	if(!exists($map{"Func.wgEncodeGencodeBasicV19"})){ print ">=1 keyColumn is not found, Func.refGene ";exit;}
	my $lastgene="";
	my $Nbase=0;
	my $totalN=0;
	while(my $line=<FIN>){
		chomp($line);
		@sets=split(/\t/,$line);
		my $gene=$sets[$map{"Gene.wgEncodeGencodeBasicV19"}];
		if($sets[$map{"Func.wgEncodeGencodeBasicV19"}] !~ /^exonic|^splic/){next;}
		if($gene !~ /$lastgene/ && $lastgene ne ""){
			print FOUT $lastgene."\t".$Nbase."\n";
			$Nbase=0;
			$lastgene=$gene;
		}elsif($lastgene eq "" ){
			$lastgene=$gene;
		}
			$Nbase=$Nbase+1;
			$totalN=$totalN+1;
	}
	print FOUT $lastgene."\t".$Nbase."\n";
	print FOUT "ALL\t".$totalN."\n";
	close FIN;
	close FOUT;
}


## read cadd
my ($fcadd ,$fexac)=@ARGV;
print $fcadd."\t".$fexac."\n";
my $f_singleton="test_singletonRate.txt";
#my $fcadd="";
my $popscore_out="test_popscore.txt";
my $N=1000000;

if(! -e $fcadd){
	print $fcadd." 1 cannot find!\n";
	exit;
}
if(! -e $f_singleton){
	#my $fexac="";
	if(! -e $fexac){
		print $fexac." 2 cannot find!\n";
		exit;
	}
	my $fgene="test_genesize.txt";
	if( ! -e $fgene){
		generate_gene($fcadd,$fgene);
		print "gene profile done\n";
	}
	singleton_summary($fexac,$fgene,$f_singleton);
		print "rate profile done\n";
}


my $chr=0;
simulation_sample($fcadd,$popscore_out,$f_singleton,$N,$chr);
	print "popscore profile done\n";