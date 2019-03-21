use strict;
use warnings;

## get the singleton rate of the coding region of gene, k bases where k decided by number of sites occured in cadd.
## where Func is in [exonic, splicing], which sites do not have exac annotation ==> get the singleton rate
## the input is ExAC dataset.
our $step=0.05;
our $maxC=70;
our $N_bin=int($maxC/$step)+1;
our $EXAC_SIZE=60706;
our $Total_Base=0;
our $Total_singleton=0;
our $COHORT_SIZE=60706;#125748;

### for exac ###
our %dataset_setting=();
$dataset_setting{"dataset_ac"}="AC_Adj";
$dataset_setting{"genename"}="Gene.wgEncodeGencodeBasicV19";
$dataset_setting{"func"}="Func.wgEncodeGencodeBasicV19";
$dataset_setting{"exonicfunc"}="ExonicFunc.wgEncodeGencodeBasicV19";
$dataset_setting{"ref"}="Ref";
$dataset_setting{"alt"}="Alt";
$dataset_setting{"pos"}="Start";
$dataset_setting{"avg_freq"}="ExAC_ALL";

### gnomAD settting  
#our %dataset_setting=();
#$dataset_setting{"dataset_ac"}="AC_Adj";
#$dataset_setting{"genename"}="Gene.wgEncodeGencodeBasicV19";
#$dataset_setting{"func"}="Func.wgEncodeGencodeBasicV19";
#$dataset_setting{"exonicfunc"}="ExonicFunc.wgEncodeGencodeBasicV19";
#$dataset_setting{"ref"}="Ref";
#$dataset_setting{"alt"}="Alt";
#$dataset_setting{"avg_freq"}="gnomAD_exome_ALL";
#$dataset_setting{"pos"}="From";




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
	if(!exists($map{$dataset_setting{"func"}})||!exists($map{$dataset_setting{"dataset_ac"}})){ 
		print ">=1 keyColumn is not found, ".$dataset_setting{"func"}. " or ".$dataset_setting{"dataset_ac"}."\n";
		exit;
	}
	while($line=<IN>){
		 chomp($line);
		 @sets=split(/\t+/,$line);
		 my $gene=$sets[$map{$dataset_setting{"genename"}}];
		 if($lastgene eq ""){$lastgene=$gene}
		 if($gene !~ /$lastgene/) {
			 print OUT $lastgene."\t".$gene_singleton/$geneMap{$lastgene}."\n";
		 	 $gene_singleton=0;
			 $lastgene=$gene;
		 }else{
			 if($sets[$map{$dataset_setting{"func"}}] !~ /^exonic|^splic/){next}
			 if(length($sets[$map{$dataset_setting{"ref"}}]) != length($sets[$map{$dataset_setting{"alt"}}])){next;}
			 if($sets[$map{$dataset_setting{"dataset_ac"}}]==0){next}
			 if($sets[$map{$dataset_setting{"dataset_ac"}}]==1){
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



sub freq_0 {
	my ($freqs)=@_;
	## it is an array of the probablity  of variants to be heterozygous
	my $f0=1;
	foreach my $f(@$freqs){
		$f0=$f0*(1-$f);
	}
	return ($f0);
}

sub freq_all {
	my $freqs=$_[0];
#	print join(":",@{$freqs})."\n";
	## it is an array of the probablity  of variants to be heterozygous
	my $f0=1;
	foreach my $f(@{$freqs}){
		$f0=$f0*($f);
	}

	return $f0;
}


sub freq_1 {
	my ($freqs)=@_;
	## it is an array of the probablity  of variants to be heterozygous
	my $f1=0;
	for (my $i=0;$i<@$freqs;$i++){
		my $sfreq=1;
		for(my $j=0;$j<@$freqs;$j++){
			if($i==$j){
				$sfreq=$sfreq*@$freqs[$j];
			}else{
				$sfreq=$sfreq*(1-@$freqs[$j]);
			}
		}
		$f1=$f1+$sfreq;
	}
	return ($f1);
}


sub join_freq {
	my ($f0_ad,$f1_ad)=@_;
	my @f0=@$f0_ad;
	my @f1=@$f1_ad;
	my $f=0;
	for(my $i=0;$i<@f0;$i++){
		my $sf=1;
		for(my $j=0;$j<@f1;$j++){
			if($i==$j){
				$sf=$sf*$f1[$j];
			}else{
				$sf=$sf*$f0[$j];
			}	
		}
		$f=$f+$sf;
	}
	return ($f);
}



sub print_string {
	my %table=%{$_[0]};
	my $outstr="";
	my @keys=sort { $b <=> $a } keys %table;
	my $lastk=0;
	my  $sum=0;
	for(my $i=$N_bin-1;$i>-1;$i--){
		if(!exists($table{$i})){
			$table{$i}=$lastk;
			$outstr=$sum."\t".$outstr;
		}else{
			$lastk=$table{$i};
			$sum=$sum+$table{$i};
			$outstr=$sum."\t".$outstr;
		}
	}
 #@keys=sort { $b <=> $a } keys %table;
# my $lastk=1;

 #my $outstr="";
 
 return ($outstr);
}



sub popscore_enume {
	my ($cadd_het_ad,$cadd_hom_ad)=@_;
	my %cadd_het=%{$cadd_het_ad};
	my %cadd_hom=%{$cadd_hom_ad};
	my @keys = sort {$b <=> $a} keys %cadd_het;
	my @f0=();
	my @f1=();
	
	for (my $k=0;$k<@keys;$k++){
		my @fc=@{$cadd_het{$keys[$k]}};
		push(@f0,$fc[0]);
		push(@f1,$fc[1]);
	}
	my %pop_het=();
	my %pop_chet=();
	$pop_het{$keys[0]}=1-$f0[0];
	$pop_chet{$keys[0]}=1-$f0[0]-$f1[0];
	#print join(":",@f0)."\n";

	if(@keys >1){
		for(my $k=1;$k<@keys;$k++){
			my @freq_0=@f0[0..($k-1)];
			my @freq_1=@f1[0..($k-1)];
			my $chet_f_lg1=join_freq(\@freq_0,\@freq_1)*(1-$f0[$k]); ## only once true
			my $chet_f_lg0=freq_all(\@freq_0)*(1-$f0[$k]-$f1[$k]);
	
			my $chet=$chet_f_lg1+$chet_f_lg0;
			my $het=freq_all(\@freq_0)*(1-$f0[$k]);
			$pop_het{$keys[$k]}=$het;
			$pop_chet{$keys[$k]}=$chet;
			#$pop_het{$keys[$k]}=$het+$pop_het{$keys[$k-1]};
			#$pop_chet{$keys[$k]}=$chet$pop_chet{$keys[$k-1]};
		}
	}
 #  exit;
	
	my @f0_hom=();
	my @f1_hom=();	
	for (my $k=0;$k<@keys;$k++){
		my $fc=$cadd_hom{$keys[$k]};
		push(@f0_hom,$fc);
	#	push(@f1_hom,$fc[1]);
	}
	my %pop_hom=();
	$pop_hom{$keys[0]}=1-$f0_hom[0];
	for(my $k=1;$k<@keys;$k++){
		my @freq_0=@f0_hom[0..($k-1)];
		#my @freq_1=@f1[0..($k-1)];
		my $hom=freq_all(\@freq_0)*(1-$f0_hom[$k]);
		$pop_hom{$keys[$k]}=$hom;
		#$pop_hom{$keys[$k]}=$hom+$pop_hom{$keys[$k-1]};
	}
	my @rs=(\%pop_het,\%pop_chet,\%pop_hom);
	#print "popscore:CADD_index\tCADD\thet\tCHET\tHOMO\n";
	#my @kk=sort {$a <=>$b }(keys %pop_het) ;
	#foreach my $k(@kk){
#		print $k."\t".($k*0.05)."\t".$pop_het{$k}."\t".$pop_chet{$k}."\t".$pop_hom{$k}."\n";#
#	}
 
 
	return \@rs;
}

sub cadd_enume {
	 my ($fhet,$fhom)=@_;
	 my %het=%{$fhet};
	 my %hom=%{$fhom};

	 ## possion test 1 times for dominant
	 ## possion test 2 times for dominant and 0 in the larger CADD, or 1 times dominant in the current cadd and 1 time dominant in larger cadd 
	 my %cadd_het=();
	 foreach my $key(keys %het){
		 my @cfreq=@{$het{$key}};
		 ## none larger cadd is selcted, at least one current cadd is select
		 my $f0=freq_0(\@cfreq);
		 ## only one is selected
		 my $f1=freq_1(\@cfreq);
		 my @rs=($f0,$f1);
		 $cadd_het{$key}=\@rs;
	 }
	 my %cadd_hom=();
	 foreach my $key(keys %hom){
			my $f0=freq_0($hom{$key});
			$cadd_hom{$key}=$f0;
	 }
	 
	 my @rs=(\%cadd_het,\%cadd_hom);
	# print "CADD\tHet_none\tHet_1\n";
	# my @kk = sort {$a<=> $b} keys %cadd_het;
	# foreach my $k(@kk){
	#	 print $k."\t".($k*0.05)."\t".join("\t",@{$cadd_het{$k}})."\n";
	# }
	 
	 return \@rs;
}

sub load_header {
	my %map=();
	my $file=$_[0];
	open IN_TEMP, "tabix  $file   -H |";
	my $line=<IN_TEMP>;
	chomp($line);
	my @sets=split(/\t/,$line);
	for(my $i=0;$i<@sets;$i++){
		$map{$sets[$i]}=$i;
	#	print $sets[$i]."\t".$i."\n";
	}
	close IN_TEMP;
	return (\%map);
}


sub load_freq {
	my ($ffreq,$region,$ad_freq_header)=@_;
	my %freq_header=%{$ad_freq_header};
	my %var_freq=();
	my $singleton=0;
	my $mis=0;
	my $freq_lof=0;
	my $freq_het_lof=0;
	my $freq_hom_lof=0;
	open IN_FREQ, "tabix  $ffreq $region |";
	while(my $fline=<IN_FREQ>){
		chomp($fline);
		my @fsets=split(/\t/,$fline);
		my $v_freq=$fsets[$freq_header{$dataset_setting{"avg_freq"}}];
	#	print $dataset_setting{"avg_freq"}."\t".$freq_header{$dataset_setting{"avg_freq"}}."\t".$v_freq."\n"; 

		if($fsets[$freq_header{$dataset_setting{"dataset_ac"}}]==1){$singleton=$singleton+1;}
	#	my @freqs=($fsets[$freq_header{"gnomad_exome_all"}],$fsets[$freq_header{"gnomad_exome_nfe"}],$fsets[$freq_header{"gnomad_exome_het"}],$fsets[$freq_header{"gnomad_exome_hom"}])
		$var_freq{$fsets[$freq_header{$dataset_setting{"pos"}}].$fsets[$freq_header{$dataset_setting{"alt"}}]}=$v_freq;
		if($v_freq eq "."){
			$mis=$mis+1;
		}else{
			if(length($fsets[$freq_header{$dataset_setting{"ref"}}]) != length($fsets[$freq_header{$dataset_setting{"alt"}}])){
				if($fsets[$freq_header{$dataset_setting{"func"}}] =~ /^exonic|^splic/ && $fsets[$freq_header{$dataset_setting{"exonicfunc"}}] =~ /^frame|^stop|^splic/ ){
					$freq_het_lof=$freq_het_lof+2*$v_freq*(1-$v_freq);
					$freq_hom_lof=$freq_hom_lof+$v_freq*($v_freq);
				}
			}
		}
		## add a indel frame shift key. then need annotation for each variant
		## get the singleton number
		## get the map of variants => freq
	}
	close IN_FREQ;
	my @rs=();
	if($mis>0 && $singleton==0){
		@rs=(-1,-1);
	}else{
		@rs=(\%var_freq,$singleton);
	}
	return (\@rs);
}

sub simulation_sample {
	my ($fin,$ffreq,$region,$map_lof)=@_;
	## per line per variant, per cadd, annotation, per frequency
	#my %map=(); # gene=> singleton
	## read lof, gene=> cadd
	## read freq, variants=> freq
	## read fin, variants=> cadd  => get total number of bases in a gene
	## region, a genomic region, like a gene or a transcript.
	## read fsingle, gene=> singleton rate, summary from freq file
	
	### get the header of ffreq 
#	
	print "tabix  $ffreq  -H\n"; 
	my %freq_header=%{load_header($ffreq)};
#   print join("\t", keys %freq_header)."\n";
	my %cadd_header=%{load_header($fin)};
	##per varaint
	## get the cadd score
	## get their  frequency


	my ($ad_var_freq,$singleton)=@{load_freq($ffreq,$region,\%freq_header)};
	if($singleton < 0) {next;}
	my %var_freq=%{$ad_var_freq};
	my $mis=0;
	my %cadd_het_freq=(); ## cadd to an arry of frequency
	my %cadd_hom_freq=();## cadd to an arry of frequency
	my %extra_cadd_freq=();
	my %cadd_chet_reccurent=();
	my $nbase=0;
	print "tabix $fin $region \n";
	open IN_CADD, "tabix $fin $region |";
	while(my $line=<IN_CADD>){
		chomp($line);
		my @csets=split(/\t/,$line);
		my $qkey=$csets[$cadd_header{"Start"}].$csets[$cadd_header{"Alt"}];
		my $rawkey=$csets[$cadd_header{"CADD_Phred"}];
		if($rawkey eq "."){next}
		my $key=int($rawkey*20)+1;
	#	print "1:\t".$key."\t".$rawkey."\t$qkey\t"."\n";
		if(!exists($var_freq{$qkey})|| $var_freq{$qkey}==0 ){
			if(exists($extra_cadd_freq{$key})){
				$extra_cadd_freq{$key}=$extra_cadd_freq{$key}+1;
			}else{
				$extra_cadd_freq{$key}=1;
			}
		#	print "0:\t".$key."\t".$rawkey."\t$qkey\t"."\n";
		}else{
			my $freq_alt=$var_freq{$qkey};
			my @het_freqs=();
			my @hom_freqs=();
			if(exists($cadd_het_freq{$key})){
				@het_freqs=@{$cadd_het_freq{$key}};
				@hom_freqs=@{$cadd_hom_freq{$key}};
			}
			push(@het_freqs,$freq_alt*(1-$freq_alt)*2);
			push(@hom_freqs,$freq_alt*$freq_alt);
		   $cadd_het_freq{$key}=\@het_freqs;
			$cadd_hom_freq{$key}=\@hom_freqs;
		#	print "2:\t".$key."\t".$rawkey."\t$qkey\t($freq_alt*$freq_alt)"."\n";
		}
		
		$nbase=$nbase+1;
	}
	
	close IN_CADD;
	if($nbase==0){next;}
	print "GeneSummary: ".$singleton."\t".$nbase."\t".$COHORT_SIZE."\n";	
	my $singleton_rate=$singleton/($nbase*$COHORT_SIZE);
	foreach my $k(keys %extra_cadd_freq){
		my @vfreq_het=();
		my @vfreq_hom=();
		if(exists($cadd_het_freq{$k})){
			@vfreq_het=@{$cadd_het_freq{$k}};
			@vfreq_hom=@{$cadd_hom_freq{$k}};
		}
		push(@vfreq_het,($singleton_rate*(1-$singleton_rate)*2) x $extra_cadd_freq{$k} );
		push(@vfreq_hom,($singleton_rate*($singleton_rate)) x $extra_cadd_freq{$k} );
		$cadd_het_freq{$k}=\@vfreq_het;
		$cadd_hom_freq{$k}=\@vfreq_hom;
	}
	
	#print "CADD\tvaraintFreq\n";
	#my @kk = sort { $a <=> $b} keys %cadd_hom_freq;
	#foreach my $k(@kk){
	#	print $k."\t".$k*0.05."\t".join("\t",@{$cadd_hom_freq{$k}})."\n";
	#}
	
	my @rs=@{cadd_enume(\%cadd_het_freq, \%cadd_hom_freq)};
	return \@rs;	
 }
 ### /home/local/ARCS/nz2274/PSAP_Setup/Data/gnomAD_exome.het.singleton.anno.txt.gz
 
 ### /home/local/ARCS/nz2274/PSAP_Setup/Data/gnomAD_exome_singleton.bed.gz
 
 
 ##  Data/gnomAD_exome.homo.singleton.anno.bed
 ## /home/local/ARCS/nz2274/PSAP_Setup/Data/gencode.v27.sim2.anno.bed.gz
 
 ## read cadd
	my ($fcadd ,$fdataset,$ftrans)=@ARGV;
	print $fcadd."\t".$fdataset."\t".$ftrans."\n";
	my $f_singleton="test_enume_singletonRate.txt";
	#my $fcadd="";
	my $popscore_out="test_enume_popscore.txt";
	# my $fout="test_"
	my $N=1000000;

	if(! -e $fcadd){
		print $fcadd." 1 cannot find!\n";
		exit;
	}
	if(! -e $f_singleton){
         #my $fexac="";
			if(! -e $fdataset){
				print $fdataset." 2 cannot find!\n";
				exit;
			}
			my $fgene="test_enume_genesize.txt";
			if( ! -e $fgene){
				generate_gene($fcadd,$fgene);
				print "gene profile done\n";
			}
			singleton_summary($fdataset,$fgene,$f_singleton);
			print "rate profile done\n";
	}


# my $chr=0;
# simulation_sample($fcadd,$popscore_out,$f_singleton,$N,$chr);
#         print "popscore profile done\n";
			
			
	#my $fsingle="";
	#my $flof="";
	my %map_lof=();
	## load f lof to a map
	 #my $ftrans="";
	open INTR, "$ftrans";
	open OUT, ">$popscore_out";
	my $oh="Transcript\tMODEL\t";
	for(my $i=0;$i<$N_bin+1;$i++){$oh=$oh."\t<".($i*0.05);}
	print OUT "$oh\n";
	while(my $line=<INTR>){
		chomp($line);
		my @sets=split(/\t/,$line); ## get the transcript name and the region
		print "start\n";
		my ($cadd_het,$cadd_hom)=@{simulation_sample($fcadd,$fdataset,join(" ",@sets[1..(@sets-1)]),\%map_lof)};
		my ($fhet,$fchet,$fhom)=@{popscore_enume($cadd_het,$cadd_hom)};
		#print OUT "$sets[0]\t";
		my $outstr=print_string($fhet);
		print OUT $sets[0]."\tDominant\t1\t".$outstr."\n";
		$outstr=print_string($fchet);
		print OUT $sets[0]."\tREC_CHET\t1\t".$outstr."\n";
		$outstr=print_string($fhom);
		print OUT $sets[0]."\tREC_HOM\t1\t".$outstr."\n";
#		exit;
	}
	close OUT;
	close INTR;
	
	 #$fin,$ffreq,$region,$map_lof