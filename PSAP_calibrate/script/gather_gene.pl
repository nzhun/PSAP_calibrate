#/usr/bin/perl
## for gnomAD
use strict;
use warnings;
my $totalN=123136-40860;
my $N_sites=104210332;
my $global_het=19218/($totalN*$N_sites);
my $global_hom=4644/($totalN*$N_sites);


my $het_file="/home/local/ARCS/nz2274/PSAP_Setup/data/gnomAD_het_sington.gene.rate.csv";
my $hom_file="/home/local/ARCS/nz2274/PSAP_Setup/data/gnomAD_homo_sington.gene.rate.csv";
my %het_map=();

open IN_Het, "$het_file" or die "cannot open $het_file";
my $l=<IN_Het>;
while($l=<IN_Het>){
  chomp($l);
  my @sets=split(/,/,$l);
#  print $sets[0]."\t".$sets[2]."\n";
  $het_map{$sets[0]}=$sets[2];
}
close IN_Het;
print "load 1\n";

my %hom_map=();
open IN_HOM, "$hom_file" or die "cannot open $het_file";;
$l=<IN_HOM>;
while($l=<IN_HOM>){
  chomp($l);
  my @sets=split(/,/,$l);
  #  print $sets[0]."\t".$sets[2]."\n";
  $hom_map{$sets[0]}=$sets[2];
}
close IN_HOM;
print "load 2\n";


my $file=$ARGV[0]; #"/home/local/ARCS/nz2274/PSAP_Setup/test/gencode.v27.simulation.table.temp.bed.gz"
open IN, "less $file| " or die "can not open $file!";
open OUT2, ">".$ARGV[1]; # "/home/local/ARCS/nz2274/PSAP_Setup/test/gencode.v27.simulation.table.out.txt";
open OUT3, ">".$ARGV[1].".allele.txt";
my %map=();
my %map2=();
my $lastgene="";
while(my $line=<IN>){
  chomp($line);
  my @sets=split(/\s+/,$line);
  my $gene=$sets[0];
  if($gene eq "Gene.refGene"){next;}
  my $score=$sets[1];
  my $freq=$sets[2];
  my $freq_het2=$sets[2];
  my $freq_hom2=$freq_het2;
  my $freq_hom=$sets[3];
  if($freq eq "." || $freq eq "0" || $freq eq "NA"){
    $freq=$global_het;
    if(exists($het_map{$gene})){$freq=$het_map{$gene}}
    $freq_het2=$freq;
    $freq_hom2=$global_hom;
    if(exists($hom_map{$gene})){$freq_hom2=$hom_map{$gene}}
  }else{
     $freq_het2=$sets[2]+$sets[3];
     $freq_hom2=$freq_het2;

  }
  if($freq_hom eq "." || $freq_hom eq "0" || $freq_hom eq "NA"){
    $freq_hom=$global_hom;
    if(exists($hom_map{$gene})){$freq_hom=$hom_map{$gene}}
  }
  my $key=$gene."\t".$freq."\t".$freq_hom."\t"."$score";
  my $key2=$gene."\t".$freq_het2."\t".$freq_hom2."\t"."$score";

  if($lastgene eq ""){$lastgene=$gene;}
  if($lastgene ne $gene){
    foreach my $k(keys %map){
       print OUT2 $k."\t".$map{$k}."\n";
    }
    %map=();

    foreach my $k(keys %map2){
       print OUT3 $k."\t".$map2{$k}."\n";
    }
    %map2=();

    $lastgene=$gene;
  }
  if(exists($map{$key})){
        $map{$key}=$map{$key}+1;
  }else{
      $map{$key}=1;
  }


  if(exists($map2{$key2})){
        $map2{$key2}=$map2{$key2}+1;
  }else{
      $map2{$key2}=1;
  }
}
close IN;
close OUT2;
close OUT3;
