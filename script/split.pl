#!/usr/bin/perl;
use strict;
use warnings;
my %hap=();
open IN, "zcat KnownCannoical.gz|";
while(my $lin=<IN>){
   my @sets= split(/\s+/,$lin);

   $hap{$sets[4]}=substr($sets[0],3);
}
close IN;
my %handlers=();
my @chrs=(1..22,'X','U');
foreach my $i(@chrs){
  my $file="/home/local/ARCS/nz2274/PSAP_Setup/data/gencode.v27.simulation.gnomad_CADD13.chr$i.table.txt";
  open my $OUT, ">$file";
  $handlers{$i}=$OUT;
}

my $lastg="";
my $OUT=0;
open IN,"/home/local/ARCS/nz2274/PSAP_Setup/data/gencode.v27.simulation.gnomad_CADD13.table.txt";
my $lin=<IN>;
while($lin=<IN>){
   my @sets=split(/\s+/,$lin);

   if($sets[0] ne $lastg){
     my $key='U';
     if(exists($hap{$sets[0]})){
        $key=$hap{$sets[0]};
     }else{
       $key='U'
     }
     if(!exists($handlers{$key})){
       $key='U';
     }
     $OUT=$handlers{$key};

     $lastg=$sets[0];
  #      print $sets[0]."\t$key\t"."\n";
   }

     print $OUT $lin;


}
close IN;


foreach my $key(keys %handlers){
  close $handlers{$key};
}
