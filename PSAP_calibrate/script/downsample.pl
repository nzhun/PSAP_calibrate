#! /usr/bin/perl;
use strict;
use warnings;
my $file=$ARGV[1];  #"/home/local/users/jw/resources/gnomad/Exome/gnomad.exomes.r2.0.1.sites.vcf.gz.2.bed.gz";
open IN,"zcat $file|" or die "cannot open $file!\n";
open OUT,">$file.downsample.$ARGV[0]EUR.txt";
my $line=<IN>;
my %map=();
my @arr=split(/\s+/,$line);
for(my $i=0;$i<@arr;$i++){
 $map{$arr[$i]}=$i;
}
print OUT join("\t",@arr[0..6])."\tfreq_het"."\tfreq_hom\n";

my $AEUR=55860*2;
my $EUR=$ARGV[0]*2;
my $AFR=7652*2;
my $AMR=16791*2;
my $ASJ=4925*2;
my $EAS=8624*2;
my $FIN=11150*2;
my $SAS=15391*2;
my $OTH=2743*2;
my @pops=("NFE","AFR","AMR","ASJ","EAS","FIN","OTH","SAS");
while($line=<IN>){
  my @sets=split(/\s+/,$line);
  #my $eur_het=$sets[$map{"AC_NFE"}]-2*$sets[$map{"Hom_NFE"}]-$sets[$map{"Hemi_NFE"}];
  #my $afr_het=$sets[$map{"AC_AFR"}]-2*$sets[$map{"Hom_AFR"}]-$sets[$map{"Hemi_AFR"}];
  #my $amr_het=$sets[$map{"AC_AMR"}]-2*$sets[$map{"Hom_AMR"}]-$sets[$map{"Hemi_AMR"}];
  #my $asj_het=$sets[$map{"AC_ASJ"}]-2*$sets[$map{"Hom_ASJ"}]-$sets[$map{"Hemi_ASJ"}];
  #my $eas_het=$sets[$map{"AC_EAS"}]-2*$sets[$map{"Hom_EAS"}]-$sets[$map{"Hemi_EAS"}];
  #my $FIN_het=$sets[$map{"AC_FIN"}]-2*$sets[$map{"Hom_FIN"}]-$sets[$map{"Hemi_FIN"}];
  #my $OTH_het=$sets[$map{"AC_FIN"}]-2*$sets[$map{"Hom_OTH"}]-$sets[$map{"Hemi_OTH"}];
  #my $FIN_het=$sets[$map{"AC_SAS"}]-2*$sets[$map{"Hom_SAS"}]-$sets[$map{"Hemi_SAS"}];
  my $tc_het=0;
  my $tc_hom=0;
  my $t_N=0;
  foreach my $pop(@pops){
    my $hemi=$sets[$map{"Hemi_".$pop}];
    my $het=$sets[$map{"AC_".$pop}];
    my $hom=$sets[$map{"Hom_".$pop}];
    my $AN=0;
    if($het eq "."){$het=0;}
    if($hom eq "."){$hom=0;}
    if($hemi eq "."){$hemi=0;}
    $het=$sets[$map{"AC_".$pop}]-2*$sets[$map{"Hom_".$pop}]-$hemi;
    $hom=$sets[$map{"Hom_".$pop}];
    $AN=$sets[$map{"AN_".$pop}];
  #  print $AN."\t".$het."\t".$hom."\t".$sets[1]."\t".$map{"AN_".$pop}."\t".$sets[$map{"AN_".$pop}]."\t".$sets[44]."\t".join("\t",@sets[0..5])."\n";
    if($het<0){$het=0;}
    if($pop eq "NFE"){
      if($AN>$EUR){
        $t_N=$t_N+$EUR;
        $tc_het=$tc_het+$het*($EUR/$AN);
        $tc_hom=$tc_hom+$hom*($EUR/$AN);
      }else{
        $t_N=$t_N+$AN;
        $tc_het=$tc_het+$het;
        $tc_hom=$tc_hom+$hom;
      }
    }else{

      $t_N=$t_N+$AN;
      $tc_het=$tc_het+$het;
      $tc_hom=$tc_hom+$hom;
    }

  }
  if($t_N==0){next;}
  my $fhet=$tc_het/$t_N;
  my $fhom=$tc_hom*2/$t_N;
  print OUT join("\t",@sets[0..6])."\t".$fhet."\t".$fhom."\n";
}

close OUT;
close IN;
