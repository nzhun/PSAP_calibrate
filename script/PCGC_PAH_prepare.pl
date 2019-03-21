
## spike in the HGMD or clinvar variants
## the probablity of picking the variant out of an exome
##get the largest CADD score from all genes, all individuals.

## randomly chose the variant from clinvar, got its cadd
#!/usr/bin/perl;
use strict;
use warnings;
open IN,"/home/local/ARCS/nz2274/PSAP_Setup/vcf/PAH_control_VQSR.vcf.gz.2.txt";
open OUT,">/home/local/ARCS/nz2274/PSAP_Setup/vcf/PAH_control_VQSR.bed";
my $line=<IN>;
my %header=();
chomp($line);
my @sets=split(/\s+/,$line);
print OUT join("\t",@sets[(0..5,13..18,24..29,93,94)])."\n";
my $i=0;
foreach my $item(@sets){
  $header{$item}=$i;
  $i=$i+1;
}


while($line=<IN>){
    chomp($line);
    my @sets=split(/\s+/,$line);
    my $varfunc=$sets[$header{"VarFunc"}];

    if($varfunc ne "exonic" && $varfunc ne "splicing"){next}

    my $filter=$sets[$header{"FILTER"}];
    if($filter eq "VQSRTrancheSNP99.90to100.00" ||
        $filter eq "VQSRTrancheSNP99.80to99.90" ||
        $filter eq "VQSRTrancheINDEL99.90to100.00" ||
        $filter eq "VQSRTrancheINDEL99.80to99.90" ||
        $filter eq "VQSRTrancheINDEL99.70to99.80"  ) {
        next;
    }

    my $segdup=$sets[$header{"SegDup"}];
    if($segdup ne "."){
        my @scores=split(/:/,$segdup);
    #    print $segdup."\n";
        if(  $scores[1] =~ /^[0-9,.E]+$/ && $scores[1]>0.95){next;}
   }
    #    print $varfunc."\n";
   print OUT join("\t",@sets[(0..5,13..18,24..29,93,94)])."\n";
}
close OUT;
close IN;
