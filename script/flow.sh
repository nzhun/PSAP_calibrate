
N=10000
for i in {1..22} X; do perl ../script/downsample.pl $N ~/Resources/gnomad/chr$i.gnomad.exomes.r2.0.1.sites.vcf.gz.2.bed.gz & done

for i in {1..22} X; do mv ~/Resources/gnomad/chr$i.gnomad.exomes.r2.0.1.sites.vcf.gz.2.bed.gz.downsample.${N}EUR.txt  data/chr$i/chr$i.gnomad.exomes.r2.0.1.sites..downsample.${N}EUR.bed & done

for i in {1..22} X; do bgzip  data/chr$i/chr$i.gnomad.exomes.r2.0.1.sites..downsample.${N}EUR.bed ;tabix data/chr$i/chr$i.gnomad.exomes.r2.0.1.sites..downsample.${N}EUR.bed.gz; done

 for i in {1..22} X ;do  bedtools intersect -a chr$i/gencode.hg19.anno.chr$i.base0.bed.gz  -b chr$i/chr$i.gnomad.exomes.r2.0.1.sites..downsample.${N}EUR.bed.gz -wao  -header|awk 'BEGIN{FS="\t";OFS="\t"}{if($1 ~/^#/){$NF=$NF"\tfreq_het\tfreq_hom";print;next;} if($1==$22 &&$2==$23 && $4==$25 && $5==$26){$22=$29;$23=$30}else{$22=$23=0; } print}'| cut -f 1-23 > chr$i/gencode.hg19.anno.chr$i.downsample.${N}EUR.bed & done
 


for i in {1..22} X; do 
awk 'BEGIN{FS="\t";OFS="\t"}{b=$9;if($9 ~/;/){n=split($9,a,";");for(i=1;i<n+1;i++){ print a[i],$(NF-2),$(NF-1),$NF;}}else{print b,$(NF-2),$(NF-1),$NF}}' chr$i/gencode.hg19.anno.chr$i.downsample.${N}EUR.bed |sort -k1,1d -k2,2n -k3,3n -k4,4n > chr$i/gencode.hg19.anno.chr$i.downsample.${N}EUR.gene.table.bed &
done 


for i in {1..22} X
do
	
perl ../script/gather_gene.pl chr$i/gencode.hg19.anno.chr$i.downsample.${N}EUR.gene.table.bed  chr$i/gencode.hg19.anno.chr$i.downsample.${N}EUR.table.simu.bed &

done