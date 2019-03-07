setwd("/home/local/ARCS/nz2274/PSAP_Setup")
#setwd("~/")
file="data/gencode.v27.simulation.gnomad_CADD13.table.txt";
tab<-read.table(file,header=1,stringsAsFactors = F)
genes<-unique(tab$Gene)
for(i in seq(1,length(genes),by = 50)){
  sets<-genes[i:min((i+49),length(genes))]
  write.table(tab[which(tab$Gene%in%sets),],row.names = F,file = paste("data/temp/gencode.v27.simulation.gnomad_CADD13.table.",i,".txt",sep=""),quote = F,sep="\t")
}