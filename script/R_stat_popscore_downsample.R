setwd("~/PSAP_Setup/") 
#setwd ("/home/local/ARCS/nz2274/PSAP_Setup/")
popscore<-function(dat){
  popscore<-rep(0,1402)
  for(i in 1:length(dat)){
    index<-as.integer(dat[i]/0.05)+1
    popscore[1:index]<-popscore[1:index]+1;
  }
  popscore<-(popscore)/(length(dat))
  return(popscore)
}
## all:  finputs="PSAP_Setup/test/gencode.v27.simulation.table.out.txt" #" c("PSAP_Setup/test/2_CADD_ExAC_.txt","PSAP_Setup/test/2_CADD13_gnomAD_.txt")
#finput="test/2_CADD_ExAC_.txt"
args = commandArgs(trailingOnly=TRUE)
finput= args[1]  #"data/temp/gencode.v27.simulation.table.out.txt" ## arg input
#for(finput in finputs){
#finput="PSAP_Setup/test/2_CADD13_gnomAD_.txt" #"PSAP_Setup/Data/source/CADD13_gnomAD_0.txt"
prefixs<-unlist(strsplit(finput,"/",fixed=T))
prefix<-prefixs[length(prefixs)]
dat<-read.table(finput,header = F,stringsAsFactors = F)
names(dat)<-c("Gene","het_freq","homo_freq","CADD","Ntime")
#  names(dat)<-c("Gene"  ,"het_freq"        ,"homo_freq"       ,"CADD"    ,"Ntime")
simu<-1402 #1000000
#N<-1402
all_het<-c()
all_hom<-c()
all_chet<-c()
genes<- unique(dat$Gene)
#genes<-"SOX17"
for(g in genes){
  # vecs<-rep(0,simu)
  if(g=="Gene" || g=="Gene.refGene"){next;}
  
  het_popscore<-rep(0,simu)
  homo_popscore<-rep(0,simu)
  chet_popscore<-rep(0,simu)
  
  print(g)
  ghet<-rep(0,simu);
  ghom<-rep(0,simu)
  gchet<-rep(0,simu)
  gdat<-dat[which(dat$Gene==g & !is.na(as.numeric(dat$CADD))),]
  if(dim(gdat)[1]<1){next;}
  gdat$CADD<-as.integer(as.numeric(gdat$CADD)/0.05)*0.05
  gdat<-gdat[order(gdat$CADD),]
  print(Sys.time())
  for(score in unique(gdat$CADD) ){
    test<-gdat[which(gdat$CADD==score),]  
    phet<-0;
    pchet<-0;
    phom<-0;
    ### get the probability of a person carried het/chet/hom with cadd=c 
    for(i in 1:dim(test)[1]){
      f_het<-as.numeric(test$het_freq[i])
      f_hom<-as.numeric(test$homo_freq[i])
      n=as.numeric(test$Ntime[i])
      phet<-phet+(1-ppois(0,lambda = f_het*n,lower.tail = T));
      pchet<-pchet+(1-ppois(1,lambda = f_het*n,lower.tail = T));
      if(f_het!=f_hom){
        phom<-phom+(1-ppois(0,lambda = f_hom*n,lower.tail = T));
      }else{
        f_hom<-pbinom(1,2,f_het,lower.tail = F)
        phom<-phom+(1-ppois(0,lambda = f_hom*n,lower.tail = T));
      }
    }
    ghet[as.integer(score/0.05)]=phet
    gchet[as.integer(score/0.05)]=pchet
    ghom[as.integer(score/0.05)]=phom
    
  }
  print(Sys.time())
  eid<-max(which(ghet>0))
  for(i in eid:1){
    het_popscore[i]<-het_popscore[i+1]+ghet[i]*(1-het_popscore[i+1])
    homo_popscore[i]<-homo_popscore[i+1]+ghom[i]*(1-homo_popscore[i+1])
    chet_popscore[i]<-chet_popscore[i+1]+gchet[i]*(1-het_popscore[i+1])+(ghet[i]-gchet[i])*(het_popscore[i+1]-chet_popscore[i+1])
  }
  all_het<-rbind(all_het,c(g,het_popscore))
  all_chet<-rbind(all_chet,c(g,chet_popscore))
  all_hom<-rbind(all_hom,c(g,homo_popscore))
  print(Sys.time()) 
  
}
## merge to all
write.csv(all_het,file = paste("result/",prefix,"_het_popscore_s.csv",sep=""),row.names = F,quote = F)
write.csv(all_chet,file = paste("result/",prefix,"_chet_popscore_s.csv",sep=""),row.names = F,quote = F)
write.csv(all_hom,file = paste("result/",prefix,"_hom_popscore_s.csv",sep=""),row.names = F,quote = F)
#}


############

