#setwd("~/") 
setwd ("/home/local/ARCS/nz2274/PSAP_Setup/")
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
#finputs="PSAP_Setup/test/2_CADD_ExAC_.txt"
args = commandArgs(trailingOnly=TRUE)
finput= args[1]  #"data/temp/gencode.v27.simulation.table.out.txt" ## arg input
#for(finput in finputs){
  #finput="PSAP_Setup/test/2_CADD13_gnomAD_.txt" #"PSAP_Setup/Data/source/CADD13_gnomAD_0.txt"
  prefixs<-unlist(strsplit(finput,"/",fixed=T))
  prefix<-prefixs[length(prefixs)]
  dat<-read.table(finput,header = F,stringsAsFactors = F)
  names(dat)<-c("Gene"  ,"het_freq"        ,"homo_freq"       ,"CADD"    ,"Ntime")
  simu<-1000000
  
  all_het<-c()
  all_hom<-c()
  all_chet<-c()
  genes<- unique(dat$Gene)
  #genes<-"SOX17"
  for(g in genes){
    # vecs<-rep(0,simu)
    if(g=="Gene"){next;}
    print(g)
    het<-rep(0,simu);
    homo<-rep(0,simu)
    chet<-rep(0,simu)
    gdat<-dat[which(dat$Gene==g),]
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
       vec<-rbinom(simu,1,phet)
       index<-intersect(which(vec==1),which(het!=0))
       chet[index]<-het[index]
       het[which(vec==1)]<-score
       vec<-rbinom(simu,1,pchet)
       chet[which(vec==1)]<-score
       vec<-rbinom(simu,1,phom)
       homo[which(vec==1)]<-score
       
    }
    ### end a gene 
    print(Sys.time()) 
    het_popscore<-popscore(het)
    homo_popscore<-popscore(homo)
    chet_popscore<-popscore(chet)
    print(Sys.time()) 
    fhet<- file(paste("result/",prefix,"_het_popscore_a.csv",sep=""),'a')
    write.table(data.frame(matrix(c(g,het_popscore),nrow=1),stringsAsFactors = F),file = fhet,row.names = F,quote = F, col.names=F)
    fchet<- file(paste("result/",prefix,"_chet_popscore_a.csv",sep=""),'a')
    close(fhet)
    
    write.table(data.frame(matrix(c(g,chet_popscore),nrow=1),stringsAsFactors = F),file = fchet,row.names = F,quote = F, col.names=F)
    close(fchet)
    fhom<- file(paste("result/",prefix,"_hom_popscore_a.csv",sep=""),'a')
    write.table((data.frame(matrix(c(g,homo_popscore),nrow=1),stringsAsFactors = F)),file = fhom,row.names = F,quote = F, col.names=F)
    close(fhom)
    print(Sys.time()) 
    all_het<-rbind(all_het,c(g,het_popscore))
    all_chet<-rbind(all_chet,c(g,chet_popscore))
    all_hom<-rbind(all_hom,c(g,homo_popscore))
    print(Sys.time()) 
    
  }
 ## merge to all
  write.csv(all_het,file = paste("result/",prefix,"_het_popscore.csv",sep=""),row.names = F,quote = F)
  write.csv(all_chet,file = paste("result/",prefix,"_chet_popscore.csv",sep=""),row.names = F,quote = F)
  write.csv(all_hom,file = paste("result/",prefix,"_hom_popscore.csv",sep=""),row.names = F,quote = F)
#}
  
 
  ############
  
