## Julia Bischof
## 16-12-2016

#library(gplots)
#library(doParallel)
#library(parallel)


compare.geneUsage<-function(gene.list=NULL,level=c("subgroup","gene","allele"),abundance=c("relative","absolute"),names=NULL, nrCores=1){  
  levelx<-level
  if(length(gene.list)<2 || is.list(gene.list)==F){
    stop("--> Need a list of at least 2 vectors of gene usage information")
  }
  if(length(level)!=1 || !(level %in% c("subgroup","gene","allele"))){
    stop("--> Gene level (subgroup, gene, allele) is missing")
  }
  if(length(abundance)!=1 || !(abundance %in% c("relative","absolute"))){
    abundance<-"relative"
  }
  if(as.numeric(nrCores)>as.numeric(detectCores())){
    stop(paste("--> nrCores is higher than available number of cores (only ",as.numeric(detectCores())," cores available)",sep=""))
  }
  
  if(is.factor(gene.list[[1]])==T){
    gene.list<-lapply(gene.list,function(x){as.vector(x)})
  }
  
  cl<-makeCluster(nrCores)
  registerDoParallel(cl)

  genes<-gene.list[[1]]
  if(length(grep(" ",genes[which(genes!="")][1]))==0){
    family<-substr(genes[which(genes!="")][1],1,4)
  }else{
    family<-substr(strsplit(genes[which(genes!="")][1],split=" ")[[1]][2],1,4)
  }
  
  genes.all<-unique(unlist(apply(data.frame(unique(unlist(gene.list))),1,function(x){strsplit(x,split=" |,|[.]|/|_|;")})))
  if(length(grep(" ", gene.list[[1]][which(gene.list[[1]]!="")[1]]))>0){
    genes.all<-sort(genes.all[grep(substr(strsplit(gene.list[[1]][which(gene.list[[1]]!="")[1]],split=" ")[[1]][grep("[A-Z]",strsplit(gene.list[[1]][which(gene.list[[1]]!="")[1]],split=" ")[[1]])][2],1,4),genes.all)])
  }else{
    genes.all<-sort(genes.all)
  }
  if(level=="subgroup"){
    genes.all<-sort(unique(unlist(apply(data.frame(genes.all),1,function(x){strsplit(x,split="S|-")[[1]][1]}))))
  }
  if(level=="gene"){
    genes.all<-sort(unique(unlist(apply(data.frame(genes.all),1,function(x){strsplit(x,split="[*]")[[1]][1]}))))
  }
  
    i<-NULL
    temp<-vector()
    genelist<-vector()
    comp.tab<-foreach(i=1:length(gene.list), .combine = "rbind") %dopar% {
      genelist<-unlist(apply(data.frame(gene.list[[i]]),1,function(x){strsplit(x,split=" |,|[.]|;|[|]|_|/")[[1]]}))
      genelist<-sort(genelist[grep(family,genelist)])
      if(level=="subgroup"){
        genelist<-sort(unlist(apply(data.frame(genelist),1,function(x){strsplit(x,split="S|-|[*]")[[1]][1]})))
      }else if(level=="gene"){
        genelist<-sort(unlist(apply(data.frame(genelist),1,function(x){strsplit(x,split="[*]")[[1]][1]})))
      }
      
      temp<-t(data.frame(apply(data.frame(genes.all),1,function(x){length(which(genelist==x))})))
      if(abundance=="relative"){
        temp<-temp/length(genelist)
      }
      return(temp)
    }
  colnames(comp.tab)<-genes.all
  if(length(names)==length(gene.list)){
    rownames(comp.tab)<-names
  }else{
    rownames(comp.tab)<-paste("Sample",1:length(gene.list),sep="")
  }
  
  stopCluster(cl)
  
  return(data.frame(comp.tab,check.names = F))
}
