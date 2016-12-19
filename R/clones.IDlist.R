## Julia Bischof
## 16-12-2016

#library(doParallel)
#library(parallel)

clones.IDlist<-function(clones.seqID=NULL, summarytab.seqID=NULL){
  if(length(clones.seqID)==0 || !is.vector(clones.seqID)){
    stop("--> Vector containing sequence ID's is missing or empty or is no vector")
  }

  # match clone ID's to sequence ID's
  clones.nr<-cbind(paste("clone",seq(1,length(clones.seqID),1),sep="_"), clones.seqID)
  out<-cbind(unlist(lapply(clones.nr[,2],function(x){strsplit(x,split=", ")[[1]]})),unlist(apply(clones.nr,1,function(x){rep(x[1],(length(gregexpr(", ",x[2])[[1]])+1))})))
  
  # check for duplicated ID's:
  if(length(which(duplicated(out[,1])))>0){
    add.out<-vector()
    dupli.set<-out[which(duplicated(out[,1])==T),1]
    for(i in 1:length(dupli.set)){  
      add.out<-rbind(add.out, c(dupli.set[i],do.call(paste, c(as.list(unique(out[which(out[,1]==dupli.set[i]),2])), sep=", ")))) 
    }
    out<-out[-which(out[,1] %in% dupli.set),]
    out<-rbind(out,add.out)
  }
  
  # order data frame in case of summary.tab != NULL
  if(length(summarytab.seqID)>0){
    out<-out[match(summarytab.seqID,out[,1]),]
    out[,1]<-summarytab.seqID
    out[which(is.na(out[,2])),2]<-"no_clone"
  }else{
    out<-out[order(out[,1]),]
  }
  out<-data.frame(out, row.names = NULL)
  colnames(out)<-c("Sequence_ID","Clone_ID")
  
  return(out)
}