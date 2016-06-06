## Julia Bischof
## 2016-02-24

plotSequencesMutationBase<-function(mutationBaseTab=NULL, plotEnvironment = FALSE, plotMutation = TRUE, 
                                    colHeatmap = c("white","darkblue"),title=NULL, PDF=NULL){
  if(length(mutationBaseTab)==0){
    stop("--> Output of sequences.mutation.base() is missing")
  }
  
  if(plotEnvironment == T){
    if(is.list(mutationBaseTab)==T){
      tab1<-mutationBaseTab$Environment
    }else{
      tab1<-mutationBaseTab
    }
    if(length(PDF)>0){
      pdf(paste(PDF,"_Base-mutation_environment.pdf",sep=""), width = 7, height=7, pointsize = 10)
    }
    par(mfrow=c(2,2),mar=c(5,5,4,2),oma=c(2,2,4,2))
    
    for(i in c("a_to","c_to","g_to","t_to")){
      barplot(as.matrix(tab1[grep(i, rownames(tab1)),]), col=c("darkgreen","darkred","darkgray","darkblue"), yaxt="n", ylim=c(0,1), xlim=c(0,10),
              names.arg = c("-3","-2","-1","0","+1","+2","+3"), xlab="Position", ylab="Percentage", main=paste(strsplit(i,split="_")[[1]][1]," -> n",sep=""))
      axis(2,at = seq(0,1,0.25), seq(0,100,25))
      legend("right", legend = c("t", "g","c","a"), col=c("darkblue", "darkgray","darkred","darkgreen"), y.intersp = 1.2, cex=1.1, pt.cex=2, pch=15)
    }
    title(title, outer=T, cex=1.8)
    
    if(length(PDF)>0){
      dev.off()
    }
  }
  
  if(plotMutation==T){
    if(is.list(mutationBaseTab)==T){
      tab2<-mutationBaseTab$Mutated_position
    }else{
      tab2<-mutationBaseTab
    }
    if(length(PDF)>0){
      pdf(paste(PDF,"_Base-mutation_mutated-position.pdf",sep=""), width = 9, height=7, pointsize = 10)
    }
    par(mfrow=c(1,1),oma=c(3,1,1,3), las=1)   
    heatmap.2(as.matrix(tab2),col=colorRampPalette(colHeatmap), Rowv=F,Colv=F,key=T,density.info="none", labRow = gsub("_"," ",rownames(tab2)), labCol = gsub("_"," ",colnames(tab2)),
              trace="none",cexCol=1.5,cexRow=1.5,dendrogram="none", keysize = 1.5,key.xlab="proportion")
    title(title, outer=T, cex=1.8)
    
    if(length(PDF)>0){
      dev.off()
    }
  }
  
}