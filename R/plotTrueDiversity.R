## Julia Bischof
## 10-09-2015

plotTrueDiversity<-function(trueDiversity.tab=NULL,color="black",PDF=NULL,...){
  if(length(trueDiversity.tab)==0){
    stop("--> True diversity tab is missing or empty")
  }
  if(length(PDF)>0){
    pdf(file = paste(PDF,"_True-diversity_q",trueDiversity.tab$True_diversity_order,".pdf",sep=""),
        width = ceiling(sqrt(length(trueDiversity.tab$True_diversity)))*5,
        height = floor(sqrt(length(trueDiversity.tab$True_diversity)))*6,
        pointsize = if(sqrt(length(trueDiversity.tab$True_diversity)+1)*4<26){26}else{sqrt(length(trueDiversity.tab$True_diversity)+1)*4},onefile = F)
  }
  par(mfrow=if(ceiling(sqrt(length(trueDiversity.tab$True_diversity)))*floor(sqrt(length(trueDiversity.tab$True_diversity)))>length(trueDiversity.tab$True_diversity)){c(ceiling(sqrt(length(trueDiversity.tab$True_diversity))),floor(sqrt(length(trueDiversity.tab$True_diversity))))}else{c(ceiling(sqrt(length(trueDiversity.tab$True_diversity))),ceiling(sqrt(length(trueDiversity.tab$True_diversity))))},
      ps=22,mar=c(5,5,5,3),oma=c(2,2,5,2))
  for(i in 1:length(trueDiversity.tab$True_diversity)){
    plot(as.numeric(trueDiversity.tab$True_diversity[[i]][1,]),col=color,type="b",
         ylab=if(trueDiversity.tab$True_diversity_order==0){"Richness"}else{"Diversity"},
         xlab="Position",axes=F,ylim=c(0,max(unlist(trueDiversity.tab$True_diversity),na.rm=T)+5),
         main=names(trueDiversity.tab$True_diversity)[i])
    axis(1,at=seq(1,ncol(trueDiversity.tab$True_diversity[[i]]),1),seq(1,ncol(trueDiversity.tab$True_diversity[[i]]),1))
    axis(2,at=seq(0,max(unlist(trueDiversity.tab$True_diversity),na.rm=T)+5,5),seq(0,max(unlist(trueDiversity.tab$True_diversity),na.rm=T)+5,5))
  }
  title(paste("True diversity, order ",as.numeric(trueDiversity.tab$True_diversity_order),sep=""),outer=T)
  if(length(PDF)>0){
    dev.off()
  }
}


