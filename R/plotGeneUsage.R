## Julia Bischof
## 10-09-2015

plotGeneUsage<-function(geneUsage.tab=NULL,plotFunctionality=FALSE,plotJunctionFr=FALSE,
                              color=c("orange","darkblue","gray"),title=NULL,PDF=NULL,...){
  
if(length(geneUsage.tab)==0){
  stop("--> geneUsage() output table is missing")
}  

if(length(PDF)>0){
      pdf(file = paste(PDF,"_Gene-usage.pdf",sep=""),width = if(ncol(geneUsage.tab$gene_usage)*0.3<7){7}else{ncol(geneUsage.tab$gene_usage)*0.3},height = 7,pointsize = if(ncol(geneUsage.tab$gene_usage)*0.5<16){16}else{ncol(geneUsage.tab$gene_usage)*0.5})
    }
if(max(geneUsage.tab$gene_usage)<=1){
  newnames<-colnames(geneUsage.tab$gene_usage)
  geneUsage.tab$gene_usage<-as.numeric(geneUsage.tab$gene_usage)*100
  barplot(geneUsage.tab$gene_usage,col=color[1],ylab="Percentage",main=title,
          ylim=c(0,max(geneUsage.tab$gene_usage)+10),las=3, names=newnames) 
}else{
  barplot(geneUsage.tab$gene_usage,col=color[1],ylab="Quantity",main=title,
          ylim=c(0,max(geneUsage.tab$gene_usage)+(max(geneUsage.tab$gene_usage)/5)),las=3) 
}
     if(length(PDF)>0){
      dev.off()
    }
    # gene usage vs. functionality
    if(plotFunctionality==TRUE){
      if(length(PDF)>0){
        pdf(file = paste(PDF,"_Gene-usage_vs_Functionality.pdf",sep=""),width = if(ncol(geneUsage.tab$gene_usage_vs_functionality)*0.2<7){7}else{ncol(geneUsage.tab$gene_usage_vs_functionality)*0.2},height = 7,pointsize = if(ncol(geneUsage.tab$gene_usage_vs_functionality)*0.5<16){16}else{ncol(geneUsage.tab$gene_usage_vs_functionality)*0.5})
      }
      barplot(as.matrix(geneUsage.tab$gene_usage_vs_functionality),col=color,ylim=if(max(geneUsage.tab$gene_usage_vs_functionality,na.rm=T)>1){c(0,max(geneUsage.tab$gene_usage_vs_functionality)+(max(geneUsage.tab$gene_usage_vs_functionality)/5))}else{c(0,1.6)},main=title,
              ylab=if(max(geneUsage.tab$gene_usage_vs_functionality,na.rm=T)<=1){""}else{"Quantity"},axes=if(max(geneUsage.tab$gene_usage_vs_functionality,na.rm=T)<=1){F}else{T},las=3)
      if(max(geneUsage.tab$gene_usage_vs_functionality,na.rm=T)<=1){
        axis(2,at=seq(0,1,0.25),seq(0,100,25))
        mtext(side = 2, line = 2.5, at = 0.5, text = "Percentage")
      }
      legend("topright",col=color,c("productive","unproductive","unknown"),title="Functionality",pch=15,y.intersp=1.2)
      if(length(PDF)>0){
        dev.off()
      }
    }
    
    # gene usage vs. JUNCTION_frame
    if(plotJunctionFr==TRUE){
      if(length(PDF)>0){
        pdf(file = paste(PDF,"_Gene-usage_vs_Junction-frame.pdf",sep=""),width = if(ncol(geneUsage.tab$gene_usage_vs_junction_frame)*0.2<7){7}else{ncol(geneUsage.tab$gene_usage_vs_junction_frame)*0.2},height = 7,pointsize = if(ncol(geneUsage.tab$gene_usage_vs_junction_frame)*0.5<16){16}else{ncol(geneUsage.tab$gene_usage_vs_junction_frame)*0.5})
      }
      barplot(as.matrix(geneUsage.tab$gene_usage_vs_junction_frame),col=color,ylim=if(max(geneUsage.tab$geneUsage.tab$gene_usage_vs_junction_frame,na.rm=T)<=1){c(0,max(geneUsage.tab$gene_usage_vs_junction_frame)+(max(geneUsage.tab$gene_usage_vs_junction_frame)/5))}else{c(0,1.6)},main=title,
              ylab=if(max(geneUsage.tab$geneUsage.tab$gene_usage_vs_junction_frame,na.rm=T)<=1){""}else{"Quantity"},axes=if(max(geneUsage.tab$geneUsage.tab$gene_usage_vs_junction_frame,na.rm=T)<=1){F}else{T},las=3)
      if(max(geneUsage.tab$geneUsage.tab$gene_usage_vs_junction_frame,na.rm=T)<=1){
        axis(2,at=seq(0,1,0.25),seq(0,100,25))
        mtext(side = 2, line = 2.5, at = 0.5, text = "Percentage")
      }
      legend("topright",col=color,c("in-frame","out-of-frame","unknown"),title="Junction frame",pch=15,y.intersp=1.2)
      if(length(PDF)>0){
        dev.off()
      }
    }
}


