## Julia Bischof
## 10-09-2015

plotClonesCopyNumber<-function(copyNumber=NULL, color="gray", title=NULL, PDF=NULL,...){
  if(length(copyNumber)==0){
    stop("--> Copy number vector is missing")
  }else{
    copyNumber<-as.numeric(copyNumber)
  }
  
  if(length(PDF)>0){
    pdf(file = paste(PDF,"_Clone-copy-number.pdf",sep=""),width = 7,height = 7,pointsize = 18)
  }
  boxplot(copyNumber,col=color,xlab=" ",ylab="Copy number",main=title)
  if(length(PDF)>0){
    dev.off()
  }
}


