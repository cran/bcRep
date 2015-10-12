## Julia Bischof
## 10-09-2015

plotAADistribution<-function(aaDistribution.tab=NULL,plotSeqN=FALSE,PDF=NULL,...){
    if(length(aaDistribution.tab)==0){
      stop("--> Amino acid distribution tab is missing or empty")
    }
    color<-c("lightblue2","red","lightgreen","lightgray","yellow","orange","purple4","cyan","darkblue","darkred",
             "darkgreen","gray33","blueviolet","aquamarine","bisque4","violetred","deepskyblue4","black","olivedrab","olivedrab3",
             "orchid1")
    if(length(PDF)>0){
      pdf(file = paste(PDF,"_Amino-acid-distribution.pdf",sep=""),
          width = ceiling(sqrt(length(aaDistribution.tab$Amino_acid_distribution)+1))*5,
          height = floor(sqrt(length(aaDistribution.tab$Amino_acid_distribution)+1))*6,
          pointsize = if(sqrt(length(aaDistribution.tab$Amino_acid_distribution)+1)*3.5<26){26}else{sqrt(length(aaDistribution.tab$Amino_acid_distribution)+1)*3.5},onefile = F)
    }
    par(mfrow=if(ceiling(sqrt(length(aaDistribution.tab$Amino_acid_distribution)+1))*floor(sqrt(length(aaDistribution.tab$Amino_acid_distribution)+1))>length(aaDistribution.tab$Amino_acid_distribution)+1){c(ceiling(sqrt(length(aaDistribution.tab$Amino_acid_distribution)+1)),floor(sqrt(length(aaDistribution.tab$Amino_acid_distribution)+1)))}else{c(ceiling(sqrt(length(aaDistribution.tab$Amino_acid_distribution)+1)),ceiling(sqrt(length(aaDistribution.tab$Amino_acid_distribution)+1)))},
        mar=c(4,4,3,1),oma=c(2,2,7,3))
    for(i in 1:length(aaDistribution.tab$Amino_acid_distribution)){
      barplot(aaDistribution.tab$Amino_acid_distribution[[i]],col=color,axes=F,ylab="Percentage",
              names=seq(1,ncol(aaDistribution.tab$Amino_acid_distribution[[i]]),1),xlab="Position",
              main=names(aaDistribution.tab$Amino_acid_distribution)[i])
      axis(2,at=seq(0,1,0.25),seq(0,100,25))
    }
    plot(x=seq(0.5,6.5,1),y=rep(5,7),pch=18,cex=2.3,col=color[1:7],axes=F,xlab="",ylab="",xlim=c(0.2,7.2),ylim=c(0,6))
    points(x=seq(0.5,6.5,1),y=rep(3,7),pch=18,cex=2.3,col=color[8:14])
    points(x=seq(0.5,6.5,1),y=rep(1,7),pch=18,cex=2.3,col=color[15:21])
    text(x=seq(1,7,1),y=5,rownames(aaDistribution.tab$Amino_acid_distribution[[1]])[1:7],cex=1.1)
    text(x=seq(1,7,1),y=3,rownames(aaDistribution.tab$Amino_acid_distribution[[1]])[8:14],cex=1.1)
    text(x=seq(1,7,1),y=1,rownames(aaDistribution.tab$Amino_acid_distribution[[1]])[15:21],cex=1.1)
    
    title("Amino acid distribution",outer=T)
    if(length(PDF)>0){
      dev.off()
    }
  if(plotSeqN==T){    
    if(nrow(aaDistribution.tab$Number_of_sequences_per_length)==0){
      stop("--> 'Number of sequences' tab is missing or empty")
    }
    if(length(PDF)>0){
      pdf(file = paste(PDF,"_Number-of-sequences.pdf",sep=""),
          width = sqrt(nrow(aaDistribution.tab$Number_of_sequences_per_length))*5,
          height = sqrt(nrow(aaDistribution.tab$Number_of_sequences_per_length))*6, 
          pointsize = if(sqrt(nrow(aaDistribution.tab$Number_of_sequences_per_length))*3.5<26){26}else{sqrt(nrow(aaDistribution.tab$Number_of_sequences_per_length))*3.5},onefile = F)
    }
    par(mfrow=c(ceiling(sqrt(nrow(aaDistribution.tab$Number_of_sequences_per_length))),ceiling(sqrt(nrow(aaDistribution.tab$Number_of_sequences_per_length)))),
        mar=c(4,4,4,1),oma=c(2,2,7,3))
    for(i in 1:nrow(aaDistribution.tab$Number_of_sequences_per_length)){
      plot(1,1,col="white",ylab="Number of sequences",xlab=" ",main=rownames(aaDistribution.tab$Number_of_sequences_per_length)[i],
           ylim=c(0,max(aaDistribution.tab$Number_of_sequences_per_length[i],na.rm=T)),xaxt="n",axes=if(max(aaDistribution.tab$Number_of_sequences_per_length[i],na.rm=T)<=5){F}else{T})
      if(max(aaDistribution.tab$Number_of_sequences_per_length[i],na.rm=T)<=5){
        axis(2,at=seq(0,max(aaDistribution.tab$Number_of_sequences_per_length[i],na.rm=T),1),seq(0,max(aaDistribution.tab$Number_of_sequences_per_length[i],na.rm=T),1))
      }
      abline(h=aaDistribution.tab$Number_of_sequences_per_length[i,1],col="red")
    }
    title("Number of sequences per length",outer=T)
    if(length(PDF)>0){
      dev.off()
    }
  }
}


