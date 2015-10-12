## Julia Bischof
## 10-09-2015

library(doParallel)
library(parallel)


clones.shared<-function(clones.tab=NULL,identity=0.85,useJ=TRUE,dispD=TRUE,dispCDR3aa=TRUE,dispCDR3nt=FALSE,dispFunctionality.list=FALSE,dispFunctionality.ratio=FALSE,
                        dispJunctionFr.list=FALSE,dispJunctionFr.ratio=FALSE,dispTotalSeq=FALSE,nrCores=1){
  if(length(clones.tab)==0){
    stop("--> Clone file is missing or empty")
  }
  if(as.numeric(nrCores)>as.numeric(detectCores())){
    stop(paste("--> nrCores is higher than available number of cores (only ",as.numeric(detectCores())," cores available)",sep=""))
  }
  
  clones.uniV<-unique(unlist(apply(data.frame(clones.tab$V_gene),1,function(x){strsplit(x,split=" |, ")[[1]]})))
  #clones.uniV<-unique(unlist(apply(data.frame(unlist(clones.V)),1,function(x){strsplit(x,split=" ")[[1]]})))
  clones.uniV<-clones.uniV[grep(substr(clones.tab$V_gene[1],1,4),clones.uniV)]
  clones.uniV.gene<-unique(apply(data.frame(clones.uniV),1,function(x){strsplit(x,split="[*]")[[1]][1]}))
  
  registerDoParallel(cores=nrCores)
  
  sharedclone.temp<-vector()
  clonelist<-foreach(i=1:length(clones.uniV.gene), .export="sharedclone.temp") %dopar%{
    cat(paste(i,"/",length(clones.uniV.gene),"\n"))
    newCDR3length=NULL
    if(length(grep(paste(clones.uniV.gene[i],"$",sep=""),clones.tab$V_gene))>1 && length(unique(clones.tab[,1][grep(clones.uniV.gene[i],clones.tab$V_gene)]))>1){
      sameV<-grep(paste(clones.uniV.gene[i],"$",sep=""),clones.tab$V_gene)
      
      if(useJ==T){
        clones.J<-apply(data.frame(clones.tab$J_gene[grep(clones.uniV.gene[i],clones.tab$V_gene)]),1,function(x){strsplit(x,split=", ")[[1]]})
        clones.uniJ<-unique(unlist(apply(data.frame(unlist(clones.J)),1,function(x){strsplit(x,split=" ")[[1]]})))
        clones.uniJ<-clones.uniJ[grep("IGHJ",clones.uniJ)]
        clones.uniJ.gene<-unique(apply(data.frame(clones.uniJ),1,function(x){strsplit(x,split="[*]")[[1]][1]}))
        if(length(clones.uniJ.gene)>1){ # diff. J's
          sameVJ<-vector()
          CDR3length<-vector()
          for(j in 1:length(clones.uniJ.gene)){
            sameVJ<-grep(clones.uniJ.gene[j],clones.tab$J_gene[sameV])
            CDR3length<-clones.tab$CDR3_length_AA[sameV[sameVJ]]
            if(length(unique(clones.tab[,1][sameV[sameVJ]]))>1){
              if(length(unique(CDR3length))==1){ # all same CDR3 length
                newCDR3length<-as.numeric(unique(CDR3length))
                tr<-floor(newCDR3length*(1-as.numeric(identity)))
                CDR3.adist.temp<-unlist(apply(data.frame(clones.tab$unique_CDR3_sequences_AA[sameV[sameVJ]]),1,function(x){strsplit(x,split=", |; |_")[[1]]}))
                CDR3.adist<-unique(unlist(apply(data.frame(clones.tab$unique_CDR3_sequences_AA[sameV[sameVJ]]),1,function(x){strsplit(x,split=", |; |_")[[1]]})))
                ind.adist<-unlist(apply(data.frame(clones.tab[sameV[sameVJ],]),1,function(x){rep(x[1],length(strsplit(x[2],split=", |; |_")[[1]]))}))
                a.dist<-adist(CDR3.adist)#adist(sameCDR3l)  
                diag(a.dist)<-NA
                a.dist[upper.tri(a.dist)]<-NA
                for(a in 1:nrow(a.dist)){
                  clones.tab.new<-vector()
                  grepCDR3<-grep(gsub("[*]","-",do.call(paste, c(as.list(as.character(unique(CDR3.adist[c(a,which(a.dist[a,]<=tr))]))), sep="|"))),gsub("[*]","-",CDR3.adist.temp),perl=T)
                  if(length(which(a.dist[a,]<=tr))>0 && length(unique(ind.adist[grepCDR3]))>1){
                    clones.tab.new<-clones.tab[grep(gsub("[*]","-",do.call(paste, c(as.list(as.character(unique(CDR3.adist[c(a,which(a.dist[a,]<=tr))]))), sep="|"))),gsub("[*]","-",clones.tab$unique_CDR3_sequences_AA),perl=T),]
                    
                    if(length(clones.tab.new)>0 && length(unique(clones.tab.new[,1]))>1){
                      sharedclone.temp<-rbind(sharedclone.temp,c(length(unique(clones.tab.new[,1])),
                                                                 do.call(paste, c(as.list(as.character(unique(clones.tab.new[,1]))), sep="; ")),
                                                                 nchar(CDR3.adist[a]),
                                                                 do.call(paste, c(as.list(as.character(unique(CDR3.adist[c(a,which(a.dist[a,]<=tr))]))), sep="; ")),
                                                                 length(unique(CDR3.adist[c(a,which(a.dist[a,]<=tr))])),
                                                                 do.call(paste, c(as.list(as.character(clones.tab.new$unique_CDR3_sequences_AA)), sep="; ")),
                                                                 do.call(paste, c(as.list(as.character(clones.tab.new$sequence_count_per_CDR3)), sep="; ")),
                                                                 clones.uniV.gene[i],
                                                                 if(useJ==T && length(clones.uniJ.gene[j])>0){do.call(paste, c(as.list(unique(clones.uniJ.gene[j])), sep=", "))}else if(useJ==T && length(clones.uniJ.gene[j])==0){"no J"}else{do.call(paste, c(as.list(unique(clones.tab.new$J_gene)), sep="; "))},
                                                                 
                                                                 if(dispD==T){do.call(paste, c(as.list(unique(clones.tab.new$D_gene)), sep="; "))},
                                                                 if(dispCDR3aa==T){do.call(paste, c(as.list(unique(clones.tab.new$all_CDR3_sequences_AA)), sep="; "))},
                                                                 if(dispCDR3nt==T){do.call(paste, c(as.list(unique(clones.tab.new$all_CDR3_sequences_nt)), sep="; "))},
                                                                 if(dispFunctionality.list==T){do.call(paste, c(as.list(unique(clones.tab.new$Functionality_all_sequences)), sep="; "))},
                                                                 if(dispFunctionality.ratio==T){length(unlist(gregexpr("^productive",clones.tab.new$Functionality_all_sequences))[which(unlist(gregexpr("^productive",clones.tab.new$Functionality_all_sequences))!=(-1))])/length(unlist(apply(data.frame(clones.tab.new$Functionality_all_sequences),1,function(x){strsplit(x,split=", ")[[1]]})))},
                                                                 if(dispFunctionality.ratio==T){length(unlist(gregexpr("^unproductive",clones.tab.new$Functionality_all_sequences))[which(unlist(gregexpr("^unproductive",clones.tab.new$Functionality_all_sequences))!=(-1))])/length(unlist(apply(data.frame(clones.tab.new$Functionality_all_sequences),1,function(x){strsplit(x,split=", ")[[1]]})))},
                                                                 if(dispJunctionFr.list==T){do.call(paste, c(as.list(clones.tab.new$Junction_frame_all_sequences), sep="; "))},
                                                                 if(dispJunctionFr.ratio==T){length(unlist(gregexpr("in-frame",clones.tab.new$Junction_frame_all_sequences))[which(unlist(gregexpr("in-frame",clones.tab.new$Junction_frame_all_sequences))!=(-1))])/length(unlist(apply(data.frame(clones.tab.new$Junction_frame_all_sequences),1,function(x){strsplit(x,split=", ")[[1]]})))},
                                                                 if(dispJunctionFr.ratio==T){length(unlist(gregexpr("out-of-frame",clones.tab.new$Junction_frame_all_sequences))[which(unlist(gregexpr("out-of-frame",clones.tab.new$Junction_frame_all_sequences))!=(-1))])/length(unlist(apply(data.frame(clones.tab.new$Junction_frame_all_sequences),1,function(x){strsplit(x,split=", ")[[1]]})))},
                                                                 if(dispTotalSeq==T){do.call(paste, c(as.list(unique(clones.tab.new$Total_sequences_nt)), sep="; "))}))
                    }
                  }
                }                        
              }else{ #different CDR3 length
                for(k in 1:length(unique(CDR3length))){
                  newCDR3length<-as.numeric(unique(CDR3length)[k])
                  tr<-floor(newCDR3length*(1-as.numeric(identity)))                                
                  sameCDR3l<-clones.tab$unique_CDR3_sequences_AA[sameV[sameVJ]][which(CDR3length==unique(CDR3length)[k])]
                  if(length(sameCDR3l)>1){
                    CDR3.adist.temp<-unlist(apply(data.frame(clones.tab$unique_CDR3_sequences_AA[sameV[sameVJ]][which(CDR3length==unique(CDR3length)[k])]),1,function(x){strsplit(x,split=", |; |_")[[1]]}))
                    CDR3.adist<-unique(unlist(apply(data.frame(clones.tab$unique_CDR3_sequences_AA[sameV[sameVJ]][which(CDR3length==unique(CDR3length)[k])]),1,function(x){strsplit(x,split=", |; |_")[[1]]})))
                    ind.adist<-unlist(apply(data.frame(clones.tab[sameV[sameVJ][which(CDR3length==unique(CDR3length)[k])],]),1,function(x){rep(x[1],length(strsplit(x[2],split=", |; |_")[[1]]))}))
                    a.dist<-adist(CDR3.adist)#adist(sameCDR3l)  
                    diag(a.dist)<-NA
                    a.dist[upper.tri(a.dist)]<-NA
                    for(a in 1:nrow(a.dist)){
                      clones.tab.new<-vector()
                      grepCDR3<-grep(gsub("[*]","-",do.call(paste, c(as.list(as.character(unique(CDR3.adist[c(a,which(a.dist[a,]<=tr))]))), sep="|"))),gsub("[*]","-",CDR3.adist.temp),perl=T)
                      if(length(which(a.dist[a,]<=tr))>0 && length(unique(ind.adist[grepCDR3]))>1){
                        clones.tab.new<-clones.tab[grep(gsub("[*]","-",do.call(paste, c(as.list(as.character(unique(CDR3.adist[c(a,which(a.dist[a,]<=tr))]))), sep="|"))),gsub("[*]","-",clones.tab$unique_CDR3_sequences_AA),perl=T),]
                        if(length(clones.tab.new)>0 && length(unique(clones.tab.new[,1]))>1){
                          sharedclone.temp<-rbind(sharedclone.temp,c(length(unique(clones.tab.new[,1])),
                                                                     do.call(paste, c(as.list(as.character(unique(clones.tab.new[,1]))), sep="; ")),
                                                                     nchar(CDR3.adist[a]),
                                                                     do.call(paste, c(as.list(as.character(unique(CDR3.adist[c(a,which(a.dist[a,]<=tr))]))), sep="; ")),
                                                                     length(unique(CDR3.adist[c(a,which(a.dist[a,]<=tr))])),
                                                                     do.call(paste, c(as.list(as.character(clones.tab.new$unique_CDR3_sequences_AA)), sep="; ")),
                                                                     do.call(paste, c(as.list(as.character(clones.tab.new$sequence_count_per_CDR3)), sep="; ")),
                                                                     clones.uniV.gene[i],
                                                                     if(useJ==T && length(clones.uniJ.gene[j])>0){do.call(paste, c(as.list(unique(clones.uniJ.gene[j])), sep=", "))}else if(useJ==T && length(clones.uniJ.gene[j])==0){"no J"}else{do.call(paste, c(as.list(unique(clones.tab.new$J_gene)), sep="; "))},
                                                                     
                                                                     if(dispD==T){do.call(paste, c(as.list(unique(clones.tab.new$D_gene)), sep="; "))},
                                                                     if(dispCDR3aa==T){do.call(paste, c(as.list(unique(clones.tab.new$all_CDR3_sequences_AA)), sep="; "))},
                                                                     if(dispCDR3nt==T){do.call(paste, c(as.list(unique(clones.tab.new$all_CDR3_sequences_nt)), sep="; "))},
                                                                     if(dispFunctionality.list==T){do.call(paste, c(as.list(unique(clones.tab.new$Functionality_all_sequences)), sep="; "))},
                                                                     if(dispFunctionality.ratio==T){length(unlist(gregexpr("^productive",clones.tab.new$Functionality_all_sequences))[which(unlist(gregexpr("^productive",clones.tab.new$Functionality_all_sequences))!=(-1))])/length(unlist(apply(data.frame(clones.tab.new$Functionality_all_sequences),1,function(x){strsplit(x,split=", ")[[1]]})))},
                                                                     if(dispFunctionality.ratio==T){length(unlist(gregexpr("^unproductive",clones.tab.new$Functionality_all_sequences))[which(unlist(gregexpr("^unproductive",clones.tab.new$Functionality_all_sequences))!=(-1))])/length(unlist(apply(data.frame(clones.tab.new$Functionality_all_sequences),1,function(x){strsplit(x,split=", ")[[1]]})))},
                                                                     if(dispJunctionFr.list==T){do.call(paste, c(as.list(clones.tab.new$Junction_frame_all_sequences), sep="; "))},
                                                                     if(dispJunctionFr.ratio==T){length(unlist(gregexpr("in-frame",clones.tab.new$Junction_frame_all_sequences))[which(unlist(gregexpr("in-frame",clones.tab.new$Junction_frame_all_sequences))!=(-1))])/length(unlist(apply(data.frame(clones.tab.new$Junction_frame_all_sequences),1,function(x){strsplit(x,split=", ")[[1]]})))},
                                                                     if(dispJunctionFr.ratio==T){length(unlist(gregexpr("out-of-frame",clones.tab.new$Junction_frame_all_sequences))[which(unlist(gregexpr("out-of-frame",clones.tab.new$Junction_frame_all_sequences))!=(-1))])/length(unlist(apply(data.frame(clones.tab.new$Junction_frame_all_sequences),1,function(x){strsplit(x,split=", ")[[1]]})))},
                                                                     if(dispTotalSeq==T){do.call(paste, c(as.list(unique(clones.tab.new$Total_sequences_nt)), sep="; "))}))
                        }
                      }
                    }
                  }
                }
              }       
            }
          }
        }else{ # same J/no J
          CDR3length<-clones.tab$CDR3_length_AA[sameV]
          if(length(clones.tab[,1][sameV])>1){
            if(length(unique(CDR3length))==1){ # all same CDR3 length
              newCDR3length<-as.numeric(unique(CDR3length))
              tr<-floor(newCDR3length*(1-as.numeric(identity)))
              CDR3.adist.temp<-unlist(apply(data.frame(clones.tab$unique_CDR3_sequences_AA[sameV]),1,function(x){strsplit(x,split=", |; |_")[[1]]}))
              CDR3.adist<-unique(unlist(apply(data.frame(clones.tab$unique_CDR3_sequences_AA[sameV]),1,function(x){strsplit(x,split=", |; |_")[[1]]})))
              ind.adist<-unlist(apply(data.frame(clones.tab[sameV,]),1,function(x){rep(x[1],length(strsplit(x[2],split=", |; |_")[[1]]))}))
              a.dist<-adist(CDR3.adist)#adist(sameCDR3l)  
              diag(a.dist)<-NA
              a.dist[upper.tri(a.dist)]<-NA
              for(a in 1:nrow(a.dist)){
                clones.tab.new<-vector()
                grepCDR3<-grep(gsub("[*]","-",do.call(paste, c(as.list(as.character(unique(CDR3.adist[c(a,which(a.dist[a,]<=tr))]))), sep="|"))),gsub("[*]","-",CDR3.adist.temp),perl=T)
                if(length(which(a.dist[a,]<=tr))>0 && length(unique(ind.adist[grepCDR3]))>1){
                  clones.tab.new<-clones.tab[grep(gsub("[*]","-",do.call(paste, c(as.list(as.character(unique(CDR3.adist[c(a,which(a.dist[a,]<=tr))]))), sep="|"))),gsub("[*]","-",clones.tab$unique_CDR3_sequences_AA),perl=T),]
                  
                  if(length(clones.tab.new)>0 && length(unique(clones.tab.new[,1]))>1){
                    sharedclone.temp<-rbind(sharedclone.temp,c(length(unique(clones.tab.new[,1])),
                                                               do.call(paste, c(as.list(as.character(unique(clones.tab.new[,1]))), sep="; ")),
                                                               nchar(CDR3.adist[a]),
                                                               do.call(paste, c(as.list(as.character(unique(CDR3.adist[c(a,which(a.dist[a,]<=tr))]))), sep="; ")),
                                                               length(unique(CDR3.adist[c(a,which(a.dist[a,]<=tr))])),
                                                               do.call(paste, c(as.list(as.character(clones.tab.new$unique_CDR3_sequences_AA)), sep="; ")),
                                                               do.call(paste, c(as.list(as.character(clones.tab.new$sequence_count_per_CDR3)), sep="; ")),
                                                               clones.uniV.gene[i],
                                                               if(useJ==T && length(clones.uniJ.gene[j])>0){do.call(paste, c(as.list(unique(clones.uniJ.gene[j])), sep=", "))}else if(useJ==T && length(clones.uniJ.gene[j])==0){"no J"}else{do.call(paste, c(as.list(unique(clones.tab.new$J_gene)), sep="; "))},
                                                               
                                                               if(dispD==T){do.call(paste, c(as.list(unique(clones.tab.new$D_gene)), sep="; "))},
                                                               if(dispCDR3aa==T){do.call(paste, c(as.list(unique(clones.tab.new$all_CDR3_sequences_AA)), sep="; "))},
                                                               if(dispCDR3nt==T){do.call(paste, c(as.list(unique(clones.tab.new$all_CDR3_sequences_nt)), sep="; "))},
                                                               if(dispFunctionality.list==T){do.call(paste, c(as.list(unique(clones.tab.new$Functionality_all_sequences)), sep="; "))},
                                                               if(dispFunctionality.ratio==T){length(unlist(gregexpr("^productive",clones.tab.new$Functionality_all_sequences))[which(unlist(gregexpr("^productive",clones.tab.new$Functionality_all_sequences))!=(-1))])/length(unlist(apply(data.frame(clones.tab.new$Functionality_all_sequences),1,function(x){strsplit(x,split=", ")[[1]]})))},
                                                               if(dispFunctionality.ratio==T){length(unlist(gregexpr("^unproductive",clones.tab.new$Functionality_all_sequences))[which(unlist(gregexpr("^unproductive",clones.tab.new$Functionality_all_sequences))!=(-1))])/length(unlist(apply(data.frame(clones.tab.new$Functionality_all_sequences),1,function(x){strsplit(x,split=", ")[[1]]})))},
                                                               if(dispJunctionFr.list==T){do.call(paste, c(as.list(clones.tab.new$Junction_frame_all_sequences), sep="; "))},
                                                               if(dispJunctionFr.ratio==T){length(unlist(gregexpr("in-frame",clones.tab.new$Junction_frame_all_sequences))[which(unlist(gregexpr("in-frame",clones.tab.new$Junction_frame_all_sequences))!=(-1))])/length(unlist(apply(data.frame(clones.tab.new$Junction_frame_all_sequences),1,function(x){strsplit(x,split=", ")[[1]]})))},
                                                               if(dispJunctionFr.ratio==T){length(unlist(gregexpr("out-of-frame",clones.tab.new$Junction_frame_all_sequences))[which(unlist(gregexpr("out-of-frame",clones.tab.new$Junction_frame_all_sequences))!=(-1))])/length(unlist(apply(data.frame(clones.tab.new$Junction_frame_all_sequences),1,function(x){strsplit(x,split=", ")[[1]]})))},
                                                               if(dispTotalSeq==T){do.call(paste, c(as.list(unique(clones.tab.new$Total_sequences_nt)), sep="; "))}))
                  }
                }
              }
            }else{ # diff CDR3
              for(k in 1:length(unique(CDR3length))){
                newCDR3length<-as.numeric(unique(CDR3length)[k])
                tr<-floor(newCDR3length*(1-as.numeric(identity)))                                
                sameCDR3l<-clones.tab$unique_CDR3_sequences_AA[sameV][which(CDR3length==unique(CDR3length)[k])]
                if(length(sameCDR3l)>1){
                  CDR3.adist.temp<-unlist(apply(data.frame(clones.tab$unique_CDR3_sequences_AA[sameV[which(CDR3length==unique(CDR3length)[k])]]),1,function(x){strsplit(x,split=", |; |_")[[1]]}))
                  CDR3.adist<-unique(unlist(apply(data.frame(clones.tab$unique_CDR3_sequences_AA[sameV[which(CDR3length==unique(CDR3length)[k])]]),1,function(x){strsplit(x,split=", |; |_")[[1]]})))
                  ind.adist<-unlist(apply(data.frame(clones.tab[sameV[which(CDR3length==unique(CDR3length)[k])],]),1,function(x){rep(x[1],length(strsplit(x[2],split=", |; |_")[[1]]))}))
                  a.dist<-adist(CDR3.adist)#adist(sameCDR3l)  
                  diag(a.dist)<-NA
                  a.dist[upper.tri(a.dist)]<-NA
                  for(a in 1:nrow(a.dist)){
                    clones.tab.new<-vector()
                    grepCDR3<-grep(gsub("[*]","-",do.call(paste, c(as.list(as.character(unique(CDR3.adist[c(a,which(a.dist[a,]<=tr))]))), sep="|"))),gsub("[*]","-",CDR3.adist.temp),perl=T)
                    if(length(which(a.dist[a,]<=tr))>0 && length(unique(ind.adist[grepCDR3]))>1){
                      clones.tab.new<-clones.tab[grep(gsub("[*]","-",do.call(paste, c(as.list(as.character(unique(CDR3.adist[c(a,which(a.dist[a,]<=tr))]))), sep="|"))),gsub("[*]","-",clones.tab$unique_CDR3_sequences_AA),perl=T),]
                      if(length(clones.tab.new)>0 && length(unique(clones.tab.new[,1]))>1){
                        sharedclone.temp<-rbind(sharedclone.temp,c(length(unique(clones.tab.new[,1])),
                                                                   do.call(paste, c(as.list(as.character(unique(clones.tab.new[,1]))), sep="; ")),
                                                                   nchar(CDR3.adist[a]),
                                                                   do.call(paste, c(as.list(as.character(unique(CDR3.adist[c(a,which(a.dist[a,]<=tr))]))), sep="; ")),
                                                                   length(unique(CDR3.adist[c(a,which(a.dist[a,]<=tr))])),
                                                                   do.call(paste, c(as.list(as.character(clones.tab.new$unique_CDR3_sequences_AA)), sep="; ")),
                                                                   do.call(paste, c(as.list(as.character(clones.tab.new$sequence_count_per_CDR3)), sep="; ")),
                                                                   clones.uniV.gene[i],
                                                                   if(useJ==T && length(clones.uniJ.gene[j])>0){do.call(paste, c(as.list(unique(clones.uniJ.gene[j])), sep=", "))}else if(useJ==T && length(clones.uniJ.gene[j])==0){"no J"}else{do.call(paste, c(as.list(unique(clones.tab.new$J_gene)), sep="; "))},
                                                                   
                                                                   if(dispD==T){do.call(paste, c(as.list(unique(clones.tab.new$D_gene)), sep="; "))},
                                                                   if(dispCDR3aa==T){do.call(paste, c(as.list(unique(clones.tab.new$all_CDR3_sequences_AA)), sep="; "))},
                                                                   if(dispCDR3nt==T){do.call(paste, c(as.list(unique(clones.tab.new$all_CDR3_sequences_nt)), sep="; "))},
                                                                   if(dispFunctionality.list==T){do.call(paste, c(as.list(unique(clones.tab.new$Functionality_all_sequences)), sep="; "))},
                                                                   if(dispFunctionality.ratio==T){length(unlist(gregexpr("^productive",clones.tab.new$Functionality_all_sequences))[which(unlist(gregexpr("^productive",clones.tab.new$Functionality_all_sequences))!=(-1))])/length(unlist(apply(data.frame(clones.tab.new$Functionality_all_sequences),1,function(x){strsplit(x,split=", ")[[1]]})))},
                                                                   if(dispFunctionality.ratio==T){length(unlist(gregexpr("^unproductive",clones.tab.new$Functionality_all_sequences))[which(unlist(gregexpr("^unproductive",clones.tab.new$Functionality_all_sequences))!=(-1))])/length(unlist(apply(data.frame(clones.tab.new$Functionality_all_sequences),1,function(x){strsplit(x,split=", ")[[1]]})))},
                                                                   if(dispJunctionFr.list==T){do.call(paste, c(as.list(clones.tab.new$Junction_frame_all_sequences), sep="; "))},
                                                                   if(dispJunctionFr.ratio==T){length(unlist(gregexpr("in-frame",clones.tab.new$Junction_frame_all_sequences))[which(unlist(gregexpr("in-frame",clones.tab.new$Junction_frame_all_sequences))!=(-1))])/length(unlist(apply(data.frame(clones.tab.new$Junction_frame_all_sequences),1,function(x){strsplit(x,split=", ")[[1]]})))},
                                                                   if(dispJunctionFr.ratio==T){length(unlist(gregexpr("out-of-frame",clones.tab.new$Junction_frame_all_sequences))[which(unlist(gregexpr("out-of-frame",clones.tab.new$Junction_frame_all_sequences))!=(-1))])/length(unlist(apply(data.frame(clones.tab.new$Junction_frame_all_sequences),1,function(x){strsplit(x,split=", ")[[1]]})))},
                                                                   if(dispTotalSeq==T){do.call(paste, c(as.list(unique(clones.tab.new$Total_sequences_nt)), sep="; "))}))
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }else{ #############################################useJ==F
        CDR3length<-clones.tab$CDR3_length_AA[sameV]
        if(length(clones.tab[,1][sameV])>1){
          if(length(unique(CDR3length))==1){ # all same CDR3 length
            newCDR3length<-as.numeric(unique(CDR3length))
            tr<-floor(newCDR3length*(1-as.numeric(identity)))
            CDR3.adist.temp<-unlist(apply(data.frame(clones.tab$unique_CDR3_sequences_AA[sameV]),1,function(x){strsplit(x,split=", |; |_")[[1]]}))
            CDR3.adist<-unique(unlist(apply(data.frame(clones.tab$unique_CDR3_sequences_AA[sameV]),1,function(x){strsplit(x,split=", |; |_")[[1]]})))
            ind.adist<-unlist(apply(data.frame(clones.tab[sameV,]),1,function(x){rep(x[1],length(strsplit(x[2],split=", |; |_")[[1]]))}))
            a.dist<-adist(CDR3.adist)#adist(sameCDR3l)  
            diag(a.dist)<-NA
            a.dist[upper.tri(a.dist)]<-NA
            for(a in 1:nrow(a.dist)){
              clones.tab.new<-vector()
              grepCDR3<-grep(gsub("[*]","-",do.call(paste, c(as.list(as.character(unique(CDR3.adist[c(a,which(a.dist[a,]<=tr))]))), sep="|"))),gsub("[*]","-",CDR3.adist.temp),perl=T)
              if(length(which(a.dist[a,]<=tr))>0 && length(unique(ind.adist[grepCDR3]))>1){
                clones.tab.new<-clones.tab[grep(gsub("[*]","-",do.call(paste, c(as.list(as.character(unique(CDR3.adist[c(a,which(a.dist[a,]<=tr))]))), sep="|"))),gsub("[*]","-",clones.tab$unique_CDR3_sequences_AA),perl=T),]
                
                if(length(clones.tab.new)>0 && length(unique(clones.tab.new[,1]))>1){
                  sharedclone.temp<-rbind(sharedclone.temp,c(length(unique(clones.tab.new[,1])),
                                                             do.call(paste, c(as.list(as.character(unique(clones.tab.new[,1]))), sep="; ")),
                                                             nchar(CDR3.adist[a]),
                                                             do.call(paste, c(as.list(as.character(unique(CDR3.adist[c(a,which(a.dist[a,]<=tr))]))), sep="; ")),
                                                             length(unique(CDR3.adist[c(a,which(a.dist[a,]<=tr))])),
                                                             do.call(paste, c(as.list(as.character(clones.tab.new$unique_CDR3_sequences_AA)), sep="; ")),
                                                             do.call(paste, c(as.list(as.character(clones.tab.new$sequence_count_per_CDR3)), sep="; ")),
                                                             clones.uniV.gene[i],
                                                             if(useJ==T && length(clones.uniJ.gene[j])>0){do.call(paste, c(as.list(unique(clones.uniJ.gene[j])), sep=", "))}else if(useJ==T && length(clones.uniJ.gene[j])==0){"no J"}else{do.call(paste, c(as.list(unique(clones.tab.new$J_gene)), sep="; "))},
                                                             
                                                             if(dispD==T){do.call(paste, c(as.list(unique(clones.tab.new$D_gene)), sep="; "))},
                                                             if(dispCDR3aa==T){do.call(paste, c(as.list(unique(clones.tab.new$all_CDR3_sequences_AA)), sep="; "))},
                                                             if(dispCDR3nt==T){do.call(paste, c(as.list(unique(clones.tab.new$all_CDR3_sequences_nt)), sep="; "))},
                                                             if(dispFunctionality.list==T){do.call(paste, c(as.list(unique(clones.tab.new$Functionality_all_sequences)), sep="; "))},
                                                             if(dispFunctionality.ratio==T){length(unlist(gregexpr("^productive",clones.tab.new$Functionality_all_sequences))[which(unlist(gregexpr("^productive",clones.tab.new$Functionality_all_sequences))!=(-1))])/length(unlist(apply(data.frame(clones.tab.new$Functionality_all_sequences),1,function(x){strsplit(x,split=", ")[[1]]})))},
                                                             if(dispFunctionality.ratio==T){length(unlist(gregexpr("^unproductive",clones.tab.new$Functionality_all_sequences))[which(unlist(gregexpr("^unproductive",clones.tab.new$Functionality_all_sequences))!=(-1))])/length(unlist(apply(data.frame(clones.tab.new$Functionality_all_sequences),1,function(x){strsplit(x,split=", ")[[1]]})))},
                                                             if(dispJunctionFr.list==T){do.call(paste, c(as.list(clones.tab.new$Junction_frame_all_sequences), sep="; "))},
                                                             if(dispJunctionFr.ratio==T){length(unlist(gregexpr("in-frame",clones.tab.new$Junction_frame_all_sequences))[which(unlist(gregexpr("in-frame",clones.tab.new$Junction_frame_all_sequences))!=(-1))])/length(unlist(apply(data.frame(clones.tab.new$Junction_frame_all_sequences),1,function(x){strsplit(x,split=", ")[[1]]})))},
                                                             if(dispJunctionFr.ratio==T){length(unlist(gregexpr("out-of-frame",clones.tab.new$Junction_frame_all_sequences))[which(unlist(gregexpr("out-of-frame",clones.tab.new$Junction_frame_all_sequences))!=(-1))])/length(unlist(apply(data.frame(clones.tab.new$Junction_frame_all_sequences),1,function(x){strsplit(x,split=", ")[[1]]})))},
                                                             if(dispTotalSeq==T){do.call(paste, c(as.list(unique(clones.tab.new$Total_sequences_nt)), sep="; "))}))
                }
              }
            }
          }else{ # diff CDR3
            for(k in 1:length(unique(CDR3length))){
              newCDR3length<-as.numeric(unique(CDR3length)[k])
              tr<-floor(newCDR3length*(1-as.numeric(identity)))                                
              sameCDR3l<-clones.tab$unique_CDR3_sequences_AA[sameV][which(CDR3length==unique(CDR3length)[k])]
              if(length(sameCDR3l)>1){
                CDR3.adist.temp<-unlist(apply(data.frame(clones.tab$unique_CDR3_sequences_AA[sameV][which(CDR3length==unique(CDR3length)[k])]),1,function(x){strsplit(x,split=", |; |_")[[1]]}))
                CDR3.adist<-unique(unlist(apply(data.frame(clones.tab$unique_CDR3_sequences_AA[sameV][which(CDR3length==unique(CDR3length)[k])]),1,function(x){strsplit(x,split=", |; |_")[[1]]})))
                ind.adist<-unlist(apply(data.frame(clones.tab[sameV[which(CDR3length==unique(CDR3length)[k])],]),1,function(x){rep(x[1],length(strsplit(x[2],split=", |; |_")[[1]]))}))
                a.dist<-adist(CDR3.adist)#adist(sameCDR3l)  
                diag(a.dist)<-NA
                a.dist[upper.tri(a.dist)]<-NA
                for(a in 1:nrow(a.dist)){
                  clones.tab.new<-vector()
                  grepCDR3<-grep(gsub("[*]","-",do.call(paste, c(as.list(as.character(unique(CDR3.adist[c(a,which(a.dist[a,]<=tr))]))), sep="|"))),gsub("[*]","-",CDR3.adist.temp),perl=T)
                  if(length(which(a.dist[a,]<=tr))>0 && length(unique(ind.adist[grepCDR3]))>1){
                    clones.tab.new<-clones.tab[grep(gsub("[*]","-",do.call(paste, c(as.list(as.character(unique(CDR3.adist[c(a,which(a.dist[a,]<=tr))]))), sep="|"))),gsub("[*]","-",clones.tab$unique_CDR3_sequences_AA),perl=T),]
                    if(length(clones.tab.new)>0 && length(unique(clones.tab.new[,1]))>1){
                      sharedclone.temp<-rbind(sharedclone.temp,c(length(unique(clones.tab.new[,1])),
                                                                 do.call(paste, c(as.list(as.character(unique(clones.tab.new[,1]))), sep="; ")),
                                                                 nchar(CDR3.adist[a]),
                                                                 do.call(paste, c(as.list(as.character(unique(CDR3.adist[c(a,which(a.dist[a,]<=tr))]))), sep="; ")),
                                                                 length(unique(CDR3.adist[c(a,which(a.dist[a,]<=tr))])),
                                                                 do.call(paste, c(as.list(as.character(clones.tab.new$unique_CDR3_sequences_AA)), sep="; ")),
                                                                 do.call(paste, c(as.list(as.character(clones.tab.new$sequence_count_per_CDR3)), sep="; ")),
                                                                 clones.uniV.gene[i],
                                                                 if(useJ==T && length(clones.uniJ.gene[j])>0){do.call(paste, c(as.list(unique(clones.uniJ.gene[j])), sep=", "))}else if(useJ==T && length(clones.uniJ.gene[j])==0){"no J"}else{do.call(paste, c(as.list(unique(clones.tab.new$J_gene)), sep="; "))},
                                                                 
                                                                 if(dispD==T){do.call(paste, c(as.list(unique(clones.tab.new$D_gene)), sep="; "))},
                                                                 if(dispCDR3aa==T){do.call(paste, c(as.list(unique(clones.tab.new$all_CDR3_sequences_AA)), sep="; "))},
                                                                 if(dispCDR3nt==T){do.call(paste, c(as.list(unique(clones.tab.new$all_CDR3_sequences_nt)), sep="; "))},
                                                                 if(dispFunctionality.list==T){do.call(paste, c(as.list(unique(clones.tab.new$Functionality_all_sequences)), sep="; "))},
                                                                 if(dispFunctionality.ratio==T){length(unlist(gregexpr("^productive",clones.tab.new$Functionality_all_sequences))[which(unlist(gregexpr("^productive",clones.tab.new$Functionality_all_sequences))!=(-1))])/length(unlist(apply(data.frame(clones.tab.new$Functionality_all_sequences),1,function(x){strsplit(x,split=", ")[[1]]})))},
                                                                 if(dispFunctionality.ratio==T){length(unlist(gregexpr("^unproductive",clones.tab.new$Functionality_all_sequences))[which(unlist(gregexpr("^unproductive",clones.tab.new$Functionality_all_sequences))!=(-1))])/length(unlist(apply(data.frame(clones.tab.new$Functionality_all_sequences),1,function(x){strsplit(x,split=", ")[[1]]})))},
                                                                 if(dispJunctionFr.list==T){do.call(paste, c(as.list(clones.tab.new$Junction_frame_all_sequences), sep="; "))},
                                                                 if(dispJunctionFr.ratio==T){length(unlist(gregexpr("in-frame",clones.tab.new$Junction_frame_all_sequences))[which(unlist(gregexpr("in-frame",clones.tab.new$Junction_frame_all_sequences))!=(-1))])/length(unlist(apply(data.frame(clones.tab.new$Junction_frame_all_sequences),1,function(x){strsplit(x,split=", ")[[1]]})))},
                                                                 if(dispJunctionFr.ratio==T){length(unlist(gregexpr("out-of-frame",clones.tab.new$Junction_frame_all_sequences))[which(unlist(gregexpr("out-of-frame",clones.tab.new$Junction_frame_all_sequences))!=(-1))])/length(unlist(apply(data.frame(clones.tab.new$Junction_frame_all_sequences),1,function(x){strsplit(x,split=", ")[[1]]})))},
                                                                 if(dispTotalSeq==T){do.call(paste, c(as.list(unique(clones.tab.new$Total_sequences_nt)), sep="; "))}))
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    return(data.frame(sharedclone.temp,stringsAsFactors = F))
  }
  
  sharedclone<-do.call(rbind.data.frame, clonelist)
  
  if(length(sharedclone)>0){
    if(length(which(as.numeric(sharedclone[,1])==1))>0){
      sharedclone<-sharedclone[which(as.numeric(sharedclone[,1])>1),]
    }
    if(nrow(sharedclone)>0){
      colnames(sharedclone)<-c("number_individuals",
                               "individuals",
                               "CDR3_length_AA",
                               "shared_CDR3",
                               "number_shared_CDR3",
                               "CDR3_sequences_per_individual",
                               "sequence_count_per_CDR3",
                               "V_gene",
                               "J_gene", 
                               if(dispD==T){"D_gene"},
                               if(dispCDR3aa==T){"all_CDR3_sequences_AA"},
                               if(dispCDR3nt==T){"all_CDR3_sequences_nt"},
                               if(dispFunctionality.list==T){"Functionality_all_sequences"},
                               if(dispFunctionality.ratio==T){c("Func_productive_sequences","Func_unproductive_sequences")},
                               if(dispJunctionFr.list==T){"Junction_frame_all_sequences"}, 
                               if(dispJunctionFr.ratio==T){c("JF_in_frame","JF_out_of_frame")},
                               if(dispTotalSeq==T){"Total_sequences_nt"})    
      
      sharedclone<-sharedclone[intersect(which(duplicated(sharedclone$individuals)==F),which(duplicated(sharedclone$"CDR3_sequences_per_individual")==F)),]
      
      if(nrow(sharedclone)>1){
        sharedclone<-sharedclone[order(sharedclone$number_shared_CDR3),]
        out<-vector()
        for(i in 1:nrow(sharedclone)){
          if(length(grep(gsub("[*]","-",sharedclone[i,"shared_CDR3"]),gsub("[*]","-",sharedclone[(i+1):nrow(sharedclone),"shared_CDR3"]),perl=T))>1 && length(unique(sharedclone$individuals[grep(gsub("[*]","-",sharedclone[i,"shared_CDR3"]),gsub("[*]","-",sharedclone[(i+1):nrow(sharedclone),"shared_CDR3"]),perl=T)]))==1){
            out<-c(out,i)
          }
        }
        out<-unique(c(out,grep("NA",sharedclone[,"individuals"])))
        if(length(out)>0){
          sharedclone<-sharedclone[-out,]
        }
      }
    }else{
      sharedclone<-vector()
    }
  }
  return(sharedclone)
}


