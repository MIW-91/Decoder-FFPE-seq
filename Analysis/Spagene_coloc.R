library(Seurat)
library(readxl)
library(ggplot2)
library(SpaGene)
library(dplyr)
library(tidyverse)
library(RColorBrewer)

def_plotLR<-function(expr,location,normalize=TRUE,topn=floor(0.2*dim(location)[1]),knn=8,LRpair=c("Ptn","Ptprz1"),pt.size=2,alpha.min=0.1,max.cut=0.95){
  if (sum(rownames(expr) %in% LRpair)!=2) { stop("ligand or receptor are not expressed")}
  nnmatrix<-RANN::nn2(location,k=knn)$nn.idx
  countsum<-Matrix::colSums(expr)
  
  ncell<-dim(expr)[2]
  if (normalize==TRUE) {
    expr<-Matrix::t(log(Matrix::t(expr)/countsum*median(countsum)+1))
  }
  
  
  ligand<-expr[LRpair[1],]
  receptor<-expr[LRpair[2],]
  LRexp<-rbind(ligand,receptor)
  neighexp<-apply(nnmatrix,1,function(x){apply(LRexp[,x[2:knn]],1,max)})
  
  #LRexp<-t(scale(t(LRexp)))
  #neighexp<-t(scale(t(neighexp)))
  #LRexp[LRexp<0]<-0
  #neighexp[neighexp<0]<-0
  LRadd<-pmax(LRexp[1,]*neighexp[2,],LRexp[2,]*neighexp[1,])
  LRadd_max<-quantile(LRadd,probs=max.cut)
  LRadd[LRadd>LRadd_max]<-LRadd_max
  if (sum(ligand>0)>topn) {n1<-order(ligand,sample(ncell,ncell),decreasing=T)[1:topn]} else{n1<-which(ligand>0)}
  if (sum(receptor>0)>topn) {n2<-order(receptor,sample(ncell,ncell),decreasing=T)[1:topn]} else{n2<-which(receptor>0)}
  expcol<-rep(0,ncell)
  expcol[n1]<-1
  expcol[n2]<-2
  expcol[intersect(n1,n2)]<-3
  tmp<-data.frame(x=location[,1],y=location[,2],Exp=as.factor(expcol))
  tmpLRadd<-data.frame(x=location[,1],y=location[,2],LR=LRadd)
  
  alpha=(LRadd-min(LRadd))/(max(LRadd)-min(LRadd))*(1-alpha.min)+alpha.min
  
  p1<-ggplot(tmp,aes(x=x,y=y,col=Exp))+
    geom_point(size=pt.size)+
    scale_color_manual(values=c("#e8eff4ff","#3e8a86","#5392f5","#b81422"),labels=c("Both low","lncRNA high","mRNA High","Both High"))+
    ggtitle(paste0(LRpair,collapse="_"))+
    xlab("")+ylab("")+
    theme_bw()+
    theme(axis.line=element_blank(),
          axis.text.x=element_blank(), axis.text.y=element_blank(),
          axis.ticks.x=element_blank(),axis.ticks.y=element_blank())
  p2<-ggplot(tmpLRadd,aes(x=x,y=y,col=LR))+
    geom_point(size=pt.size,alpha=alpha)+
    scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(9,"RdYlBu")))(100))+
    xlab("")+ylab("")+
    theme_bw()+
    theme(axis.line=element_blank(),
          axis.text.x=element_blank(), axis.text.y=element_blank(),
          axis.ticks.x=element_blank(),axis.ticks.y=element_blank())+
    labs(color = "colocalization level")
  p1+p2&scale_y_reverse()
}



mRNA = c('TLR4', 'CCL2', 'NLRP3', 'SFPQ', 'CXCL8','IL1B', 'NFKBIA', 'MMP1','CSF3','G0S2','CXCL3')
lnc = c('NEAT1','MALAT1','SFPQ', 'CXCL8')
pair_matrix <- expand.grid(lnc = lnc, mRNA = mRNA)
rownames(pair_matrix) <- paste(pair_matrix$lnc, pair_matrix$mRNA, sep = " - ")

expr<-data[['SCT']]$data
geometry<-GetTissueCoordinates(data,cols = c("imagerow", "imagecol"), scale = NULL)
geometry$cell<-NULL
geometry<-geometry[match(colnames(expr),rownames(geometry)),]

OSCC_lr<-SpaGene_LR(expr,geometry,normalize = F,LRpair=pair_matrix) 
OSCC_lr<-filter(OSCC_lr,adjp<0.05)
write.csv(OSCC_lr,'sig.lnc_mRNA.SCT.csv',row.names=T)

sig.coloc<-rownames(OSCC_lr)
lapply(sig.coloc,function(x){
  str<-strsplit(x, " - ")[[1]]
  lncg=str[1]
  mg=str[2]
  # p=def_plotLR(count,geometry,LRpair=c(lncg,mg),alpha.min=0.1)
  p=def_plotLR(expr,geometry,LRpair=c(lncg,mg),alpha.min=0.1,normalize = F)
  pdf(paste0(x,'.coloc.SCT.pdf'),12,6)
  print(p)
  dev.off()
})