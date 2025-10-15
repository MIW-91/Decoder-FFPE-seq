library(dplyr)
library(Seurat)
library(ggplot2)
library(biomaRt)
library(data.table)
library(readxl)
##########mouse GPCR gene
mus38g<-read_xlsx('MouseHumanRatGPCRs.xlsx')
colnames(mus38g)<-mus38g[2,]
mus38g<-mus38g[-c(1:2),]
mus38g<-filter(mus38g,!Subgroup %in% c('Opsin receptors','Taste 1 receptors','Taste 2 receptors','Vomeronasal receptors')) %>%
    filter(!grepl('Mrgpr|Taar', `Gene Symbol (Mouse)`))

##########load and subset MERFISH
merfish<-read.csv('../datasets_mouse_brain_map_BrainReceptorShowcase_Slice1_Replicate1_cell_by_gene_S2R1.csv',header = T,row.names = 1,check.names = F)
position<-read.csv('datasets_mouse_brain_map_BrainReceptorShowcase_Slice2_Replicate1_cell_metadata_S2R1.csv',header = T,row.names = 1,check.names = F)
merfish<-CreateSeuratObject(counts = t(merfish),assay = 'Spatial')
position<-position[,c('center_x','center_y')]
colnames(position) = paste0("Spatial_",1:ncol(position))
merfish[["spatial"]] <- CreateDimReducObject(embeddings = as.matrix(position), 
                                          key = "Spatial",assay = "Spatial")
coord<-data.frame(merfish@reductions$spatial@cell.embeddings)
coord<-coord[coord$Spatial_2>=3500 & coord$Spatial_1>=4200 & coord$Spatial_1<=7000,]
merfish<-subset(merfish,cells = rownames(coord))


mfgene<-Features(mf)[!grepl('Blank',Features(merfish))]
merfish<-subset(merfish,features = mfgene)

mfexpr<-GetAssayData(merfish,assay = 'Spatial',layer = 'counts')
mfexpr<-data.frame(mfexpr)

mfexpr<-data.frame(X=rownames(mfexpr),total_expression=rowSums(mfexpr),mean_expression=rowMeans(mfexpr))
mfexpr['ADGRF3','X']='Adgrf3'
mfexpr<-filter(mfexpr,X %in% mus38g$`Gene Symbol (Mouse)`)

##########load SRT data
platform<-'Decoder-FFPE-Seq'
data<-Load10X_Spatial('Decoder-FFPE-seq_8um_brain')
STexpr<-GetAssayData(data,assay = 'Spatial',layer = 'counts')

expr<-data.frame(total=rowSums(STexpr),mean=rowMeans(STexpr))
expr$X=rownames(expr)
expr<-filter(expr,X %in% c(mus38g$`Gene Symbol (Mouse)`))
write.csv(expr,paste0(platform,'_nonsensory_GPCR.csv'),row.names = F)

##########calculate correlation             
expr<-left_join(mfexpr,expr,by = 'X')
length(which(is.na(expr$total)))
expr[which(is.na(expr$total)),'total']=0
expr[which(is.na(expr$mean)),'mean']=0
expr[which(is.na(expr$total_expression)),'total_expression']=0
expr[which(is.na(expr$mean_expression)),'mean_expression']=0

expr.st<-expr$mean
expr.sc<-expr$mean_expression
res<-cor.test(expr.st,expr.sc,method = 'spearman')
corr <-round(res$estimate,4)
pval <-max(res$p.value,2.2e-16)
fplot<-data.frame(ST=expr.st,sc=expr.sc)

ggplot(fplot,aes(x=log(ST+1e-6),y=log(sc+1e-6)))+
  geom_point(alpha=1,size=1) +
  scale_x_continuous(n.breaks = 7)+
  scale_y_continuous(n.breaks = 6)+
  geom_smooth(method = lm,linetype=1,color='red',se=F)+
  guides(alpha='none')+
  labs(x=paste0('Log(mean counts)\n',platform),
       y='Log(mean counts)\nMERFISH')+
  theme_bw()+
  theme(legend.position = 'right',
        axis.text.x = element_text(hjust = 1,colour = 'black',size=14),
        axis.text.y = element_text(colour = 'black',size=14,face = 'plain'),
        axis.title.x = element_text(colour = 'black',size=18,face = 'plain'),
        axis.title.y = element_text(colour = 'black',size=18,face = 'plain'),
        panel.border = element_blank(),
        axis.line = element_line(colour = 'black',size=1),
        legend.text = element_text(colour = 'black',size=14,face = 'plain'),
        legend.title = element_text(colour = 'black',size=18,face = 'plain'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  annotate("text",x=-10,y=5,label=paste0('R = ',corr,"\np < ",pval),color='black',size=6,hjust=0)

##########compare with other platforms
expr<-read.csv('Decoder-FFPE-seq_nonsensory_GPCR.csv',header = T)
colnames(expr)<-c('Decoder-FFPE-seq_total','Decoder-FFPE-seq_mean','gene')
for(i in c('Visium HD','Stereo-N-FFPE')){
  tmp<-read.csv(paste0(i,'_nonsensory_GPCR.csv'),header = T)
  colnames(tmp)<-c(paste0(i,'_total'),paste0(i,'_mean'),'gene')
   expr<-full_join(expr,tmp,by='gene')
}

library(tidyr)

fplot<-na.omit(expr)

expr[is.na(expr)] =0
fplot <- expr 

rownames(fplot)<-fplot$gene
fplot$gene<-NULL

##non_zero
fplot <- fplot %>%
  filter(across(all_of(colnames(fplot)), ~ .x > 0))

fplot<-data.frame(t(colSums(fplot)),check.names = F)


fplot <- fplot %>%
  pivot_longer(
    cols = colnames(fplot),  
    names_to = c("technology", "statistic"),  
    names_sep = "_",  
    values_to = "Total UMI" 
  )
write.csv(fplot,'nonsensory_GPCRs_total_UMI.csv',row.names = F)


fplot$UMI_um<-ifelse(fplot$technology=='Stereo-N-FFPE',fplot$`Total UMI`/(10^2), fplot$`Total UMI`/(8^2))
fplot$technology<-factor(fplot$technology,levels = c("Decoder-FFPE-seq","Visium HD", "Stereo-N-FFPE"))

fplot2<-filter(fplot,statistic == 'mean')
p=ggplot(fplot2, aes(x = technology, y = `UMI_um`,fill = technology)) +
  geom_bar(stat = "identity",  width = 0.9) +
  labs(
    title = "Sum of average expression of GPCR genes in MERFISH panel",                      
    x = "Platform",                        
    y = "Total UMIs per um"                           
  #scale_y_log10(n.breaks = 6)+
  scale_fill_manual(values = c(#'MERFISH'="#FED9A6",
    "Decoder-FFPE-seq"="#FBB4AE","Visium HD" = "#B3CDE3", "Stereo-N-FFPE"= "#DECBE4"))+
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), 
    panel.grid.major = element_blank(),                
    panel.grid.minor = element_blank(),                                
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45,hjust=1)
  ) 
p
pdf('Sum of average expr_GPCR_um.pdf',7,5)
print(p)
dev.off()
fplot2$`Total UMI`<-NULL
write.csv(fplot2,'nonsensory_GPCRs.4polt.csv',row.names = F)


fplot1<-filter(fplot,statistic == 'total')
p=ggplot(fplot1, aes(x = technology, y = `Total UMI`,fill = technology)) +
  geom_bar(stat = "identity",  width = 0.9) +
  labs(
    title = "Total UMIs of GPCR genes in MERFISH panel",                      
    x = "Platform",                     
    y = "Total UMIs"                           
  ) +
  scale_y_log10(n.breaks=6)+
  scale_fill_manual(values = c('MERFISH'="#FED9A6","Decoder-FFPE-seq"="#FBB4AE","Visium HD" = "#B3CDE3", "Stereo-N-FFPE"= "#DECBE4"))+
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), 
    panel.grid.major = element_blank(),            
    panel.grid.minor = element_blank(),                               
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45,hjust=1)
  ) 
p
pdf('Total UMI_359g.pdf',7,5)
print(p)
dev.off()
