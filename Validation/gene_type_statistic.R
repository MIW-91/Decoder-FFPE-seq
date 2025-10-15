library(rtracklayer);library(dplyr);library(Seuart)

# gtf <- import('/home/disk/miw/GRCm38_bundle/refdata-gex-mm10-2020-A/genes/genes.gtf', format="gtf")
gtf.ori <- import('E:\\decoderFFPE\\Mus_musculus.GRCm38.102.gtf', format="gtf")
gtf.ori<-as.data.frame(gtf.ori)
data<-Load10X_Spatial('E:/decoderFFPE/seqData/RJ-GL-241221-3')


# dim(gtf)
# unique(gtf$gene_type)
gtf<-gtf.ori[,c('gene_id','gene_name','gene_biotype')]

gtf$gene_type<-gtf$gene_biotype
gtf$gene_type[gtf$gene_type%in%c('protein_coding','TR_V_gene','TR_D_gene','TR_J_gene','TR_C_gene','IG_LV_gene','IG_V_gene','IG_J_gene','IG_C_gene','IG_D_gene')]<-'protein_coding'
gtf$gene_type[gtf$gene_type%in%c('lincRNA','3prime_overlapping_ncRNA','antisense','macro_lncRNA','sense_intronic','sense_overlapping')]<-'lncRNA'
gtf$gene_type[gtf$gene_type%in%c('miRNA')]<-'miRNA'
gtf$gene_type[gtf$gene_type%in%c('rRNA')]<-'rRNA'
gtf$gene_type[gtf$gene_type%in%c('processed_pseudogene','transcribed_unprocessed_pseudogene','unprocessed_pseudogene','transcribed_processed_pseudogene','unitary_pseudogene','pseudogene','polymorphic_pseudogene','transcribed_unitary_pseudogene','translated_unprocessed_pseudogene','TR_V_pseudogene','TR_J_pseudogene','IG_V_pseudogene','IG_C_pseudogene','IG_D_pseudogene','IG_pseudogene')]<-'pseudogene'
gtf$gene_type<-ifelse(gtf$gene_type%in%c('pseudogene','rRNA','miRNA','lncRNA','protein_coding'),gtf$gene_type,'other_types')
table(gtf$gene_type)

gtf$gene_name<-ifelse(is.na(gtf$gene_name),gtf$gene_id,gtf$gene_name)
gtf<-gtf[!duplicated(gtf$gene_id),]


gene_df<-data.frame(gene_name=rownames(tdata))

gene_df<-left_join(gene_df,gtf,by='gene_name')

gene_df %>% distinct(gene_id,gene_name,.keep_all = TRUE) -> gene_df

if(sum(is.na(gene_df$gene_name))==0){
  gene_df %>% distinct(gene_name,.keep_all = TRUE) -> gene_df
}else{
  gene_df$gene_name<-gene_df$gene_id
}

gene_df$gene_type[which(is.na(gene_df$gene_type))]<-'other_types'

stat.d<-as.data.frame(table(gene_df$gene_type))

stat.d$Prop<-stat.d$Freq/sum(stat.d$Freq)
stat.d$Prop.r<-round(stat.d$Prop,4)
colnames(stat.d)<-c('gene_type','frequency','prop_prec','proportion')

usecol=c('#DBDD8D',"#E68FAC","#A1CAF1",'#F9BA79','#976BA6',"#FBE426")
names(usecol)=c('protein_coding','lncRNA','miRNA','rRNA','pseudogene','other_types')

ggplot(stat.d,aes(x=2,y=proportion,fill=gene_type))+
  geom_bar(stat='identity')+
  coord_polar("y", start = 0)+
  # geom_text(aes(y = pos, label = paste0(prop*100,'%')), color = "black")+
  scale_fill_manual(values = usecol,breaks = stat.d$gene_type,labels=paste0(stat.d$gene_type,': ',stat.d$proportion*100,'%')) +
  theme_void()+
  theme(legend.position='right')+
  xlim(0.5, 2.5)
