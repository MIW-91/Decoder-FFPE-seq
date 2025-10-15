commot<-lapply(c('8_um-post-pNR-1','8_um-pre-pNR-1','8_um-post-MPR-1','8_um-pre-MPR-1'),function(x){
  fname=paste0('commot_filepath/',x,'_lrpair_commot.csv')
  tmp<-read.csv(fname,header = T)
  tmp<-tmp %>%
    filter(celltype_target %in% c('CD8_Tn','CD8_EarlyActiv','CD8_Tem','CD8_Tex') & p.values<0.001) %>%
    mutate(sample=x)
  return(tmp)
}) %>% bind_rows()


commot$class<-ifelse(commot$sample=='8_um-post-MPR-1','responder_post',
                     ifelse(commot$sample=='8_um-pre-MPR-1','responder_pre',
                            ifelse(commot$sample=='8_um-post-pNR-1','non-responder_post','non-responder_pre')))

print(table(commot$class))


calculate_ratio<-function(data,condition){
  data_c<-filter(data,class==condition)
  lr=unique(data_c$LR.pair)
  res<-lapply(lr,function(x){
    c_ik=sum(data_c$Communication_score[data_c$LR.pair==x])
    c_i=sum(data$Communication_score[data$LR.pair==x])
    C_k=sum(data_c$Communication_score)
    C_m=sum(data$Communication_score)
    ratio=(c_ik/c_i)/(C_k/C_m)
    return(ratio)
  })
  names(res)<-lr
  res<-unlist(res)
  return(res)
}

calculate_countratio<-function(data,condition){
  data_c<-filter(data,class==condition)
  lr=unique(data_c$LR.pair)
  res<-lapply(lr,function(x){
    c_ik=sum(data_c$LR.pair==x)
    c_i=sum(data$LR.pair==x)
    C_k=nrow(data_c)
    C_m=nrow(data)
    ratio=(c_ik/c_i)/(C_k/C_m)
    return(ratio)
  })
  names(res)<-lr
  res<-unlist(res)
  return(res)
}

comp<-list(
  c('responder_post','responder_pre'),
  c('non-responder_post','non-responder_pre'),
  c('non-responder_pre','responder_pre'),
  c('non-responder_post','responder_post'),
  rev(c('responder_post','responder_pre')),
  rev(c('non-responder_post','non-responder_pre')),
  rev(c('non-responder_pre','responder_pre')),
  rev(c('non-responder_post','responder_post'))
)
class<-unlist(lapply(comp,function(k) paste0(k[1],'.',k[2])))
print(class)


Ratio<-sapply(unique(commot$class),function(k){
  result<-calculate_ratio(data=commot,condition = k)
  return(result)
}) %>% bind_rows() %>%
  t() %>%
  data.frame()
colnames(Ratio)<-unique(commot$class)
write.csv(Ratio,'enrich_interaction_by_strength.csv',row.names = T)

Ratio_c<-sapply(unique(commot$class),function(k){
  result<-calculate_countratio(data=commot,condition = k)
  return(result)
}) %>% bind_rows() %>%
  t() %>%
  data.frame()
colnames(Ratio_c)<-unique(commot$class)
write.csv(Ratio_c,'enrich_interaction_by_frequency.csv',row.names = T)

library(dplyr)
library(tibble)


Ratio_long <- Ratio %>%
  rownames_to_column("LRpair") %>%
  pivot_longer(cols = -LRpair, names_to = "class", values_to = "value")

# 分组计算前10%
results<- Ratio_long %>%
  group_by(class) %>%
  filter(!is.na(value)) %>%
  arrange(class, desc(value)) %>%
  slice_head(prop = 0.25) %>%  
  ungroup()

defcol1<-c('#1C79B7','#F38329','#2DA248','#DC403E','#976BA6','#8D574C','#D07DB0','#BDBF3B',
           '#27BDD0','#B0C7E5','#F9BA79','#A0CC8A','#F09EA1','#C7B0D1','#D4BCB8','#F7CCDD',
           '#DBDD8D','#A3D7E5','#B04E4F','#A38A58','#ED5351','#0C8945','#0F71A7','#A82764',
           '#F8DAE4','#7A4B1F','#5E6BAE','#8AC997','#DAC9AC','#0F5148','#A0B5DC','#9F858E',
           '#5C181D','#7B8380','#E8DDDA','#264220','#5AB747','#5169AE','#4B3D58','#CD428A','#62615A','#B82129','#66762E')

Tc= c('CD8_Tn','CD8_EarlyActiv','CD8_Tem','CD8_Tex')
n=sum(is.na(Ratio$responder_pre))
signal=c('SELPLG-SELL','CCL5-CCR4','CCL5-CCR5','CCL19-CCR7','MIF-CD74_CXCR4','SEMA4D-PLXNB2','GZMA-F2R','CD99-CD99L2')
lrp1<-read.csv('8_um-post-MPR-1_lrpair_commot.csv',header=T)

lr1<-filter(lrp1,celltype_target %in% c('CD8_EarlyActiv') & celltype_source %in% c( 'B_cells','Th1','CD4_Tn', 'CD8_Tn', 'CD8_Tem') & 
             LR.pair %in% signal
           # p.values<0.05
           )
a=unique(lr1$pathway)
lr1$ccc<-paste0(lr1$celltype_source,' -> ', lr1$celltype_target)
lr1$logp<--log10(lr1$p.values+1e-6)
p=ggplot(lr1, aes(y = LR.pair, x = ccc, color=pathway)) +
  geom_point(aes(size = Communication_score, alpha = logp)) +
  scale_color_manual(values = defcol) +
  theme_bw()+
  scale_x_discrete(position = "top")+
  theme(axis.text.x = element_text(angle = 90))

pdf('~/DecoderFFPE/result/commot/8_um-post-MPR-1.CD8Ta.new.pdf',4.5,4.5)  
print(p)
dev.off()


signal=c('FN1-ITGA5_ITGB1','THBS2-CD47','FN1-CD44','CSF1-CSF1R','CD99-CD99','CD46-JAG1','COL4A1-SDC1','COL1A1-CD44')
lrp2<-read.csv('~/DecoderFFPE/result/commot/8_um-post-pNR-1_lrpair_commot.csv',header=T)
lr<-filter(lrp2,celltype_target %in% 'CD8_Tex' & celltype_source %in%  c( 'CAF','CD4_Tn', 'CD8_Tem','Monocytes','DC','Mast_cells','Tfh','Th1','Plasma')&
              LR.pair %in% signal
             # p.values<0.05
             )
lr$ccc<-paste0(lr$celltype_source,' -> ', lr$celltype_target)
lr$logp<--log10(lr$p.values+1e-6)
b=unique(lr$pathway)
defcol<-defcol1[1:length(unique(c(a,b)))]
names(defcol)<-unique(c(a,b))

p=ggplot(lr, aes(y = LR.pair, x = ccc, color=pathway)) +
  geom_point(aes(size = Communication_score, alpha = logp)) +
  scale_color_manual(values = defcol) +
  theme_bw()+
  scale_x_discrete(position = "top")+
  theme(axis.text.x = element_text(angle = 90))

pdf('~/DecoderFFPE/result/commot/8_um-post-pNR-1.CD8Tex.new.pdf',5,4)  
print(p)
dev.off()

