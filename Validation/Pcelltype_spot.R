library(ggplot2)
acell<-read.csv('st_cell2location_res.csv',header = T,check.names = F,row.names = 1)
colnames(acell)<-gsub('q05cell_abundance_w_sf_','',colnames(acell))



x=0.1
# x=0.5

stat<-apply(acell,1,function(y){
  sum(y>=x)
})
table(stat>0)

stat<-ifelse(stat<1,1,stat)

stat1<-data.frame(num=c('1','2','3','≥4'),freq=c(sum(stat==1),sum(stat==2),sum(stat==3),sum(stat>=4)))
stat1$perc<-stat1$freq/length(stat)
sum(stat1$perc)

p1 <- sum(rowSums(acell >= x) == 1) / sum(rowSums(acell >= x) >= 1)
# 计算概率 p2
p2 <- sum(rowSums(acell >= x) == 2) / sum(rowSums(acell >= x) >= 1)
# 计算概率 p3
p3 <- sum(rowSums(acell >= x) == 3) / sum(rowSums(acell >= x) >= 1)
p4 <- sum(rowSums(acell >= x) >= 4) / sum(rowSums(acell >= x) >= 1)

stat1<-data.frame(num=c('1','2','3','≥4'),perc=c(p1,p2,p3,p4))
sum(stat1$perc)

stat1$num<-factor(stat1$num,levels=c('1','2','3','≥4'))

ggplot(stat1,aes(x=num,y=round(perc*100,0)))+
  geom_bar(stat="identity",alpha=0.8,fill = "#A1CAF1")+
  geom_label(  aes(label = paste0(round(perc*100,2),'%')), data=stat1,  nudge_y = 1,  alpha = 0.5)+
  theme_bw() +
  scale_y_continuous(limits = c(0,100))+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())+
  labs( x="numbers of cell types per spot", y="Percentage")+
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size = 16))


