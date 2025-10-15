####obtain a pseudo bulk expression data from a SRT data by k-means clustering based on the coordinate matrix
defcol1<-c('#1C79B7','#F38329','#2DA248','#DC403E','#976BA6','#8D574C','#D07DB0','#BDBF3B',
           '#27BDD0','#B0C7E5','#F9BA79','#A0CC8A','#F09EA1','#C7B0D1','#D4BCB8','#F7CCDD',
           '#DBDD8D','#A3D7E5','#B04E4F','#A38A58','#ED5351','#0C8945','#0F71A7','#A82764',
           '#F8DAE4','#7A4B1F','#5E6BAE','#8AC997','#DAC9AC','#0F5148','#A0B5DC','#9F858E',
           '#5C181D','#7B8380','#E8DDDA','#264220','#5AB747','#5169AE','#4B3D58','#CD428A','#62615A','#B82129','#66762E')

# data: rds file path of a SRT data
# n: number of spots in a cluster
# out.table: path of the output pseudo bulk data table
# out.figure: path of the output boxplot

Pseudobulk_expr<-function(data,n = 50,out.table,out.figure){
  srt<-readRDS(data)
  coords <- GetTissueCoordinates(srt,cols = c("imagerow", "imagecol"), scale = NULL)  
  set.seed(123)
  k <- ceiling(nrow(coords) / n)  
  clusters <- kmeans(coords[, 1:2], centers = k)$cluster
  coords$group <- paste0("cluster_", clusters)  
  coords<-coords[match(rownames(srt@meta.data),rownames(coords)),]
  srt$pseudo_group <- coords$group
  pseudo_bulk <- AggregateExpression( 
    srt,
    features = gene,
    assays = "SCT",              
    slot = "data",             
    group.by = "pseudo_group",    
    return.seurat = FALSE,       
    fun = "sum"                   
  )
  pseudo_matrix <- pseudo_bulk$SCT  
  pseudo_matrix<-data.frame(t(pseudo_matrix),check.names = F)     
  pseudo_matrix<-pivot_longer(pseudo_matrix,cols = colnames(pseudo_matrix),names_to = 'gene',values_to = 'expr')
  write.csv(pseudo_matrix,out.table)
  p=ggplot(pseudo_matrix,aes(x=gene,y=expr,fill=gene,col=gene))+
    geom_jitter(aes(col=gene),size=0.5,width = 0.1)+
    geom_boxplot(aes(col=gene),size=0.5,width=0.5,alpha=0,outlier.shape = NA)+
    scale_fill_manual(values = defcol1[10:20])+         
    scale_color_manual(values = defcol1[10:20])+
    xlab("")+                                                   
    ylab("Expression level")+    
    theme_bw()+
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.border=element_rect(size=1.2),                  
          axis.text.x = element_text(angle=0,size=10,vjust = 1,hjust =1,color = "black"),                                                   
          axis.text.y = element_text(size =10),
          legend.position = 'right' )
  
  pdf(out.figure,10,4)
  print(p)
  dev.off()
  
})