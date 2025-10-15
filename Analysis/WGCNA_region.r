library(WGCNA);library(hdWGCNA);library(tidyverse);library(cowplot);library(patchwork)
library(Seurat);library(dplyr);library(ggrepel);library(ggplot2);library(RColorBrewer)
defcol1<- c("#F3C300", "#875692", "#F38400", "#A1CAF1", "#BE0032", "#C2B280",
           "#008856", "#E68FAC", "#0067A5", "#F99379", "#604E97", "#F6A600", "#B3446C", 
           "#DCD300", "#882D17", "#8DB600", "#654522", "#E25822", "#2B3D26", "#5A5156", "#E4E1E3", 
           "#F6222E", "#FE00FA", "#16FF32", "#3283FE", "#FEAF16", "#B00068", "#1CFFCE", "#90AD1C",
           "#2ED9FF", "#DEA0FD", "#AA0DFE", "#F8A19F", "#325A9B", "#C4451C", "#1C8356", "#85660D",
           "#B10DA1", "#FBE426", "#1CBE4F", "#FA0087", "#FC1CBF", "#F7E1A0", "#C075A6", "#782AB6")
sample='my_sample'
fpath='my_direction'
# dir.create(fpath)
setwd(fpath)
# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)

# optionally enable multithreading
enableWGCNAThreads(nThreads = 8)

data<-readRDS(paste0(sample,'.rds'))
levels(data)
data$subcluster<-factor(data$seurat_clusters,levels = c(0,1,2,3,4,5,6,7,8,9,10)) # we have 10 clusters
Idents(data)<-data$subcluster

DefaultAssay(data)<-'Spatial'


library(rtracklayer)
gtf <- import('/Homo_sapiens.GRCh38.112.gtf', format="gtf")
gtf<-as.data.frame(gtf)

colnames(gtf)
gtf<-gtf[gtf$gene_biotype=='lncRNA',c('gene_id','gene_name','gene_biotype')]
gtf<-gtf[!duplicated(gtf$gene_id),]
head(gtf)
lnc<-na.omit(unique(gtf$gene_name,gtf$gene_id))
lnc<-intersect(lnc,rownames(data))

ppath='/data1/zhaoshutao/miw/decoder-FFPE/seqdata/lncRNA0708/'
ppath=paste0(ppath,sample,'/spatial/')
# add barcode column to Seurat obj 
data$barcode <- colnames(data)

# update this with the path for your sample.
tissue_positions <- read.csv(paste0(ppath,'tissue_positions_list.csv'),
                    header=F)
colnames(tissue_positions)=c('barcode','in_tissue','array_row','array_col','pxl_row_in_fullres','pxl_col_in_fullres')

tissue_positions <- subset(tissue_positions, barcode %in% data$barcode)
data
dim(tissue_positions)
# join the image_df with the Seurat metadata
new_meta <- dplyr::left_join(data@meta.data, tissue_positions, by='barcode')

# add the new metadata to the seurat seurat_vhd
data$row <- new_meta$array_row 
data$imagerow <- new_meta$pxl_row_in_fullres
data$col <- new_meta$array_col
data$imagecol <- new_meta$pxl_col_in_fullres

wgcna<- SetupForWGCNA(data, gene_select = "features", # the gene selection approach
                          features = Features(data), # fraction of cells that a gene needs to be expressed in order to be included
                          wgcna_name = "Tumor" # the name of the hdWGCNA experiment,we mainly detected network of tumor region
)

# construct metacells  in each group
seu_wgcna <- MetaspotsByGroups(
  wgcna,
  group.by = c("region"),
  ident.group = "region", # region contains Tumor and Stroma
  assay = 'Spatial',
  slot = 'counts'
)

seu_wgcna <- NormalizeMetacells(seu_wgcna)
m_obj <- GetMetacellObject(seu_wgcna) 

seu_wgcna <- SetDatExpr(
  seu_wgcna,
  group_name = 'Tumor', # the name of the group of interest in the group.by column
  group.by='region', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'Spatial', # using RNA assay
#   slot = 'data' # using normalized data
)
seu_wgcna <- TestSoftPowers(
  seu_wgcna,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(seu_wgcna)

dir.create('output_direction')
setwd('output_direction')
seurat_obj <- ConstructNetwork(
  seu_wgcna, soft_power=10, #set soft_power based on plot # 8 for seu_wgcna_all # 12 for HH Morans
  setDatExpr=FALSE,
  tom_name = 'Tumor',
  overwrite_tom=TRUE # name of the topoligical overlap matrix written to disk
)  #Takes long
# PlotDendrogram(seurat_obj, main='Region hdWGCNA Dendrogram')

pdf('Tumor_Spatial.dend.pdf',9,3)
PlotDendrogram(seurat_obj, main='Spatial hdWGCNA dendrogram of tumor region')
dev.off()

TOM <- GetTOM(seurat_obj)
# compute all MEs in the full single-cell dataset
seurat_obj <- ModuleEigengenes(
  seurat_obj,
  group.by.vars=NULL
)


# compute eigengene-based connectivity (kME):
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'region', group_name = 'Tumor', sparse=FALSE
)


# rename the modules
seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = "SM"
)

# harmonized module eigengenes:
hMEs <- GetMEs(seurat_obj)

# module eigengenes:
MEs <- GetMEs(seurat_obj, harmonized=FALSE)
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, hMEs)

# plot genes ranked by kME for each module
pkME <- PlotKMEs(seurat_obj, ncol=5)

ModuleNetworkPlot(
    seurat_obj, 
    outdir='ModuleNetworks_Spatial', # new folder name
    n_inner = 5, # number of genes in inner ring
    n_outer = 10, # number of genes in outer ring
    n_conns = Inf, # show all of the connections
    plot_size=c(13,13), # larger plotting area
    vertex.label.cex=1 # font size
)

modules <- GetModules(seurat_obj) %>% 
  subset(module != 'grey') %>% 
  mutate(module = droplevels(module))

hub_genes <- GetHubGenes(seurat_obj,n_hubs = 30) %>% .$gene_name
hub_lnc=intersect(hub_genes,lnc)
print(hub_lnc)
modules[modules$gene_name%in% hub_lnc,c('module','gene_name')]

write.table(modules,"modules.tumor.txt", sep="\t", row.names = TRUE, quote = FALSE)

library(Seurat)

# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)

# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)

# gene enrichment packages
library(enrichR)
library(GeneOverlap)

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)
dir.create('Enrichment')
setwd('Enrichment')
dbs <- c('GO_Biological_Process_2023','MSigDB_Hallmark_2020')

seurat_obj <- RunEnrichr(
  seurat_obj,
  dbs=dbs, # character vector of enrichr databases to test
  max_genes = 500 # number of genes per module to test
)

# retrieve the output table
enrich_df <- GetEnrichrTable(seurat_obj)
save(enrich_df,file='ML1_SM_enrich.Tumor.RData')

lapply(unique(modules$module[modules$gene_name%in% hub_lnc]),function(sm){
    p=EnrichrDotPlot(
    seurat_obj,
    wgcna_name='Tumor',
    mods = sm, # use all modules (this is the default behavior)
    database = 'GO_Biological_Process_2023', # this has to be one of the lists we used above!!!
    n_terms=10 # number of terms for each module
    )
    pdf(paste0('GO.',sm,'_Tumor.pdf'),8,12)
    print(p)
    dev.off()
    p=EnrichrDotPlot(
    seurat_obj,
    wgcna_name='Tumor',
    mods = sm, # use all modules (this is the default behavior)
    database = 'MSigDB_Hallmark_2020', # this has to be one of the lists we used above!!!
    n_terms=10 # number of terms for each module
    )
    pdf(paste0('H.',sm,'_Tumor.pdf'),8,12)
    print(p)
    dev.off()
})

saveRDS(seurat_obj,'my_output_name.rds')

seurat_obj@meta.data <- cbind(seurat_obj@meta.data, hMEs)
lapply(unique(modules$module[modules$gene_name%in% hub_lnc]),function(x){
    pdf(paste0(x,'_Tumor_Spatial.pdf'),5,5)
    print(SpatialFeaturePlot(
    seurat_obj,
    features = x,
    # restrict_range=FALSE
     pt.size.factor = 2
    ))
    dev.off()
    })



modules <- GetModules(seurat_obj)
mods <- levels(modules$module); mods <- mods[mods != 'grey']

# add the MEs to the seurat metadata so we can plot it with Seurat functions


p <- DotPlot(seurat_obj, features=mods, group.by = 'region', dot.min=0.1)

# flip the x/y axes, rotate the axis labels, and change color scheme:
p <- p +
  coord_flip() +
  RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue') +
  xlab('') + ylab('')

pdf(paste0(sample,'_region.Tdot.pdf'),4,12)
print(p)
dev.off()


# plot with Seurat's DotPlot function
p <- DotPlot(seurat_obj, features=mods, group.by = 'SpaCluster', dot.min=0.1)

# flip the x/y axes, rotate the axis labels, and change color scheme:
p <- p +
  coord_flip() +
  RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue') +
  xlab('') + ylab('')

pdf(paste0(sample,'_SpaCluster.Tdot.pdf'),15,12)
print(p)
dev.off()
