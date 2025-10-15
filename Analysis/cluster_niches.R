suppressPackageStartupMessages({
library(compositions)
library(Seurat)
library(tidyverse)
library(clustree)
library(uwot)
library(scran)
library(cluster)
library(patchwork)
library(ISCHIA)
library(robustbase)
library(data.table)
})

colours <- c("#00CC99", "#F0E685", "#99CC00", "#C75127", "#5050FF", "#CDDEB7","#EE82EE",
"#924822","#FFC20A","#33CC00","#0099CC","#7A65A5","#AE1F63","#D595A7",
"#A9A9A9","#E4AF69","#990080","#14FFB1","#003399", "#FF1463","#D60047",
"#CE3D32","#996600","#809900","#008099","#0A47FF","#660099","#FFD147",
"#339900","#5A655E","#4775FF","#802268","#8B0000","#BBFFFF","#00868B",
"#00008B","#CDBA96","#FF69B4","#3B1B53","#D2B48C","#00FFFF","#836FFF",
"#BA6338","#991A00","#D58F5C","#E7C76F","#CC9900","#749B58","#6BD76B",
"#009966","#00D68F","#5DB1DD","#466983","#837B8D","#612A79","#990033"
)

createDir <- function(dir){
    if(!dir.exists(dir)){
        dir.create(dir, recursive=T)
    }
}

merge_sample <- function(abundance, sample){
    list_matrices <- map2(abundance, sample, function(f, s) {
        mat <- read.csv(f, head=T, row.names=1, check.names=F)
        rownames(mat) <- paste0(s, "_", rownames(mat))
        prop_mat <- base::apply(mat, 1, function(x) {x/sum(x)})
        return(t(prop_mat))
    })

    names(list_matrices) <- sample
    integrated_compositions <- purrr::reduce(list_matrices, rbind)
    return(integrated_compositions)
}

mergeSeurat <- function(sample, outs){
    scRNA_obj_list <- list()
    for(i in 1:length(sample)){
        stRNA <- Load10X_Spatial(data.dir = outs[i], slice = sample[i])
        stRNA <- RenameCells(object = stRNA, add.cell.id = sample[i])
        stRNA$orig.ident <- sample[i]
        print(head(stRNA@meta.data))
        scRNA_obj_list[[i]] <- stRNA
    }
    if(length(sample)==1){
        return(stRNA)
    }
    stRNA <- merge(scRNA_obj_list[[1]], scRNA_obj_list[2:length(scRNA_obj_list)], project = 'Integration')
}

plotFunc <- function(pic, pdf, pdfwidth=6, pdfheight=6, png=TRUE, res=200){
    pdf(pdf, width=pdfwidth, height=pdfheight, bg='white')
    print(pic)
    dev.off()
    if(png){
        png <- gsub('pdf$', 'png', pdf)
        png(png, width=pdfwidth, height=pdfheight, res = res, unit="in")
        print(pic)
        dev.off()
    }
}
preprocessing <- function(integrated_compositions,resol){
    # Generate ILR transformation
    baseILR <- ilrBase(x = integrated_compositions,
                   method = "basic")

    cell_ilr <- as.matrix(ilr(integrated_compositions, baseILR))
    colnames(cell_ilr) <- paste0("ILR_", 1:ncol(cell_ilr))

    # Make community graph
    k_vect <- c(10,20,40)
    k_vect <- set_names(k_vect, paste0("k_",k_vect))

    cluster_info <- map(k_vect, function(k) {
        print(k)
        print("Generating SNN")
        snn_graph <- scran::buildSNNGraph(x = t(cell_ilr %>% as.data.frame() %>% as.matrix()), k = k)
        print("clustering")
#         clust.info <- igraph::cluster_louvain(snn_graph, resolution = 1.7)
        clust.info <- igraph::cluster_leiden(snn_graph, resolution_parameter = resol)
        print(length(unique(clust.info$membership)))
        clusters <- tibble(cluster = clust.info$membership,
                            spot_id = rownames(cell_ilr))
    })

    cluster_info <- cluster_info %>% 
        enframe() %>%
        unnest() %>%
        pivot_wider(names_from = name,
                    values_from = cluster)

    set.seed(2024)

    k_vect <- set_names(names(k_vect))

    subsampling_mean_ss <- map(k_vect, function(k) {
        print(k)
        
        cluster_info_summary <- cluster_info %>%
            group_by_at(k) %>%
            summarize(ncells = floor(n() * 0.3))
        
        cells <- cluster_info %>%
            dplyr::select_at(c("spot_id", k)) %>%
            group_by_at(k) %>%
            nest() %>%
            left_join(cluster_info_summary) %>%
            mutate(data = map(data, ~ .x[[1]])) %>%
            mutate(selected_cells = map2(data, ncells, function(dat,n) {
            sample(dat, n)
            })) %>%
            pull(selected_cells) %>%
            unlist()
        
        dist_mat <- dist(cell_ilr[cells, ])
        k_vect <- purrr::set_names(cluster_info[[k]], cluster_info[["spot_id"]])[cells]
        sil <- cluster::silhouette(x = k_vect, dist = dist_mat)
        mean(sil[, 'sil_width'])
    })

    subsampling_mean_ss <- enframe(subsampling_mean_ss) %>% 
        unnest() %>%
        dplyr::filter()

    niche_resolution <- dplyr::filter(subsampling_mean_ss, value == max(value)) %>% pull(name)
    print(niche_resolution)

    comp_umap <- umap(cell_ilr, 
                  n_neighbors = 30, 
                  n_epochs = 1000,
                  metric = "cosine") %>%
        as.data.frame() %>%
        mutate(barcode = rownames(cell_ilr))

    comp_umap <- comp_umap[, 1:3]
    print(head(comp_umap))
    comp_umap <- comp_umap %>% left_join(cluster_info, by = c("barcode" = "spot_id")) %>% dplyr::rename("UMAP1"="V1", "UMAP2"="V2")

    cluster_info <- comp_umap %>%
        dplyr::select(c("barcode", !!niche_resolution, UMAP1, UMAP2)) %>%
        dplyr::rename("niche" = !!niche_resolution)
    levels <- paste0("niche", sort(unique(as.numeric(cluster_info$niche))))
    cluster_info$niche <- factor(paste0("niche", cluster_info$niche), levels=levels)

    return(cluster_info)
}
run_wilcox_up <- function(prop_data) {
  
  prop_data_group <- prop_data[["niche"]] %>%
    unique() %>%
    set_names()
  
  map(prop_data_group, function(g) {
    
    test_data <- prop_data %>%
      mutate(test_group = ifelse(niche == g,
                                 "target", "rest")) %>%
      mutate(test_group = factor(test_group,
                                 levels = c("target", "rest")))
    
    wilcox.test(ct_prop ~ test_group, 
                data = test_data,
                alternative = "greater") %>%
      broom::tidy()
  }) %>% enframe("niche") %>%
    unnest()
}

reductionParams <- function(cluster){
    cluster_num <- length(cluster)
    if(NA %in% (as.numeric(cluster))){
        if(cluster_num>20){
            ncol <- ceiling(cluster_num/20)
            ratio <- 1.8
        }else{
            ncol <- 1
            ratio <- 1.6
        }
    }else{
        if(cluster_num>20){
            ncol <- ceiling(cluster_num/20)
            ratio <- 12/10
        }else{
            ncol <- 1
            ratio <- 12/11
        }
    }
    return(list(ncol=ncol, ratio=ratio))
}

visNiche <- function(integrated_compositions, cluster_info, colours){
    log_comps <- log10(integrated_compositions)
    params <- reductionParams(unique(cluster_info[, 'niche']))

    p <- cluster_info %>%
        ggplot(aes(x = UMAP1, y = UMAP2, color = niche)) +
            ggrastr::geom_point_rast(size = 0.1) +
            theme_classic() + 
            theme(aspect.ratio=1) + 
            labs(x="UMAP1", y="UMAP2", color="niche") + 
            scale_color_manual(values=colours) + 
            guides(color = guide_legend(override.aes = list(size=4), ncol=params$ncol))

    pdf <- paste0(file.path(outdir, prefix), "_umap_cluster.pdf")
    plotFunc(p, pdf, pdfwidth=7*params$ratio)

    cts <- set_names(colnames(integrated_compositions))

    createDir(file.path(outdir, "celltype_proption"))
    walk(cts, function(ct){
    
        plot_df <- cluster_info %>%
            mutate(ct_prop = log_comps[ , ct])
        
        p <- plot_df %>%
            ggplot(aes(x = UMAP1, y = UMAP2, color = ct_prop)) +
            ggrastr::geom_point_rast(size = 0.1) +
            theme_classic() + 
            theme(aspect.ratio=1, plot.title=element_text(hjust=0.5)) + 
            labs(x="UMAP1", y="UMAP2", color="proption", title=ct) + 
            scale_color_viridis()
        
        ct_string <- str_replace_all(ct, " ", "_")
        pdf <- paste0(file.path(outdir, "celltype_proption", prefix), "_umap_", ct_string, "_feature.pdf")
        plotFunc(p, pdf, pdfwidth=7*params$ratio)
    })

}

outdir='my_direction'
prefix='S10'
resol='0.6'
main <- function(){
    
    createDir(outdir)
    sample <- c('sample_1','sample_n')
    abundance <- c(
        "sample_1_CARD.Proportion.csv",
        "sample_n_CARD.Proportion.csv"
    )
    integrated_compositions <- merge_sample(abundance, sample)
    
    outs <- c(
        'sample_1_SRT/',
        'sample_n_SRT/')
    chip=rep('5mm*5mm',10)   
    stRNA <- mergeSeurat(sample, outs)
    method = 'scran'
    if(method == 'ISCHIA'){
        cluster_info <- runISCHIA(stRNA, integrated_compositions)
    }else{
        cluster_info <- preprocessing(integrated_compositions,resol)
    }
    write.csv(cluster_info, file=paste0(file.path(outdir,prefix), "_niche_annotation.csv"), row.names=F)

    cluster <- cluster_info %>% select(barcode, niche) %>% column_to_rownames('barcode')
    stRNA <- subset(stRNA, cells=rownames(cluster))
    stRNA <- AddMetaData(stRNA, cluster, "CompositionCluster_CC")
    stRNA$niche <- stRNA$CompositionCluster_CC

    names(colours) <- levels(cluster_info$niche)
    
    print(colours)
    visNiche(integrated_compositions, cluster_info, colours)
    niche_spatial(cluster_info, sample, outs, chip, colours)
    niche_summary(integrated_compositions, cluster_info)

}
main()
