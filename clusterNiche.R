suppressPackageStartupMessages({
library(argparse)
})

parser <- ArgumentParser()
parser$add_argument('--outdir', help='outdir of result, [default %(default)s]', default='.')
parser$add_argument('--prefix', help='prefix of outdir, [default %(default)s]', default='program')
parser$add_argument('--abundance', help='abundance csv file for each sample, split by ,')
parser$add_argument('--method', help='scran or ISCHIA, [default %(default)s]', default='scran')
parser$add_argument('--sample', help='sample name, split by ,')
parser$add_argument('--outs', help='outs directory of sample, which include sptial. split by ,')
parser$add_argument('--chip_type', help='eq HD 3mm*3mm 3mm*3mm. split by ,')
args <- parser$parse_args()

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

createDir <- function(dir){
    if(!dir.exists(dir)){
        dir.create(dir, recursive=T)
    }
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

merge_sample <- function(abundance, sample){
    abundance_file <- unlist(strsplit(abundance, split=","))
    sample <- unlist(strsplit(sample, split=","))
    list_matrices <- map2(abundance_file, sample, function(f, s) {
        mat <- read.csv(f, head=T, row.names=1, check.names=F)
        rownames(mat) <- paste0(s, "_", rownames(mat))
        prop_mat <- base::apply(mat, 1, function(x) {x/sum(x)})
        return(t(prop_mat))
    })

    names(list_matrices) <- sample
    integrated_compositions <- purrr::reduce(list_matrices, rbind)
    return(integrated_compositions)
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

    pdf <- paste0(file.path(args$outdir, args$prefix), "_umap_cluster.pdf")
    plotFunc(p, pdf, pdfwidth=7*params$ratio)

    cts <- set_names(colnames(integrated_compositions))

    createDir(file.path(args$outdir, "celltype_proption"))
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
        pdf <- paste0(file.path(args$outdir, "celltype_proption", args$prefix), "_umap_", ct_string, "_feature.pdf")
        plotFunc(p, pdf, pdfwidth=7*params$ratio)
    })

}

preprocessing <- function(integrated_compositions){
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
        print("Louvain clustering")
#         clust.info <- igraph::cluster_louvain(snn_graph, resolution = 1.7)
        clust.info <- igraph::cluster_leiden(snn_graph, resolution_parameter = 0.55)
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


runISCHIA <- function(stRNA, norm_weights){
    # Generate ILR transformation
    baseILR <- ilrBase(x = norm_weights,
                   method = "basic")

    cell_ilr <- as.matrix(ilr(norm_weights, baseILR))
    colnames(cell_ilr) <- paste0("ILR_", 1:ncol(cell_ilr))

    comp_umap <- umap(cell_ilr, 
                  n_neighbors = 30, 
                  n_epochs = 1000,
                  metric = "cosine") %>%
        as.data.frame() %>%
        mutate(barcode = rownames(cell_ilr))

    comp_umap <- comp_umap[, 1:3]
    print(head(comp_umap))


    # Elbow Method
    k.values <- 1:20
    wss_values <- sapply(k.values, function(k) kmeans(norm_weights, k, nstart = 10)$tot.withinss)

    png(paste0(file.path(args$outdir, args$prefix), "_elbow_plot.png"))
    print(plot(k.values, wss_values, type = "b", pch = 19, frame = FALSE, 
        xlab = "Number of clusters K", ylab = "Total within-cluster sum of squares",
        main = "Elbow Method for Optimal K"))
    dev.off()

    # ISCHIA Analysis
    png(paste0(file.path(args$outdir, args$prefix), "_composition_cluster_k_plot.png"))
    print(Composition.cluster.k(norm_weights, 20))
    dev.off()

    k <- 10
    sample.spot.clusters <- kmeans(norm_weights, k)
    cluster <- as.numeric(as.vector(sample.spot.clusters$cluster))
    levels <- paste0('niche', sort(unique(cluster)))
    cluster <- factor(paste0('niche', cluster), levels=levels)
    cluster_info <- data.frame(barcode=names(sample.spot.clusters$cluster), niche=cluster)



    comp_umap <- comp_umap %>% left_join(cluster_info, by = c("barcode")) %>% dplyr::rename("UMAP1"="V1", "UMAP2"="V2")
    cluster_info <- comp_umap %>%
        dplyr::select(c(barcode, UMAP1, UMAP2)) %>% 
        left_join(cluster_info, by='barcode')
    
    write_csv(cluster_info, file = paste0(file.path(args$outdir, args$prefix), "_niche_annotation.csv"))

    return(cluster_info)
}

get_pt_size_factor <- function(stRNA, chip_type){
    chip_point <- c(0.65, 0.65, 1.25, 1.25, 0.65, 1.25, 1.25)
    names(chip_point) <- c(
        "HD 3mm*3mm", 
        "3mm*3mm", 
        "5mm*5mm", 
        "6mm*6mm", 
        "HD 6mm*6mm", 
        "7.5mm*7.5mm", 
        "10mm*10mm")

    if(chip_type %in% names(chip_point)){
        pt_size_factor <- chip_point[chip_type]
    }else{
        print(paste0("chip_type error: ", chip_type))
        pt_size_factor <- 1.25
    }

    print(paste0("pt_size_factor selected: ", pt_size_factor))
    return(pt_size_factor)
}

niche_spatial <- function(cluster_info, sample, outs, chip, colours){
    for(i in 1:length(sample)){
        stRNA <- Load10X_Spatial(data.dir = outs[i], slice = sample[i])
        stRNA <- RenameCells(object = stRNA, add.cell.id = sample[i])
        stRNA$orig.ident <- sample[i]
        meta <- data.frame(barcode=colnames(stRNA)) %>% 
            left_join(cluster_info) %>% 
            column_to_rownames("barcode") %>% 
            select(niche) %>% 
            filter(!is.na(niche))
        stRNA <- subset(stRNA, cells=rownames(meta))
        stRNA$niche <- meta$niche
        Idents(stRNA) <- stRNA$niche
        print(head(stRNA@meta.data))

        num <- length(unique(stRNA$orig.ident))
        pt.size.factor <- get_pt_size_factor(stRNA, chip[i])

        # stRNA@images[[sample[i]]]@spot.radius <- stRNA@images[[sample[i]]]@scale.factors$spot *  stRNA@images[[sample[i]]]@scale.factors$lowres
        if (num == 1){
            p <- SpatialPlot(stRNA, cols=colours, label = TRUE, label.size = 3, pt.size.factor = pt.size.factor, crop=FALSE) + 
                guides(fill = guide_legend(override.aes = list(size=4)))
            width_num <- 1
            height_num <- 1
        }else{
            p <- SpatialPlot(stRNA, cols=colours, label = TRUE, label.size = 3, ncol = 2, pt.size.factor = pt.size.factor, crop=FALSE) & 
                guides(fill = guide_legend(override.aes = list(size=4)))
            width_num <- 2
            height_num <- ceiling(num / 2)
        }
        ggsave(p, path = args$outdir, filename = paste0(args$prefix, '.niche.spatialDimPlot_', sample[i], '.png'), device = 'png', width=7*width_num, height=7*height_num)
        ggsave(p, path = args$outdir, filename = paste0(args$prefix, '.niche.spatialDimPlot_', sample[i], '.pdf'), device = 'pdf', width=7*width_num, height=7*height_num)
    }
}

niche_summary <- function(integrated_compositions, cluster_info){
    niche_summary_pat <- integrated_compositions %>%
        as.data.frame() %>%
        rownames_to_column("barcode") %>%
        pivot_longer(-barcode,values_to = "ct_prop", 
                    names_to = "cell_type") %>%
        left_join(cluster_info) %>%
        mutate(orig.ident = strsplit(barcode, "[..]") %>%
                map_chr(., ~ .x[1]))# %>%
        # group_by(orig.ident, niche, cell_type) %>%
        # summarize(median_ct_prop = median(ct_prop))

    niche_summary <- niche_summary_pat %>%
        ungroup() %>%
        group_by(niche, cell_type) %>%
        summarise(patient_median_ct_prop = median(ct_prop))

    # Data manipulation to have clustered data
    niche_summary_mat <- niche_summary %>%
        pivot_wider(values_from = patient_median_ct_prop, 
                    names_from =  cell_type, values_fill = 0) %>%
        column_to_rownames("niche") %>%
        as.matrix()

    niche_order <- hclust(dist(niche_summary_mat))
    niche_order <- niche_order$labels[niche_order$order]

    ct_order <- hclust(dist(t(niche_summary_mat)))
    ct_order <- ct_order$labels[ct_order$order]

    wilcoxon_res <- niche_summary_pat %>%
        ungroup() %>%
        group_by(cell_type) %>%
        nest() %>%
        mutate(wres = map(data, run_wilcox_up)) %>%
        dplyr::select(wres) %>%
        unnest() %>%
        ungroup() %>%
        mutate(p_corr = p.adjust(p.value))

    wilcoxon_res <- wilcoxon_res %>%
        mutate(significant = ifelse(p_corr <= 0.15, "*", ""))

    write.table(niche_summary_pat, file = paste0(file.path(args$outdir, args$prefix), "_niche_summary.txt"), 
                col.names = T, row.names = F, quote = F, sep = "\t")

    write.table(wilcoxon_res, file = paste0(file.path(args$outdir, args$prefix), "_wilcoxon_res_cells_niches.txt"), 
                col.names = T, row.names = F, quote = F, sep = "\t")

    mean_ct_prop_plt <- niche_summary %>%
        left_join(wilcoxon_res, by = c("niche", "cell_type")) %>%
        mutate(cell_type = factor(cell_type, levels = ct_order),
                niche = factor(niche, levels = niche_order)) %>%
        ungroup() %>%
        group_by(cell_type) %>%
        mutate(scaled_pat_median = (patient_median_ct_prop - mean(patient_median_ct_prop))/sd(patient_median_ct_prop)) %>%
        ungroup() %>%
        ggplot(aes(x = cell_type, y = niche, fill = scaled_pat_median)) +
        geom_tile() +
        labs(x=NULL, y=NULL, fill="scaled prop median") +
        geom_text(aes(label = significant)) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
                legend.position = "bottom",
                legend.title = element_text(margin = margin(0, 15, 0, 0), vjust = 0.5), 
                plot.margin = unit(c(0.2, 0, 0.2, 0.2), "cm"),
                axis.text.y = element_text(size=12)) +
        scale_fill_gradient(high = "#ffd89b", low = "#19547b")
    

    # Finally describe the proportions of those niches in all the data
    cluster_counts <- cluster_info %>%
        dplyr::select_at(c("barcode", "niche")) %>%
        group_by(niche) %>%
        summarise(nspots = length(niche)) %>%
        mutate(prop_spots = nspots/sum(nspots))

    write.table(cluster_counts, file = paste0(file.path(args$outdir, args$prefix), "_niche_prop_summary.txt"),
        col.names = T, row.names = F, quote = F, sep = "\t")

    barplts <- cluster_counts %>%
        mutate(niche = factor(niche, levels = niche_order)) %>%
        ggplot(aes(y = niche, x = prop_spots)) +
        geom_bar(stat = "identity") +
        theme_classic() + ylab("") +
        theme(
            # axis.text.y = element_blank(),
                plot.margin = unit(c(0.2, 0.4, 0.2, 0), "cm"),
                axis.text.x = element_text(size=12)) 

    niche_summary_plt <- cowplot::plot_grid(mean_ct_prop_plt, barplts, align = "hv", axis = "tb")

    pdf <- paste0(file.path(args$outdir, args$prefix), "_characteristic_niches.pdf")
    plotFunc(niche_summary_plt, pdf, pdfwidth = 8, pdfheight = 4)

    # Show the compositions of cells of each niche 
    p <- niche_summary_pat %>%
        ggplot(aes(x = niche, y = ct_prop)) +
        geom_boxplot(outlier.size = 0.2) + 
        labs(x=NULL) + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
        facet_wrap(.~cell_type, ncol = 3,scales = "free_y")
    pdf <- paste0(file.path(args$outdir, args$prefix), "_niche_summary.pdf")
    plotFunc(p, pdf, pdfwidth = 10, pdfheight = 6)
}

celltypeCoocurence <- function(stRNA, integrated_compositions){

    createDir(file.path(args$outdir, "celltype_cooccurrence"))
    # Cell type co-occurrence
    for(niche in levels(stRNA$niche)){
        print(niche)
        niche_spot <- stRNA$niche[stRNA$niche==niche]
        if(length(niche_spot)==1)next
        celltype.cooccur <- spatial.celltype.cooccurence(spatial.object=stRNA, deconv.prob.mat=integrated_compositions, 
                                                            COI=niche, prob.th = 0.05, 
                                                            Condition=unique(stRNA$orig.ident))
        p <- plot.celltype.cooccurence(celltype.cooccur)
        pdf <- paste0(file.path(args$outdir, "celltype_cooccurrence", args$prefix), "_niche_celltype_cooccurrence_", niche, ".pdf")
        plotFunc(p, pdf, pdfwidth=15, pdfheight=12)
    }
}

source("/disk/pipeline/DynaSpatail/SpatialPlot.R")

colours <- c("#00CC99", "#F0E685", "#99CC00", "#C75127", "#5050FF", "#CDDEB7","#EE82EE",
"#924822","#FFC20A","#33CC00","#0099CC","#7A65A5","#AE1F63","#D595A7",
"#A9A9A9","#E4AF69","#990080","#14FFB1","#003399", "#FF1463","#D60047",
"#CE3D32","#996600","#809900","#008099","#0A47FF","#660099","#FFD147",
"#339900","#5A655E","#4775FF","#802268","#8B0000","#BBFFFF","#00868B",
"#00008B","#CDBA96","#FF69B4","#3B1B53","#D2B48C","#00FFFF","#836FFF",
"#BA6338","#991A00","#D58F5C","#E7C76F","#CC9900","#749B58","#6BD76B",
"#009966","#00D68F","#5DB1DD","#466983","#837B8D","#612A79","#990033"
)

main <- function(){
    createDir(args$outdir)
    integrated_compositions <- merge_sample(args$abundance, args$sample)

    sample <- unlist(strsplit(args$sample, split=","))
    outs <- unlist(strsplit(args$outs, split=","))
    chip <- unlist(strsplit(args$chip, split=","))
    stRNA <- mergeSeurat(sample, outs)

    if(args$method == 'ISCHIA'){
        cluster_info <- runISCHIA(stRNA, integrated_compositions)
    }else{
        cluster_info <- preprocessing(integrated_compositions)
    }
    write.csv(cluster_info, file=paste0(file.path(args$outdir, args$prefix), "_niche_annotation.csv"), row.names=F)

    cluster <- cluster_info %>% select(barcode, niche) %>% column_to_rownames('barcode')
    stRNA <- subset(stRNA, cells=rownames(cluster))
    stRNA <- AddMetaData(stRNA, cluster, "CompositionCluster_CC")
    stRNA$niche <- stRNA$CompositionCluster_CC

    names(colours) <- levels(cluster_info$niche)
    print(colours)
    visNiche(integrated_compositions, cluster_info, colours)
    niche_spatial(cluster_info, sample, outs, chip, colours)
    niche_summary(integrated_compositions, cluster_info)
    celltypeCoocurence(stRNA, integrated_compositions)
}
main()
