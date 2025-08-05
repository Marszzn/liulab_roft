setwd("/mnt/seagate1/roft_v2/final_code")
# load("spatial_shiny.RData")
packages <- c(
  "SeuratObject", "shiny", "shinyjs", "shinythemes", "Seurat",
  "plotly", "patchwork", "cowplot", "gridExtra", "RColorBrewer",
  "ggplot2", "grImport", "plyr", "ggefp", "DT", "stringi",
  "gtools", "pheatmap", "igraph", "networkD3", "WGCNA",
  "dplyr", "CellChat", "SeuratDisk", "PlantPhoneDB", "viridis",
  "topGO", "visNetwork"
)
invisible(lapply(packages, require, character.only = TRUE))
source("/mnt/seagate1/roft_v2/final_code/roft_function.R")
st0 <- readRDS("st0-2.single_seruat.Rds")
st2 <- readRDS("st1-1.single_seruat.Rds")
st3 <- readRDS("st3-2.single_seruat.Rds")
st5 <- readRDS("St5.single_seruat.Rds")
# load("seurat_with_annotation_final.RData")

#====geid rename=====================================================
st0 <- UpdateSeuratObject(st0)
st2 <- UpdateSeuratObject(st2)
st3 <- UpdateSeuratObject(st3)
st5 <- UpdateSeuratObject(st5)

st0 <- sanitize_gene_names(st0)
st2 <- sanitize_gene_names(st2)
st3 <- sanitize_gene_names(st3) 
st5 <- sanitize_gene_names(st5)
# Seurat v4 的对象更新为当前 Seurat 版本（v5）的对象结构；统一版本===

# spatial cluster annotation==========================================
st0.anno <- data.frame(V1 = sub("-1$", "", Cells(st0)),V2 = paste0("cluster", st0$seurat_clusters))
st2.anno <- read.delim("st2.annotation.txt",header = F)
st3.anno <- read.delim("st3.annotation.txt",header = F)
st5.anno <- read.delim("st5.annotation.txt",header = F)

# cell-type annotation (barcode)
st0 <- add_annotation_to_seurat(st0, st0.anno, anno_colname = "annotation")
st2 <- add_annotation_to_seurat(st2, st2.anno, anno_colname = "annotation")
st3 <- add_annotation_to_seurat(st3, st3.anno, anno_colname = "annotation")
st5 <- add_annotation_to_seurat(st5, st5.anno, anno_colname = "annotation")
rm(st0.anno,st2.anno,st3.anno,st5.anno)
# load("seurat_with_annotation_final.RData")

#======================CCI=========================================
#=================================UMAP_PLOT=======================
st_list <- list(st0=st0,st2=st2,st3=st3,st5=st5)
# all_annotations <- sort(unique(unlist(lapply(st_list, function(obj) unique(obj$annotation)))))

for (name in names(st_list)) {
  obj <- st_list[[name]]
  # obj <- st5
  p <- SpatialDimPlot(obj, 
                      label = TRUE, 
                      group.by = "annotation", 
                      label.size = 2, 
                      label.color = "gray90", 
                      pt.size = 300, 
                      image.alpha = 0.4) +
    scale_fill_manual(values = mycolor) +
    ggtitle(paste0(name, "_Cell-type"))
  
  ggsave(filename = file.path("/mnt/seagate1/roft_v2/final_code/result/", paste0(name, "_celltype.jpg")),
         plot = p, device = "jpeg", dpi = 600, width = 6, height = 6, units = "in")
}
rm(p,obj)
#=====LR_combined=========================================================
#=================== Step 1: 原始数据 ===================
load("LR_pair_ath.RDa")  # 含 Arabidopsis 的配对信息 LR_pair
blast_unique <- read.csv("fvh_blast_ara.csv")  # Fragaria -> Arabidopsis 的基因映射

#=================== Step 2: 匹配 Ligands and Receptors====================
idx_ligands <- match(blast_unique$ID2, LR_pair$Ligands)
ligands <- blast_unique[!is.na(idx_ligands), ]

idx_receptors <- match(blast_unique$ID2, LR_pair$Receptors)
receptors <- blast_unique[!is.na(idx_receptors), ]

#=================== Step 3: 构建 Fragaria vesca LR配对,以及LR score =========
LR_combined <- merge(LR_pair, ligands, by.x = "Ligands", by.y = "ID2", all.x = FALSE) %>%
  merge(receptors, by.x = "Receptors", by.y = "ID2", all.x = FALSE, suffixes = c("_ligand", "_receptor")) %>%
  dplyr::select(-Receptors, -Ligands) %>%
  dplyr::rename(Receptors = ID1_ligand,
         Ligands = ID1_receptor) %>%
  dplyr::mutate(Organism = "Fragaria vesca") %>%
  dplyr::select(Ligands, Receptors, source, Organism)
rm(objs,ligands,receptors)
# objs<- objs_list[[1]]
# cluster <- objs$annotation
# names(cluster) <- colnames(objs)
# score <- LRscore(objs@assays$SCT@data, LRdb=LR_combined, 
# cluster = cluster, min.pct = 0.1,iterations=100, method='Average')
LRscore_list <- list()
for (name in names(st_list)) {
  message("Processing: ", name)
  objs <- st_list[[name]]
  cluster <- objs$annotation
  names(cluster) <- colnames(objs)
  # 计算 LRscore
  score <- LRscore(
    objs@assays$SCT@data,
    LRdb = LR_combined,
    cluster = cluster,
    min.pct = 0.1,
    iterations = 100,
    method = 'Average'
  )
  LRscore_list[[name]] <- score
}
LR_sig <- lapply(names(LRscore_list), function(name) {
  df <- LRscore_list[[name]]
  df <- df %>% filter(Pvalue < 0.05)
  df$Sample <- name
  df
})
# save.image(LR_sig,LRscore_list,file="LR_sig0801.rdata")
rm(df,score,LR_combined,LR_pair)
names(LR_sig) <- names(LRscore_list)
# save(LR_sig, LRscore_list, file = "LR_sig0801.rdata")
# load("LR_sig0801.rdata")


#==========UMAP_plot===================================================
for (name in names(st_list)) { 
  obj <- st_list[[name]]
  p <- DimPlot(obj,group.by = 'annotation', 
               label=TRUE, 
               label.size = 3,
               label.color = "white", label.box = T, reduction='umap',
                cols=mycolor)+
  NoLegend()+ggtitle(paste0(name,"_UMAP_plot"))
  ggsave(filename = file.path("/mnt/seagate1/roft_v2/final_code/result/", paste0(name, "_UMAP_plot.jpg")),
         plot = p, device = "jpeg", dpi = 600, width = 7, height = 7, units = "in")
  rm(obj,p)
  }

#================barplot==============================================
meta_all <- bind_rows(lapply(names(st_list), function(name) {
  df <- st_list[[name]]@meta.data
  df$Stage <- name
  return(df)
}))

barplot_combined <- ggplot(meta_all, aes(x = Stage, fill = annotation)) +
  geom_bar(position = "fill", color = "white") + 
  coord_flip() +
  theme_bw() +
  ylab("") +
  ggplot2::theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 12, color = 'black'),
    legend.text = element_text(size = 8, face = "bold"),
    legend.title = element_blank(),
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),  # 加居中和放大
    axis.line = element_line(color = 'black')
  ) +
  scale_y_continuous(position = "right", expand = c(0, 0)) +
  scale_fill_manual(values = mycolor) +
  ggtitle("Cluster Proportion across Samples")

ggsave(
  filename = file.path("/mnt/seagate1/roft_v2/final_code/result/", "spatial_celltype_barplot.jpg"),
  plot = barplot_combined,
  device = "jpeg",
  dpi = 600,
  width = 9,
  height = 4.5,
  units = "in"
)
rm(meta_all,barplot_combined)
#=========================stage_image & gene expression=====================
# gene
SpatialFeaturePlot(
  st_list$st5,
  slot = "data",
  # layer = "data", 
  features = "FvH4_1g01140",
  pt.size.factor = 350,
  image.alpha = 0.3) +
  ggtitle("gene_expression") +
  scale_fill_gradientn(colors = c("#6238ff",  "#f2e9e1", "#ff220e"))
# cell 
VlnPlot(
  st_list$st5,
  group.by = "annotation",
  assay = "Spatial",
  features = "FvH4_1g01140",
  pt.size = 1
) + scale_fill_manual(values = mycolor) 

#============= find marker and GO function======================================================
# P-value < 0.5 marker gene
st_marker_list <- list()
for (names in names(st_list)) {
  obj <- st_list[[names]]
  markers <- FindAllMarkers(obj,
                        group.by = "annotation",
                        logfc.threshold=0.25,
                        min.diff.pct = 0.25,
                        max.cells.per.ident = Inf,
                        only.pos=F,
                        return.thresh = 0.01)
  st_marker_list[[names]] <- markers
}
table(st_marker_list[["st5"]][,6])

go_marker_list <- list()
for (sample_name in names(st_marker_list)) {
  cat("Processing sample:", sample_name, "\n")
  marker_df <- st_marker_list[[sample_name]]
  
  clusters <- unique(marker_df$cluster)
  sample_go_list <- list()
  
  for (clust in clusters) {
    gene_list <- unique(marker_df$gene[marker_df$cluster == clust & marker_df$p_val_adj < 0.05])
    
    cluster_go_list <- list()
    
    if (length(gene_list) > 1) {
      for (ont in c("BP", "MF", "CC")) {
        cat("  → Cluster:", clust, "| Ontology:", ont, "\n")
        go_table <- tryCatch({
          get_GO_strawberry(gene_list, type = ont)
        }, error = function(e) {
          warning(paste("Error in", sample_name, clust, ont, ":", e$message))
          return(NULL)
        })
        cluster_go_list[[ont]] <- go_table
      }
    } else {
      cluster_go_list <- list(BP = NULL, MF = NULL, CC = NULL)
    }
    
    sample_go_list[[as.character(clust)]] <- cluster_go_list
  }
  
  go_marker_list[[sample_name]] <- sample_go_list
}
rm(sample_go_list,obj,cluster_go_list,markers,marker_df,go_table)



#====cluster2cluster compare=================================
st_c2c_marker_list <- list()
for (name in names(st_list)) {
  message("Processing: ", name)
  obj <- st_list[[name]]
  st_c2c_marker_list[[name]] <- get_cluster_markers(obj, sample_name = name)
}

st_c2c_marker_list <- lapply(st_c2c_marker_list, function(stage_list) {
  lapply(stage_list, function(df) {
    df$Gene.ID <- rownames(df)
    df <- df[, c("Gene.ID", setdiff(names(df), "Gene.ID"))]
    rownames(df) <- NULL
    return(df)
  })
})
rm(obj,clust,clusters,combined_df,gene_list,name,names,ont,sample_name)
#=======GO cluster2cluster=====================================================
go_c2c_list <- lapply(st_c2c_marker_list, function(stage_list) {
  lapply(stage_list, function(df) {
    gene_list <- df$Gene.ID
    go_types <- c("BP", "MF", "CC")
    
    go_tables <- lapply(go_types, function(go_type) {
      if (length(gene_list) > 1) {
        get_GO_strawberry(gene_list, type = go_type)
      } else {
        NULL
      }
    })
    names(go_tables) <- go_types
    return(go_tables)
  })
})

save.image(file = "roft_v1_0804.rdata")
#===============================================================================
