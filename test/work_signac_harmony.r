library(Seurat,lib.loc = "/home/mjwang/.conda/envs/myenv/lib/R/library_opt/") #3.1.5
#library(Signac) #0.2.4
library(Signac, lib.loc = "/home/mjwang/.conda/envs/myenv/lib/R/library_opt/") #1.0

#library(Seurat)
library(GenomeInfoDb)
#library(EnsDb.Hsapiens.v75)
library(EnsDb.Hsapiens.v86) #hg38, GRCh38, https://rdrr.io/bioc/ensembldb/f/vignettes/ensembldb.Rmd
library(ggplot2)
set.seed(1234)
library(Matrix)
library(dplyr)
library(reshape2)
library(magrittr)
library(plyr)
library(viridis)
#library(matrixStats)
library(ComplexHeatmap)
library(viridis)
library('ggplot2')


library(harmony)



sample <- "PLA-8w-ATAC"
#sample <- "PLA-8w-ATAC-new-force7025"
#sample <- "PLA-8w-ATAC-2"
#sample <- "PLA-25w-ATAC-1"
#sample <- "PLA-213-ATAC"
#sample <- "PLA-month-ATAC-1"
#sample <- "PLA-month-ATAC-2"

color <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c',
           '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', 
           '#000075', '#808080', '#000000')

color1 <- c("#FFFFFF" ,"#1E1E1E" ,"#E7D654", "#6F1482" ,"#DC7035", "#9CC8E6", "#BC2439" ,"#C2C281" ,"#7F8084", "#63AC4E", "#D181B0" ,"#476DAD",
            "#DF8969", "#49188D", "#E2AC3B" ,"#8F1287" ,"#E7F06F" ,"#7E191A" ,"#95B44F", "#6C3719" ,"#CA362E" ,"#2B3918")

#21 colors darker but some colors near
color_good <- c("#E7D654", "#6F1482" ,"#DC7035", "#9CC8E6", "#BC2439" ,"#C2C281" ,"#7F8084", "#63AC4E", "#D181B0" ,
                "#476DAD","#DF8969", "#49188D", "#E2AC3B" ,"#8F1287" ,"#E7F06F" ,"#7E191A" ,"#95B44F", "#6C3719" ,
                "#CA362E" ,"#2B3918","#1E1E1E" )

color.red.mono <- colorRampPalette(colors = c("white","red"))(10)




library(BuenColors)
color_set0 <- jdb_color_maps #17 different colors
names(color_set0) <- NULL
#plot(1:17,1:17,pch = 19, cex = 5,col=jdb_color_maps)

#discrete colors
color_set1 <- jdb_palette("solar_extra") #9 discrete but gradient colors
color_set2 <- jdb_palette("brewer_spectra") #9 discrete but gradient colors
color_set3 <- jdb_palette("flame_light") #9 discrete but gradient colors, good!

color_set3_ext12 <- colorRampPalette(colors = as.character(color_set3))(12)
color_set3_ext17 <- colorRampPalette(colors = as.character(color_set3))(17)

color_snap = c('1'='grey','2'='#E31A1C','3'='#FFD700','4'='#771122','5'='#777711','6'='#1F78B4','7'='#68228B','8'='#AAAA44','9'='#60CC52','10'='#771155','11'='#DDDD77','12'='#774411','13'='#AA7744','14'='#AA4455','15'='#117744')
names(color_snap) <- NULL

color_signac = c(
'0'='#E6D55E','1'='#792B8A','2'='#DA703D','3'='#9DC8E5','4'='#BA273C','5'='#C2C184','6'='#7F8084','7'='#65AB53','8'='#D082AF','9'='#496EAB','10'='#DE896D','11'='#491F8B','12'='#E1AD49','13'='#8E1B85','14'='#E7EE77','15'='#7D1A1D','16'='#96B355')
names(color_signac) <- NULL

color_good <- c("#E7D654", "#6F1482" ,"navy", "#9CC8E6", "#BC2439" ,"#C2C281" ,"#7F8084", 
                "#63AC4E", "#D181B0" ,"#476DAD","#DF8969", "#49188D", "#E2AC3B" ,"#8F1287" ,
                "#E7F06F" ,"#7E191A" ,"#95B44F", "#6C3719" ,"#CA362E" ,"#2B3918","#1E1E1E" )


color <- color_good






#################read in cellranger results and generate Seurat object
#fragments <- "test_data/atac_v1_pbmc_10k_fragments.tsv.gz"
#filter.mat.h5 <- "test_data/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5"
#singlecell <- "test_data/atac_v1_pbmc_10k_singlecell.csv"

# fragments <- "../01.data_cellranger_atac/PLA-8w-ATAC-new-force7025/fragments.tsv.gz"
# filter.mat.h5 <- "../01.data_cellranger_atac/PLA-8w-ATAC-new-force7025/filtered_peak_bc_matrix.h5"
# singlecell <- "../01.data_cellranger_atac/PLA-8w-ATAC-new-force7025/singlecell.csv"


# fragments <- "../01.aggregation/PLA-8w-ATAC-aggre/fragments.tsv.gz"
# filter.mat.h5 <- "../01.aggregation/PLA-8w-ATAC-aggre/filtered_peak_bc_matrix.h5"
# singlecell <- "../01.aggregation/PLA-8w-ATAC-aggre/singlecell.csv"

fragments <- "../01.aggregation/PLA-8w-ATAC-aggre_nonorm/fragments.tsv.gz"
filter.mat.h5 <- "../01.aggregation/PLA-8w-ATAC-aggre_nonorm/filtered_peak_bc_matrix.h5"
singlecell <- "../01.aggregation/PLA-8w-ATAC-aggre_nonorm/singlecell.csv"


### readin peak annotation files #fread (fast read), readr
##peakAnno <- read.table(paste("../../01.data_cellranger_atac/",sample,"/peak_annotation.tsv",sep=''),header=T,sep='\t',stringsAsFactors = FALSE)
# separate each row into a single peak-gene-type combination, i.e. split by ";"
peakAnno <- readr::read_tsv("../01.aggregation/PLA-8w-ATAC-aggre_nonorm/peak_annotation.tsv") %>%
  tidyr::separate_rows(gene, distance, peak_type, sep = ';')
peakAnno <- peakAnno[!is.na(peakAnno$gene),]
###


###use cisTopic to get a dgCMatrix, when you did not have .h5 file.

# peak_bc_mat <- paste("../01.data_cellranger_atac/PLA-8w-ATAC-new-force7025/filtered_peak_bc_matrix/")
# metrics <- singlecell
# cisTopicObject <- cisTopic::createcisTopicObjectFrom10Xmatrix(peak_bc_mat, metrics,project.name = sample) #quick but use 10X filter method 
# counts <- as(cisTopicObject@count.matrix,'dgCMatrix') #dgTMatrix to dgCMatrix


counts <- Read10X_h5(filename = filter.mat.h5) #195098 x 14441

metadata <- read.csv(
  file = singlecell,
  header = TRUE,
  row.names = 1
)


chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = fragments,
  min.cells = 10,
  min.features = 200
)

#metadata[colnames(counts),]

##create objects
placenta <- CreateSeuratObject(
  counts = chrom_assay,
  assay = 'peaks',
  project = 'ATAC',
  #min.cells = 1,
  meta.data = metadata
)


# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86) ##cellranger atac 1.1 use genecode 28, comparable to ensembl 92, but rbioconductor the latest is EnsDb.Hsapiens.v86

# change to UCSC style since the data was mapped to hg38
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

# add the gene information to the object
Annotation(placenta) <- annotations


# ##filtering fragments (optional)
# fragment_file_filtered <- "fragments.signacv1.tsv.gz"
# #fragment_file_filtered <- "test_data/atac_v1_placenta_10k_filter.signac.fragments.tsv.gz"
# FilterFragments(
#   fragment.path = fragments,
#   cells = colnames(placenta),  #filter fragments file with keeping barcodes
#   output.path = fragment_file_filtered
# )

# placenta <- SetFragments(
#   object = placenta,
#   file = fragments 
#   #file = fragment_file_filtered
# )

#####


###########QC: TSS score and fragments count (Greenleaf Massive ... paper ), must have id ######
placenta <- NucleosomeSignal(object = placenta)

# compute TSS enrichment score per cell
placenta <- TSSEnrichment(object = placenta, fast = FALSE)
placenta$high.tss <- ifelse(placenta$TSS.enrichment > 2, 'High', 'Low')

options(repr.plot.height=7.5,repr.plot.width=7.5)
TSSPlot(placenta, group.by = 'high.tss') + NoLegend()

placenta$nucleosome_group <- ifelse(placenta$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = placenta, group.by = 'nucleosome_group')


placenta$pct_reads_in_peaks <- placenta$peak_region_fragments / placenta$passed_filters * 100
placenta$blacklist_ratio <- placenta$blacklist_region_fragments / placenta$peak_region_fragments


options(repr.plot.height=15,repr.plot.width=15)
VlnPlot(
  object = placenta,
  features = c('passed_filters','peak_region_fragments','pct_reads_in_peaks', 'nucleosome_signal','TSS.enrichment','blacklist_ratio' ),
  pt.size = 0.1,
  ncol = 3) + NoLegend()

#do filtering

saveRDS(placenta, "placenta.beforefilter.rds")

placenta_bk <- placenta
#191715 x 14441

placenta <- subset(
  x = placenta,
  subset = peak_region_fragments > 3000 &
    peak_region_fragments < 30000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)
placenta
#191715 x 11391

placenta <- RunTFIDF(placenta)
placenta <- FindTopFeatures(placenta, min.cutoff = 'q0')
#191715 mostVar?? too many variable peaks, because of no depth norm?
placenta <- RunSVD(placenta)

options(repr.plot.height=7.5,repr.plot.width=7.5)
DepthCor(placenta)

DefaultAssay(placenta) <- 'peaks'

# placenta <- RunUMAP(object = placenta, reduction = 'lsi', dims = 2:30)
# placenta <- FindNeighbors(object = placenta, reduction = 'lsi', dims = 2:30)
# placenta<- FindClusters(object = placenta, verbose = FALSE, algorithm = 2,resolution=1.2)

# options(repr.plot.height=7.5,repr.plot.width=7.5)
# DimPlot(object = placenta, label = TRUE,cols=color_good,pt.size = 0.8,label.size = 8,reduction = "umap") + NoLegend()


saveRDS(object = placenta,file = "placenta.rds")


##########recalculate with harmony remove batch effect###
V <- placenta@reductions$lsi@cell.embeddings
#11391 x 50 #14441 x 50

meta_data <- placenta@meta.data
#11391 x 26 #14441 x 28

meta_data$sample = ifelse(grepl(pattern = "-1$", x=rownames(meta_data) ), 'D1', 'D2'  )
meta_data$sample <- factor(meta_data$sample,levels=c('D1','D2'))
table(meta_data$sample)
#
#  D1   D2 
#4440 6951 

# D1   D2 
#7416 7025 

V_harmony <- HarmonyMatrix(V,meta_data, 'sample', do_pca = FALSE)
##did not converged 
#11391 x 50 

#converged after 6 iterations
#14441 x 50

placenta@reductions$lsi@cell.embeddings <- V_harmony


DefaultAssay(placenta) <- 'peaks'


options(repr.plot.height=7.5,repr.plot.width=7.5)
DepthCor(placenta)


placenta <- RunUMAP(object = placenta, reduction = 'lsi', dims = 1:50)
placenta <- FindNeighbors(object = placenta, reduction = 'lsi', dims = 1:50)
placenta<- FindClusters(object = placenta, verbose = FALSE, algorithm = 2,resolution=1)


placenta <- RunUMAP(object = placenta, reduction = 'lsi', dims = 2:50)
placenta <- FindNeighbors(object = placenta, reduction = 'lsi', dims = 2:50)
placenta<- FindClusters(object = placenta, verbose = FALSE, algorithm = 1,resolution=0.6)


options(repr.plot.height=7.5,repr.plot.width=7.5)
DimPlot(object = placenta, label = TRUE,pt.size = 0.8,cols=color_good,label.size = 8,reduction = "umap") + NoLegend()


##tuning UMAP
res=0.8
cluster_method=1 #1 ori louvain, 2 louvain with refinement 3 SLM 4 leiden
dim_use = 2:30

for (seed in sample(1:10000,10,replace = FALSE) ){  #RunUMAP default use seed 42
  cat('seed is ',seed,'\n')
  placenta <- FindNeighbors(object = placenta, reduction = 'lsi', dims = dim_use)
  placenta<- FindClusters(object = placenta, verbose = FALSE, algorithm = cluster_method,resolution=0.8,random.seed=seed,method='igraph')
  placenta <- RunUMAP(object = placenta, reduction = 'lsi', dims = dim_use,n.neighbors=15,min.dist =0.1,spread=1.2,seed.use=seed,verbose=FALSE)
  options(repr.plot.height=7.5,repr.plot.width=7.5)
  p = DimPlot(object = placenta, label = TRUE,pt.size = 0.8,cols=color_good,label.size = 8,reduction = "umap") +
    NoLegend() +
    ggtitle(paste("Signac harmony dims=",dim_use, " seed:",seed,' res:',res,' cluster:',cluster_method,sep='') ) +
      theme(
        legend.position = 'none',
        axis.text=element_blank(), 
        axis.title = element_text(size = 15, face = "bold"),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(color="black", fill = NA,size=1),
#         panel.background = element_rect(fill = "white", colour = "white", 
#                 size = rel(1)),
        #panel.border = element_blank(),
        plot.title = element_text(size = 15, face = "bold"),
        #complete = TRUE
        plot.margin = unit(c(1,1,1,1), "lines") #add margin benifit the outside box frame
       )
    print(p)
    
    
}

####################plot with customized method (from work_seurat_harmony.r)##########################

# ##add sample metadata
placenta <- AddMetaData(placenta,metadata = meta_data$sample, col.name = 'sample')


#####get cluster.df
cluster <- Idents(placenta)
umap <- placenta@reductions$umap@cell.embeddings
all.equal(names(cluster),rownames(umap)) #TRUE

cluster.df <- data.frame(cluster=cluster,umap)
metadata <- placenta@meta.data
all.equal(rownames(cluster.df),rownames(metadata))#TRUE

cluster.df.add <- cbind(cluster.df, metadata)
#11391 x 34

cluster.df.add$UMAP_1 = -1 * cluster.df.add$UMAP_1
cluster.df.add$UMAP_2 = -1 * cluster.df.add$UMAP_2

saveRDS(file = 'cluster.df.add.rds',object = cluster.df.add)


###customized way to plot umap-cluster with text halo
centers <- cluster.df.add %>% dplyr::group_by(cluster) %>% dplyr::summarize(x = median(x = UMAP_1), 
        y = median(x = UMAP_2))

centers_shift = data.frame()
##add little shift for text x and y, to plot text halo, borrow from snapATAC
theta= seq(0, 2*pi, length.out=50)
r=0.1
strwidth = 0.5 
strheight = 0.5
xo <- r*strwidth # r*strwidth('A')
yo <- r*strheight #r*strheight('A')
for (i in seq_len(nrow(centers))){
  for (j in theta) {
        centers_shift = rbind(centers_shift,
                              data.frame(
                                  cluster=as.character(unlist(centers[i,'cluster'])),
                                  x=centers[i,'x'] + cos(j)*xo, 
                                  y=centers[i,'y'] + sin(j)*yo
                                 )
                       )
      }
}

####the UMAP plot with annotation
#label right
options(repr.plot.height=5,repr.plot.width=5.5)
ggplot(cluster.df.add,aes(x=UMAP_1,y=UMAP_2,col=cluster  )) +
  geom_point(size = 1.2,show.legend = TRUE,alpha= 1 ) +
  scale_colour_manual(values = color_good)  +
  #scale_colour_gradientn(colors = viridis(6,option = 'C'))  +
  #theme_classic() +
  #theme_bw() +
  
  theme(
        legend.position = 'right',
        axis.text=element_blank(), 
        axis.title = element_text(size = 15, face = "bold"),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(color="black", fill = NA,size=1),
#         panel.background = element_rect(fill = "white", colour = "white", 
#                 size = rel(1)),
        #panel.border = element_blank(),
        plot.title = element_text(size = 15, face = "bold"),
        #complete = TRUE
        plot.margin = unit(c(1,1,1,1), "lines") #add margin benifit the outside box frame
       ) +
 #theme(legend.position = 'none',axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14, face = "bold"),axis.title.y = element_text(size = 14, face = "bold")  ) +
  ggtitle(paste(sample, "two donors, total cells:",nrow(cluster.df.add),  sep=" ") ) +
#   geom_text(data = centers_shift, #the halo
#             mapping = aes(x=x,y=y,label = cluster), 
#             colour = "white", 
#             size = 6.5) +
#   geom_text(data = centers, 
#             mapping = aes(x=x,y=y,label = cluster), 
#             colour = "black", 
#             size = 6) +
  guides(col = guide_legend(override.aes = list(size = 3))) +  ##no effect ??
  #ylim(0,1) + xlim(0,6) +  ##will overwrite scale_x/y_continuous
  #labs(x = "Fragments per barcode", y = "proportion of promoter reads")
  labs(x = "UMAP1", y = "UMAP2")
  #labs(x = "log10(pass_filters_reads)", y = "proportion of enhancer reads")

ggsave(filename = "PLA-8w-ATAC-UMAP.pdf",height=5,width=5.5,useDingbats=FALSE)

##label on cluster
options(repr.plot.height=5,repr.plot.width=5)
ggplot(cluster.df.add,aes(x=UMAP_1,y=UMAP_2,col=cluster  )) +
  geom_point(size = 1.2,show.legend = TRUE,alpha= 1 ) +
  scale_colour_manual(values = color_good)  +
  #scale_colour_gradientn(colors = viridis(6,option = 'C'))  +
  #theme_classic() +
  #theme_bw() +
  theme(
        legend.position = 'none',
        axis.text=element_blank(), 
        axis.title = element_text(size = 15, face = "bold"),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(color="black", fill = NA,size=1),
#         panel.background = element_rect(fill = "white", colour = "white", 
#                 size = rel(1)),
        #panel.border = element_blank(),
        plot.title = element_text(size = 15, face = "bold"),
        #complete = TRUE
        plot.margin = unit(c(1,1,1,1), "lines") #add margin benifit the outside box frame
       ) +
 #theme(legend.position = 'none',axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14, face = "bold"),axis.title.y = element_text(size = 14, face = "bold")  ) +
  ggtitle(paste(sample, "two donors, total cells:",nrow(cluster.df.add),  sep=" ") ) +
  geom_text(data = centers_shift, #the halo
            mapping = aes(x=x,y=y,label = cluster), 
            colour = "white", 
            size = 6.5) +
  geom_text(data = centers, 
            mapping = aes(x=x,y=y,label = cluster), 
            colour = "black", 
            size = 6) +
  ##guides(col = guide_legend(override.aes = list(size = 6))) +  ##no effect ??
  #ylim(0,1) + xlim(0,6) +  ##will overwrite scale_x/y_continuous
  #labs(x = "Fragments per barcode", y = "proportion of promoter reads")
  labs(x = "UMAP1", y = "UMAP2")
  #labs(x = "log10(pass_filters_reads)", y = "proportion of enhancer reads")

ggsave(filename = "PLA-8w-ATAC-UMAP.labelon.pdf",height=5,width=5)
###############


##by donors
options(repr.plot.height=5.5,repr.plot.width=5)
ggplot(cluster.df.add,aes(x=UMAP_1,y=UMAP_2,col= sample  )) +
  geom_point(size = 0.1,show.legend = TRUE,alpha= 0.35 ) +
  scale_colour_manual(values = c('red','navy'))  +
  #scale_colour_gradientn(colors = viridis(6,option = 'C'))  +
  #theme_classic() +
  #theme(legend.position = 'top',axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14, face = "bold"),axis.title.y = element_text(size = 14, face = "bold")  ) +
  theme(legend.position = 'top',
        axis.text=element_blank(), 
        axis.title = element_text(size = 15, face = "bold"),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(color="black", fill = NA,size=1),
        plot.title = element_text(size = 15, face = "bold"),
        #complete = TRUE
        plot.margin = unit(c(1,1,1,1), "lines")
       )+
  ggtitle(paste(sample, "Source of donor ",  sep=" ") ) +
  #ggtitle(paste(sample, "nFeature_RNA",  sep=" ") ) +
  guides(col = guide_legend(override.aes = list(size = 5))) +  ##no effect ??
  #ylim(0,1) + xlim(0,6) +  ##will overwrite scale_x/y_continuous
  #labs(x = "Fragments per barcode", y = "proportion of promoter reads")
  labs(x = "UMAP1", y = "UMAP2")
  #labs(x = "log10(pass_filters_reads)", y = "proportion of enhancer reads")
ggsave(filename = "PLA-8w-ATAC-source-or-donor.pdf",height=5.5,width=5)
####


###by depth
options(repr.plot.height=5.5,repr.plot.width=5)
ggplot(cluster.df.add,aes(x=UMAP_1,y=UMAP_2,col= log10(nCount_RNA) )) +
  geom_point(size = .5,show.legend = TRUE,alpha= 1 ) +
  #scale_colour_manual(values = c('red','navy'))  +
  scale_colour_gradientn(colors = viridis(6,option = 'C'))  +
  #theme_classic() +
  #theme(legend.position = 'top',axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14, face = "bold"),axis.title.y = element_text(size = 14, face = "bold")  ) +
  theme(legend.position = 'top',
        axis.text=element_blank(), 
        axis.title = element_text(size = 15, face = "bold"),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(color="black", fill = NA,size=1),
        plot.title = element_text(size = 15, face = "bold"),
        #complete = TRUE
        plot.margin = unit(c(1,1,1,1), "lines")
       )+
  ggtitle(paste(sample, "Sequence depth",  sep=" ") ) +
  #ggtitle(paste(sample, "nFeature_RNA",  sep=" ") ) +
  #guides(col = guide_legend(override.aes = list(size = 8))) +  ##no effect ??
  #ylim(0,1) + xlim(0,6) +  ##will overwrite scale_x/y_continuous
  #labs(x = "Fragments per barcode", y = "proportion of promoter reads")
  labs(x = "UMAP1", y = "UMAP2")
  #labs(x = "log10(pass_filters_reads)", y = "proportion of enhancer reads")
ggsave(filename = "PLA-8w-ATAC-depth.pdf",height=5.5,width=5)

########################









#########################
####plot dot-heatmap for sequence depth
options(repr.plot.height=7.5,repr.plot.width=7.5)
ggplot(cluster.df.combine ,aes(x=UMAP_1,y=UMAP_2,col=log10(nFeature_peaks))  ) +
  geom_point(size = 1,show.legend = TRUE,alpha= 1 ) +
  scale_colour_gradientn(colors = viridis(6,option = 'C'))  +
  #theme_classic() +
  #theme(legend.position = 'bottom',axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14, face = "bold"),axis.title.y = element_text(size = 14, face = "bold")  ) +
  theme(
        legend.position = 'none',
        axis.text=element_blank(), 
        axis.title = element_text(size = 15, face = "bold"),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(color="black", fill = NA,size=1),
#         panel.background = element_rect(fill = "white", colour = "white", 
#                 size = rel(1)),
        #panel.border = element_blank(),
        plot.title = element_text(size = 15, face = "bold"),
        #complete = TRUE
        plot.margin = unit(c(1,1,1,1), "lines") #add margin benifit the outside box frame
       )+
  ggtitle(paste(sample, "nFeature_peaks",  sep=" ") ) +
  #guides(shape = guide_legend(override.aes = list(size = 8))) +  ##no effect ??
  #ylim(0,1) + xlim(0,6) +  ##will overwrite scale_x/y_continuous
  #labs(x = "Fragments per barcode", y = "proportion of promoter reads")
  labs(x = "UMAP1", y = "UMAP2")
  #labs(x = "log10(pass_filters_reads)", y = "proportion of enhancer reads")


#ggplot(cluster.df.combine ,aes(x=UMAP_1,y=UMAP_2,col=cluster)  ) +
ggplot(cluster.df.combine ,aes(x=UMAP_1,y=UMAP_2,col=sample)  ) +
  geom_point(size = 0.3,show.legend = TRUE,alpha= 0.5 ) +
  #scale_colour_gradientn(colors = viridis(6,option = 'C'))  +
  scale_color_manual(values = c('red','blue')) +
  theme_classic() +
  theme(legend.position = 'bottom',axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14, face = "bold"),axis.title.y = element_text(size = 14, face = "bold")  ) +
  ggtitle(paste(sample, "Samples",  sep=" ") ) +
  #guides(shape = guide_legend(override.aes = list(size = 8))) +  ##no effect ??
  #ylim(0,1) + xlim(0,6) +  ##will overwrite scale_x/y_continuous
  #labs(x = "Fragments per barcode", y = "proportion of promoter reads")
  labs(x = "UMAP_1", y = "UMAP_2") +
  guides(col = guide_legend(override.aes = list(size = 8)))
  #labs(x = "log10(pass_filters_reads)", y = "proportion of enhancer reads")

ggplot(cluster.df.combine ,aes(x=UMAP_1,y=UMAP_2,col=TSS.enrichment)  ) +
  geom_point(size = 0.3,show.legend = TRUE,alpha= 0.5 ) +
  scale_colour_gradientn(colors = viridis(6,option = 'C'))  +
  #scale_color_manual(values = c('red','blue')) +
  theme_classic() +
  theme(legend.position = 'bottom',axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14, face = "bold"),axis.title.y = element_text(size = 14, face = "bold")  ) +
  ggtitle(paste(sample, "TSS.enrichment",  sep=" ") ) +
  #guides(shape = guide_legend(override.aes = list(size = 8))) +  ##no effect ??
  #ylim(0,1) + xlim(0,6) +  ##will overwrite scale_x/y_continuous
  #labs(x = "Fragments per barcode", y = "proportion of promoter reads")
  labs(x = "UMAP_1", y = "UMAP_2") +
  guides(col = guide_legend(override.aes = list(size = 8)))
  #labs(x = "log10(pass_filters_reads)", y = "proportion of enhancer reads")



################








####save rds###
#saveRDS(placenta, "placenta.harmony_corrected.snATAC.rds")
placenta = readRDS("placenta.harmony_corrected.snATAC.rds")




##calculate gene activity for all genes
gene.activities <- GeneActivity(placenta)

# add the gene activity matrix to the Seurat object as a new assay and normalize it
placenta[['RNA']] <- CreateAssayObject(counts = gene.activities)
placenta <- NormalizeData(
  object = placenta,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(placenta$nCount_RNA)
)



DefaultAssay(placenta) <- 'RNA'

options(repr.plot.height=12,repr.plot.width=15)
FeaturePlot(
  object = placenta,
  features = c('CGA', 'FLT1', 'PSG1', 'CDH1', 'DNMT1', 'ERVFRD-1','HLA-G','MMP2','VIM'),
  pt.size = 0.5,
  max.cutoff = 'q95',
  ncol = 3
)

options(repr.plot.height=12,repr.plot.width=15)
FeaturePlot(
  object = placenta,
  features = c('MKI67', 'PCNA', 'PECAM1', 'GCM1', 'PAGE4', 'LAIR2','MMP11','KRT7','HLA-A'),
  pt.size = 0.5,
  max.cutoff = 'q95',
  ncol = 3
)

options(repr.plot.height=9,repr.plot.width=15)
FeaturePlot(
  object = placenta,
  features = c('CSH2','PAPPA', 'FLT1', 'LEP', 'PSG1', 'PSG8'),
  pt.size = 0.5,
  max.cutoff = 'q95',
  ncol = 3
)







#####filtering noisy clusters#####






saveRDS(placenta, "placenta.rds")

