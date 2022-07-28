
library(SnapATAC) #v1.0.0
##use R CMD INSTALL /home/mjwang/pwdex/placenta_10X_combine/02.snapATAC_harmony/SnapATAC-v1.0.release2.modify.tar.gz, comment out import of plot3D scatter2D in NAMESPACE, which will call a x11 connection

library('magrittr')

library('dplyr')
library(GenomicRanges)
library(ggplot2)

library(hrbrthemes)

library('patchwork')
library(viridis)


sample <- 'PLA-8w-ATAC'
#sample <- 'placenta at 8 weeks stage '






#########customized colors########
color_snap = c('1'='grey','2'='#E31A1C','3'='#FFD700','4'='#771122','5'='#777711','6'='#1F78B4','7'='#68228B','8'='#AAAA44','9'='#60CC52','10'='#771155','11'='#DDDD77','12'='#774411','13'='#AA7744','14'='#AA4455','15'='#117744')
#names(color_snap) <- NULL

#modified for CTB with dark red colors
color_snap_mod1 = c('1'='#777711','2'='#E31A1C','3'='#68228B','4'='#771122','5'='grey','6'='#1F78B4','7'='#FFD700','8'='#AAAA44','9'='#60CC52','10'='#771155','11'='#DDDD77','12'='#774411','13'='#AA7744','14'='#AA4455','15'='#117744')

color_gradient_my <- c(
    rgb(5,48,97,maxColorValue = 255),
    rgb(42,113,178,maxColorValue = 255),
    rgb(147,198,222,maxColorValue = 255),
    rgb(239,243,245,maxColorValue = 255),
    rgb(253,219,199,maxColorValue = 255),
    rgb(214,96,77,maxColorValue = 255),
    rgb(121,5,34,maxColorValue = 255)

)

color_cellranger <-c('#820610','#C50F1E','#F42428','#F86D30','#FBB33D','#FCFB4E','#C0FB61','#87FB8F','#41FEF9','#2BAED7','#155CB1','#08238D') #for depth

#############ArchR colors############
#hmcols <- colorRamps::blue2green2red(length(bks) ) #colors
color_peak <- ArchR::paletteContinuous(set = 'solarExtra',n=256,reverse=FALSE)  
color_tfdev = ArchR::paletteContinuous(set = 'blueYellow',n=257,reverse=FALSE)                       
#color_ga <- paletteContinuous(set='solarExtra',n=257,reverse=FALSE) 
#color_ga <- paletteContinuous(set='horizon',n=257,reverse=FALSE)                      
#color_ga <- paletteContinuous(set='horizonExtra',n=257,reverse=FALSE)  #good        
#color_ga <- paletteContinuous(set='greenBlue',n=257,reverse=FALSE)
#color_ga <- paletteContinuous(set='blueYellow',n=257,reverse=FALSE)
color_ga <- ArchR::paletteContinuous(set='greyMagma',n=257,reverse=FALSE)


####read in color palette list###
color_set_yellowbrick.flat <- readRDS('/home/mjwang/progs/misc-tools/r/color_set_yellowbrick.flat.rds')
color_list_archr <- readRDS('/home/mjwang/progs/misc-tools/r/ArchR.color_list.rds')


options(repr.plot.width=12,repr.plot.height=6)
par(mfrow=c(3,3),oma=c(1,1,1,1))
for(name in names(color_list_archr)){ 
  len = length(color_list_archr[[name]])
  color = color_list_archr[[name]]
  #barplot(rep(5,len),col = color,main = name,cex.main=2)
  barplot(1:len,col = color,main = name,cex.main=2)
}

options(repr.plot.width=12,repr.plot.height=6)
par(mfrow=c(3,3),oma=c(1,1,1,1))
for(name in names(color_set_yellowbrick.flat)){ 
  len = length(color_set_yellowbrick.flat[[name]])
  color = color_set_yellowbrick.flat[[name]]
  #barplot(rep(5,len),col = color,main = name,cex.main=2)
  barplot(1:len,col = color,main = name,cex.main=2)
}



##########reload#########################
x.after.sp <- readRDS("x.after.sp.onestep_filter.final.new.rds")
#x.after.sp.bk <- readRDS("x.after.sp.onestep_filter.final.bk.rds")
#all.equal(x.after.sp,x.after.sp.bk) #TRUE
x.after.sp 
number of barcodes: 12526
number of bins: 620094
number of genes: 58344
number of peaks: 299708
number of motifs: 0

cluster.df.add <- readRDS("cluster.df.add.final.rds")

##check cluster in rds and x.after.sp
cluster.df = data.frame(cluster=x.after.sp@cluster,UMAP_1=x.after.sp@umap[,1],UMAP_2=x.after.sp@umap[,2]) 
rownames(cluster.df) = rownames(x.after.sp@metaData)

#all.equal (rownames(cluster.df),rownames(x.after.sp@metaData) )#TRUE
cluster.df.add <- cbind(cluster.df,x.after.sp@metaData)
all.equal(cluster.df.add , readRDS("cluster.df.add.final.rds")) #TRUE



##quick look the ori umap embedding and cluster
options(repr.plot.height=7.5,repr.plot.width=7.5)
plotViz( #with text label halo, point size = 1
    obj=x.after.sp,
    method="umap", 
    main=paste("snapATAC harmony",sep=''),
    point.color=x.after.sp@cluster, 
    point.size=0.3, 
    point.shape=19, 
    point.alpha=0.8, 
    text.add=TRUE,
    text.size=1,
    text.color="black",
    text.halo.add=TRUE,
    text.halo.color="white",
    text.halo.width=0.2,
    #down.sample=10000,
    legend.add=FALSE
    );

##plot a single DAR peak in all cell (too weak for one peak)
#gr.all <- x.after.sp@peak 

# ##the igv peak
# peak_region = "chr7:87414143-87414728"

# #find the true peak in pmat
# peak_gr <- GRanges(seqnames = 'chr7',ranges = IRanges(start = 87414143,end = 87414728 )  )
# hit <- findOverlaps(query = peak_gr,subject = gr.all)

# hit.df <- as.data.frame(hit)
# peak_idx <- hit.df[1,'subjectHits']
# peak_hit <- mcols(gr.all[hit.df[1,'subjectHits'] ])$name

# # plotFeatureSingle( #use plot3D, can not append ggplot theme
# #     #plotFeatureSingle_new(
# #         obj=x.after.sp,
# #         feature.value=x.after.sp@pmat[,peak_idx],
# #         method="umap", 
# #         main=peak_hit,
# #         point.size=0.2, 
# #         point.shape=19 
# #         #down.sample=10000,
# #         #quantiles=c(0.1, 1)##make sharper
# #       )


##plot the DAG gene activity score (simple gene body sum, not cicero )
##gene = 'CROT'
gene = 'DNMT1'
gene = 'PSG8'
gene = 'HLA-G'
gene = 'VIM'
plotFeatureSingle( #use plot3D, can not append ggplot theme
    #plotFeatureSingle_new(
        obj=x.after.sp,
        feature.value=x.after.sp@gmat[,gene],
        method="umap", 
        main=gene,
        point.size=0.2, 
        point.shape=19, 
        #down.sample=10000,
        quantiles=c(0.1, 1),##make sharper
        #pdf.file.name = paste("pdfs/knowMarkers/know_marker.",gene,".UMAP.pdf",sep=''),
        pdf.width=7, 
	    pdf.height=7
      )


###filter STR and cluster 7

cluster.df.add <- subset(cluster.df.add, ! cluster %in% c('12','10','11','13','14','15') ) #8140
table(cluster.df.add$cluster)

  1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
1832 1783 1090 1546 1475 1509 1019  444  395    0    0    0    0    0    0 

#####################the final filtering object###################
idx <- which (rownames(x.after.sp@metaData) %in% rownames(cluster.df.add)     )
table(rownames(x.after.sp@metaData) %in% rownames(cluster.df.add))

FALSE  TRUE 
 1433 11093 



#x.after.sp.bk <- x.after.sp
x.after.sp <- x.after.sp[idx,]

number of barcodes: 11093
number of bins: 620094
number of genes: 58344
number of peaks: 299708
number of motifs: 0

all.equal(rownames(x.after.sp@metaData),rownames(cluster.df.add) ) #TRUE
all.equal(rownames(x.after.sp@bmat),rownames(cluster.df.add) ) #diff

#drop levels
x.after.sp@cluster <- droplevels(x.after.sp@cluster)
table(x.after.sp@cluster)
  1    2    3    4    5    6    7    8    9 
1832 1783 1090 1546 1475 1509 1019  444  395 

#rename cluster id
mapid <- list('1'='1','2'='2', '3'='3','4'='4','5'='5',
'6'='6',
'7'='7',
'8'='8',
'9'='9'
 )

names(mapid)

cluster <- as.character(x.after.sp@cluster)

for(i in 1:length(cluster)){ cluster[i] = mapid[[ cluster[i] ]] }

cluster <- factor(cluster,levels=gtools::mixedsort( names(table(cluster)) ))

table(cluster)

 1    2    3    4    5    6    7    8    9 
1832 1783 1090 1546 1475 1509 1019  444  395 

#   1    2    3    4    5    6    7    8    9   10   11   12 
# 1801 1749 1539 1044 1048  613  441  346  375  194  127   59 


x.after.sp@cluster <- cluster #the same, no need to rewrite

##quick look the ori umap embedding and cluster
options(repr.plot.height=7.5,repr.plot.width=7.5)
plotViz( #with text label halo, point size = 1
    obj=x.after.sp,
    method="umap", 
    main=paste("snapATAC harmony",sep=''),
    point.color=x.after.sp@cluster, 
    point.size=0.3, 
    point.shape=19, 
    point.alpha=0.8, 
    text.add=TRUE,
    text.size=1,
    text.color="black",
    text.halo.add=TRUE,
    text.halo.color="white",
    text.halo.width=0.2,
    #down.sample=10000,
    legend.add=FALSE
    );



##add cell name and cell color


map_cellname <- list(
    '1'='STB2',
    '2'='STB3', 
    '3'='CTB2',
    '4'='STB4',
    '5'='STB naive',
    '6'='CTB1',
    '7'='STB1',
    '8'='STB5',
    '9'='Fusion-competent'
 )



#STB ()  navy gradient
#Syncitial knot: 
#CTB 
#STR greens (3)

reds <- c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#ffffbf')
blues <- c('#e0f3f8','#abd9e9','#74add1','#4575b4','#313695')

map_cellcolor <- list(
    '1'=reds[3],#'STB-2',
    '2'=reds[2],#'STB-3', 
    '3'=blues[4],#'CTB-2',
    '4'=reds[1],#'STB-4',
    '5'=reds[6],#'naive STB',
    '6'=blues[3],#'CTB-1',
    '7'=reds[4],#'STB-1',
    '8'=reds[5],#'STB-5',
    '9'='darkgreen'#'Fusion Component'
    
#     '1'=reds[5], #STB1
#     '2'=reds[6],  #STB-naive
#     '3'=reds[3],  #STB3
#     '4'=reds[2], #STB4
#     '5'=reds[4], #STB2
#     '6'=reds[1], #Syncytial knot
#     '7'='#694d9f', #STB-new
#     '8'='darkgreen' #CTB
 )


map_cellcolor <- unlist(map_cellcolor[c('6','3','9','5','7','1','8','2','4')])

saveRDS(map_cellcolor,'map_cellcolor.rds')

#saveRDS(map_cellcolor[c(6,3,9,5,7,1,2,4)],'map_cellcolor.rds')

#saveRDS(map_cellcolor[c(6,3,9,5,7,1,2,4,8)],'map_cellcolor.rds')

#saveRDS(map_cellcolor[c('6','3','9','5','7','1','8','2','4')],'map_cellcolor.rds')

write.table(unlist(map_cellcolor[c(6,3,9,5,7,1,2,4)]),'map_cellcolor.txt',sep=',',row.names = TRUE, col.names = FALSE)


##get tne final cluster.df.add
cluster.df = data.frame(cluster=x.after.sp@cluster,UMAP_1=x.after.sp@umap[,1],UMAP_2=x.after.sp@umap[,2]) 
rownames(cluster.df) = rownames(x.after.sp@metaData)

all.equal (rownames(cluster.df),rownames(x.after.sp@metaData) )#TRUE
cluster.df.add <- cbind(cluster.df,x.after.sp@metaData)


all.equal(x.after.sp@cluster,cluster.df.add$cluster)#TRUE

cellname <- as.character(x.after.sp@cluster)
cellcolor<- as.character(x.after.sp@cluster)
for(i in 1:length(cellname)){ cellname[i] = map_cellname[[ cellname[i] ]] }
for(i in 1:length(cellcolor)){ cellcolor[i] = map_cellcolor[[ cellcolor[i] ]] }

cluster.df.add[,'cellname'] <- factor(cellname)
cluster.df.add[,'cellcolor'] <- factor(cellcolor)



#########################save/reload the final object #################

#saveRDS(x.after.sp,"x.after.sp.final.final.rds")     
x.after.sp <- readRDS("x.after.sp.final.final.rds")

#write.table(cluster.df,file='snapATAC.PLA-8w-ATAC.umap.bin5k.final.txt',sep='\t',quote = FALSE, row.names = TRUE, col.names = TRUE)

#saveRDS(cluster.df.add,"cluster.df.add.final.final.rds")

cluster.df.add <- readRDS("cluster.df.add.final.final.rds")



##how many CTB in 8w atac data?
subset(cluster.df.add,subset = cluster %in% c('6','3','9') )
#2994 (compatible to rna of 2878??)




######plot the final fixed cluster, color and label with customized plot function
# centers <- cluster.df.add %>% dplyr::group_by(cellname) %>% dplyr::summarize(x = median(x = UMAP_1), 
#         y = median(x = UMAP_2))

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
                                  ##cluster=as.character(unlist(centers[i,'cellname'])),
                                  cluster=as.character(unlist(centers[i,'cluster'])),
                                  x=centers[i,'x'] + cos(j)*xo, 
                                  y=centers[i,'y'] + sin(j)*yo
                                 )
                       )
      }
}

####the UMAP plot with annotation


##label on cluster
options(repr.plot.height=5,repr.plot.width=5.5,repr.plot.res = 150)
ggplot(cluster.df.add,aes(x=UMAP_1,y=UMAP_2,col=cluster  )) +
#ggplot(cluster.df.add,aes(x=UMAP_1,y=UMAP_2,col=cellcolor  )) +
  geom_point(size = 0.2,show.legend = TRUE,alpha= 1 ) +
  #scale_colour_manual(values = color_good)  +
  scale_colour_manual(values = unlist(map_cellcolor) )  +
  #scale_colour_manual(values = color_snap_mod1)  +
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
            ##size = 4.5) +
            size = 6.5) +
  geom_text(data = centers, 
            ##mapping = aes(x=x,y=y,label = cellname), 
            mapping = aes(x=x,y=y,label = cluster), 
            colour = "black", 
            ##size = 4.5) +
            size = 6.5) +
  ##guides(col = guide_legend(override.aes = list(size = 6))) +  ##no effect ??
  #ylim(0,1) + xlim(0,6) +  ##will overwrite scale_x/y_continuous
  xlim(min(cluster.df.add$UMAP_1 )*1.5,1.5*max(cluster.df.add$UMAP_1 )) +
  #labs(x = "Fragments per barcode", y = "proportion of promoter reads")
  labs(x = "UMAP1", y = "UMAP2")
  #labs(x = "log10(pass_filters_reads)", y = "proportion of enhancer reads")

##ggsave(filename = "pdfs/PLA-8w-ATAC-UMAP.labelon.pdf",height=5,width=5.5,useDingbats=FALSE)
ggsave(filename = "pdfs/PLA-8w-ATAC-UMAP.labelon.clusterid.pdf",height=5,width=5.5,useDingbats=FALSE)
##ggsave(filename = "pdfs/PLA-8w-ATAC-UMAP.labeloff.pdf",height=5,width=5.5,useDingbats=FALSE)
###############



##by donors
options(repr.plot.height=5.5,repr.plot.width=5)
ggplot(cluster.df.add,aes(x=UMAP_1,y=UMAP_2,col= sample  )) +
  geom_point(size = .1,show.legend = TRUE,alpha= 1 ) +
  scale_colour_manual(values = c('#94C6DD', '#1273AE'))  + #'#94C6DD', '#1273AE'; 'red','navy'
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
  xlim(min(cluster.df.add$UMAP_1 )*1.5,1.5*max(cluster.df.add$UMAP_1 )) +
  #labs(x = "Fragments per barcode", y = "proportion of promoter reads")
  labs(x = "UMAP1", y = "UMAP2")
  #labs(x = "log10(pass_filters_reads)", y = "proportion of enhancer reads")

ggsave(filename = "pdfs/PLA-8w-ATAC-source-of-donor.pdf",height=5.5,width=5,useDingbats=FALSE)




###by depth
#options(repr.plot.height=5.5,repr.plot.width=5)
options(repr.plot.height=5.5,repr.plot.width=5)
#ggplot(cluster.df.add,aes(x=UMAP_1,y=UMAP_2,col= passed_filters)) +
#ggplot(cluster.df.add,aes(x=UMAP_1,y=UMAP_2,col= log(UQ) )) +
ggplot(cluster.df.add,aes(x=UMAP_1,y=UMAP_2,col= logUMI )) +
  geom_point(size = .2,show.legend = TRUE,alpha= 1 ) +
  #scale_colour_manual(values = c('red','navy'))  +
  scale_colour_gradientn(colors = rev(color_cellranger))  +
  #scale_colour_gradientn(colors = viridis(6,option = 'D'))  + 
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
  ggtitle(paste(sample, "Sequence depth (log)",  sep=" ") ) +
  #ggtitle(paste(sample, "Sequence depth ",  sep=" ") ) +
  #ggtitle(paste(sample, "nFeature_RNA",  sep=" ") ) +
  #guides(col = guide_legend(override.aes = list(size = 8))) +  ##no effect ??
  #ylim(0,1) + xlim(0,6) +  ##will overwrite scale_x/y_continuous
  xlim(min(cluster.df.add$UMAP_1 )*1.5,1.5*max(cluster.df.add$UMAP_1 )) +
  #labs(x = "Fragments per barcode", y = "proportion of promoter reads")
  labs(x = "UMAP1", y = "UMAP2")
  #labs(x = "log10(pass_filters_reads)", y = "proportion of enhancer reads")
ggsave(filename = "pdfs/PLA-8w-ATAC-depth.pdf",height=5.5,width=5,useDingbats=FALSE)


##by FRiP (promoter percentage)
options(repr.plot.height=5.5,repr.plot.width=5)
ggplot(cluster.df.add,aes(x=UMAP_1,y=UMAP_2,col= promoter_ratio )) +
  geom_point(size = .2,show.legend = TRUE,alpha= 1 ) +
  #scale_colour_manual(values = c('red','navy'))  +
  scale_colour_gradientn(colors = rev(color_cellranger))  +
  #scale_colour_gradientn(colors = viridis(6,option = 'D'))  +
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
  ggtitle(paste(sample, "Promoter ratio",  sep=" ") ) +
  #ggtitle(paste(sample, "nFeature_RNA",  sep=" ") ) +
  #guides(col = guide_legend(override.aes = list(size = 8))) +  ##no effect ??
  #ylim(0,1) + xlim(0,6) +  ##will overwrite scale_x/y_continuous
  xlim(min(cluster.df.add$UMAP_1 )*1.5,1.5*max(cluster.df.add$UMAP_1 )) +
  #labs(x = "Fragments per barcode", y = "proportion of promoter reads")
  labs(x = "UMAP1", y = "UMAP2")
  #labs(x = "log10(pass_filters_reads)", y = "proportion of enhancer reads")
ggsave(filename = "pdfs/PLA-8w-ATAC-promotorRatio.pdf",height=5.5,width=5,useDingbats=FALSE)


########################





## Heretical clustering ###
# calculate the ensemble signals for each cluster
ensemble.ls = lapply(split(seq(length(x.after.sp@cluster)), x.after.sp@cluster), function(x){
	SnapATAC::colMeans(x.after.sp[x,], mat="bmat");
	})
# cluster using 1-cor as distance  
hc = hclust(as.dist(1 - cor(t(do.call(rbind, ensemble.ls)))), method="ward.D2");

pdf( "pdfs/hclust.cluster.pdf",height=5.5,width=4,useDingbats = FALSE)
options(repr.plot.height=5.5,repr.plot.width=4)
plot(hc, hang=-1, xlab="");
dev.off()


####bin ave coverage(depth) of cluster (bmat)
ensemble.ls.df <- do.call(cbind,ensemble.ls) #620094  
#all.equal(t(do.call(rbind, ensemble.ls)),do.call(cbind,ensemble.ls)) #TRUE
ensemble.ls.df <- ensemble.ls.df[rowSums(ensemble.ls.df) != 0,] #573453 #568573 #575002
ensemble.ls.df

options(repr.plot.height=5,repr.plot.width=5)
boxplot(ensemble.ls.df,outline = FALSE,col = color_snap_mod1,las=2)


##bin sum coverage(depth) of cluster (bmat)
ensembleSum.ls = lapply(split(seq(length(x.after.sp@cluster)), x.after.sp@cluster), function(x){
	SnapATAC::colSums(x.after.sp[x,], mat="bmat");
	})
ensembleSum.ls.df <- do.call(cbind,ensembleSum.ls) #620094

ensembleSum.ls.df <- ensembleSum.ls.df[rowSums(ensembleSum.ls.df) != 0,] #568573 #574062
ensembleSum.ls.df

options(repr.plot.height=5,repr.plot.width=5)
boxplot(ensembleSum.ls.df,outline = FALSE,col = color_snap_mod1,las=2) #quick

##bin sum coverage(depth) of cluster (bmat) by donor
#all.equal(cluster.df.add[,-c(1,2,3)],x.after.sp@metaData) #TRUE
ensembleSum.ls = lapply(split(seq(length(x.after.sp@cluster)),  #group by donor and cluster
                              list(cluster=x.after.sp@cluster,
                                   lib=x.after.sp@metaData$library)
                             ), 
                        function(x){SnapATAC::colSums(x.after.sp[x,], mat="bmat")}
                       )
ensembleSum.ls.df <- do.call(cbind,ensembleSum.ls) #620094
saveRDS(ensembleSum.ls.df,'ensembleSum.ls.df.rds')


ensembleSum.ls.df.d1 <- ensembleSum.ls.df[,grep(pattern = 'D1',x=colnames(ensembleSum.ls.df))]
ensembleSum.ls.df.d2 <- ensembleSum.ls.df[,grep(pattern = 'D2',x=colnames(ensembleSum.ls.df))]



ensembleSum.ls.df.d1 <- ensembleSum.ls.df.d1[rowSums(ensembleSum.ls.df.d1) != 0,] #566361 #572203
ensembleSum.ls.df.d1
ensembleSum.ls.df.d2 <- ensembleSum.ls.df.d2[rowSums(ensembleSum.ls.df.d2) != 0,] #565418 #568099
ensembleSum.ls.df.d2


pdf( "pdfs/boxplot_depth_by_donor.pdf",height=7.5,width=10,useDingbats = FALSE)
par(mfrow = c(2,1),mar=c(3,3,3,1))
#options(repr.plot.height=10,repr.plot.width=5)
boxplot(ensembleSum.ls.df.d1,outline = FALSE,col = color_snap_mod1,las=2,ylim=c(0,40),main='donor1') #quick
boxplot(ensembleSum.ls.df.d2,outline = FALSE,col = color_snap_mod1,las=2,ylim=c(0,40),main='donor2') #quick
dev.off()



######stat d1 d2 and cluster cell number correlation######
res.stat <- table(cluster.df.add$cluster,cluster.df.add$library) 
   
      D1   D2
  1  779 1053
  2  774 1009
  3  450  640
  4  667  879
  5  631  844
  6  585  924
  7  347  672
  8  172  272
  9  152  243

#      D1   D2
#   1 1068  733
#   2  925  824
#   3  992  547
#   4  659  385
#   5  597  451
#   6  321  292
#   7  181  165
#   8  170   24

#        D1   D2
#   1  1068  733
#   2   925  824
#   3   992  547
#   4   659  385≠
#   5   597  451
#   6   321  292
#   7   243  198
#   8   181  165
#   9   295   80
#   10  170   24
#   11   81   46
#   12   36   23

as.data.frame(res.stat) %>% dplyr::group_by(Var2) %>% summarise(sum = sum(Freq))

  D1   D2 
4557 6536




res.cor <- cor(res.stat) #0.924 #0.977 #0.983 #0.988
res.stat[,1] = -1 * res.stat[,1]

res.stat.df <- as.data.frame(res.stat)
#res.stat.df <- reshape2::melt(res.stat)
colnames(res.stat.df) <- c('cluster','sample','count')

##plot horizonal barplot  #  '#94C6DD', '#1273AE'  vs '#F3A585','#C80927'
options(repr.plot.height=5.5,repr.plot.width=4)
ggplot(res.stat.df, aes(fill=sample, y=count, x=cluster )) + 
    #geom_hline(yintercept = c(0,25,50,75),linetype='solid',size=.3,col='black') +
    #geom_hline(yintercept = c(0),linetype='solid',size=1,col='black') +
    #geom_bar(position="stack", stat="identity",alpha=1,width = 0.5) +  
    geom_bar(position=position_stack(reverse=TRUE), stat="identity",width = 0.5) +  
    #xlim(100,0) +
    #scale_x_continuous(breaks = seq(100,0,-25),labels=seq(100,0,-25)) +
    #scale_y_reverse() +
    #annotate('text',x = 16, y = 1000,label = paste('spearman cor=',cor(res.stat)[1,2]) ) +
    coord_flip() +
    #scale_fill_viridis(discrete = T,option = "E") +
    scale_fill_manual(values = c('#94C6DD', '#1273AE'),labels=c('D1','D2'),name='sample' ) +
    #ggtitle("cell number count") +
    labs(title = "cell number count", subtitle=paste('spearman cor=',round(res.cor[1,2],digits=2)) ) +
    theme_ipsum(base_family = 'sans') + #to avoid Arial narrow problem in ggsave
    ylab("count")

ggsave(filename="pdfs/donor1.vs.donor2.cluster.cor.pdf",width = 4, height = 5.5)





#####gene activity score annatation with know genes



marker.genes = c(
    "DNMT1", "CDH1", "MKI67",
    "FLT1", "CSHL1", "PSG8", 
    "ERVFRD-1", "LAIR2", "PLAC8",
    'VIM','PECAM1','CD14'
  );


marker.genes = c( #try to distinguish STB terminals 
    "FLT1", 'LEP','INSIG2', #FLT1-enriched early stage
    "CSHL1",'CSH1','PAPPA',  #PAPPA-enriched late stage
    'PSG1','CGA','GCM1' #general

  );

##general STB: PSG1,PSG8,CGA,;naive STB: ERVFRD-1,GCM1; hormone-bias STB:PAPPA,CSHL1;FLT1, ENG,LEP, INSIG2
#CTB DNMT1, CDH1
#Syncytial knot, CDKN1A SPATA18, CROT, PTCHD4
#STB-new:  DDX60,LIFR,MAPKAPK3,VAC14
#total 9 or 12 subplot
marker.genes = c(
    'PSG1','PSG8','CGA','CGB3',
     'PAPPA','CSHL1', 'CSH1','CSH2',
     'FLT1', 'ENG','LEP','INSIG2',
    'ERVFRD-1','GCM1','DNMT1', 'CDH1',
     'TEAD3','TEAD4','TEAD1','GATA3',
    'CDKN1A', 'SPATA18', 'CROT', 'PTCHD4',
    'DDX60','LIFR','MAPKAPK3','VAC14'
)


#marker.genes <- c('PAPPA', 'FLT1', 'LEP', 'PSG8', 'PTCHD4', 'SH3TC2')
marker.genes <- c('PAPPA', 'FLT1', 'LEP', 'PSG8', 'ERVFRD-1', 'SH3TC2')


marker.genes <- c('DNMT1','ERVFRD-1', 'PSG8', 'CGA', 'LEP','CSHL1' , 'PAPPA','FLT1')

marker.genes %in% colnames(x.after.sp@gmat)



#######recheck gmat pmat and save renamed pmat####
gmat <- readRDS('gmat.rds') #12526 x 58344 (full)
pmat <- readRDS('pmat.raw.rds') #12526 x 299708 (full)

##check core barcode
check.df <- data.frame(colid = rownames(gmat), cellid = rownames(pmat)   )
check.df$colid <- gsub(pattern = "^placenta_donor1#|placenta_donor2#",replacement = "", x = check.df$colid )
check.df$colid <- gsub(pattern = '-1|-2$',replacement = "", x = check.df$colid )
check.df$cellid <- gsub(pattern = '-1|-2$',replacement = "", x = check.df$cellid )
all.equal(check.df$colid,check.df$cellid) #TRUE, the core barcode 

rownames(pmat) <- rownames(gmat)#12526 x 299708 (full)


all.equal(x.after.sp@gmat, gmat[rownames(x.after.sp@metaData),] ,check.attributes = FALSE) #TRUE
all.equal(x.after.sp@pmat, pmat[rownames(x.after.sp@metaData),] ,check.attributes = FALSE) #TRUE

max(pmat) #255

saveRDS(pmat,'pmat.raw.rds') #save cellid changed pmat

#####



######customized plot marker gene scatter plot#########

marker.gmat <- as.data.frame(as(x.after.sp@gmat[, marker.genes],'matrix'))
#marker.gmat <- as.data.frame(as(x.after.sp@gmat[, marker.genes],'matrix'))
rownames(marker.gmat) <- rownames(x.after.sp@metaData)
all.equal (rownames(cluster.df.add),rownames(marker.gmat) ) #TRUE
marker.gmat.df <- cbind(cluster.df.add[,c('cluster','UMAP_1','UMAP_2')],marker.gmat )


##cutoff quantiles
#low.q <- 0.1
#high.q <- 0.9

# low.q <- c('DNMT1'=0.1,'CDH1'=0.1,'PAGE4'=0.1,'FLT1'=0,'CSHL1'=0.1,'PSG8'=0,'ERVFRD-1'=0,'LAIR2'=0,'PLAC8'=0.1,'VIM'=0,'PECAM1'=0,'MKI67'=0,'CD14'=0)
# high.q <- c('DNMT1'=1,'CDH1'=1,'PAGE4'=1,'FLT1'=1,'CSHL1'=1,'PSG8'=1,'ERVFRD-1'=1,'LAIR2'=1,'PLAC8'=1,'VIM'=1,'PECAM1'=1,'MKI67'=1,'CD14'=1)

# low.q <- c('FLT1'=0,'LEP'=0,'INSIG2'=0,'CSHL1'=0,'CSH1'=0,'PAPPA'=0,'PSG1'=0,'CGA'=0,'GCM1'=0)
# high.q <- c('FLT1'=1,'LEP'=1,'INSIG2'=1,'CSHL1'=1,'CSH1'=1,'PAPPA'=1,'PSG1'=1,'CGA'=1,'GCM1'=1)


##used for viridis D 256 color
low.q <- c('PSG1'=0,'PSG8'=0,'CGA'=0,'CGB3'=0,'PAPPA'=0.1,'CSHL1'=0.1,'CSH1'=0.1,'CSH2'=0.1,'FLT1'=0.1,'ENG'=0.1,'LEP'=0.1,'INSIG2'=0.05,'ERVFRD-1'=0.1,'GCM1'=0.05,'DNMT1'=0,'CDH1'=0,'TEAD3'=0.1,'TEAD4'=0,'TEAD1'=0,'GATA3'=0.1,'CDKN1A'=0.1,'SPATA18'=0,'CROT'=0,'PTCHD4'=0,'DDX60'=0.1,'LIFR'=0.05,'MAPKAPK3'=0.1,'VAC14'=0.1)
high.q <- c('PSG1'=1,'PSG8'=1,'CGA'=1,'CGB3'=1,'PAPPA'=1,'CSHL1'=1,'CSH1'=1,'CSH2'=1,'FLT1'=1,'ENG'=1,'LEP'=1,'INSIG2'=1,'ERVFRD-1'=1,'GCM1'=1,'DNMT1'=1,'CDH1'=1,'TEAD3'=1,'TEAD4'=1,'TEAD1'=1,'GATA3'=1,'CDKN1A'=1,'SPATA18'=1,'CROT'=1,'PTCHD4'=1,'DDX60'=1,'LIFR'=1,'MAPKAPK3'=1,'VAC14'=1)


# ##used for color_gradient_my
# low.q <- c('PSG1'=0,'PSG8'=0,'CGA'=0,'CGB3'=0,'PAPPA'=0.1,'CSHL1'=0.1,'CSH1'=0.1,'CSH2'=0.1,'FLT1'=0.1,'ENG'=0.1,'LEP'=0.1,'INSIG2'=0.05,'ERVFRD-1'=0.1,'GCM1'=0.05,'DNMT1'=0,'CDH1'=0,'TEAD3'=0.1,'TEAD4'=0,'TEAD1'=0,'GATA3'=0.1,'CDKN1A'=0.1,'SPATA18'=0,'CROT'=0,'PTCHD4'=0,'DDX60'=0.06,'LIFR'=0.02,'MAPKAPK3'=0.02,'VAC14'=0.08)
# high.q <- c('PSG1'=1,'PSG8'=1,'CGA'=1,'CGB3'=1,'PAPPA'=1,'CSHL1'=1,'CSH1'=1,'CSH2'=1,'FLT1'=1,'ENG'=1,'LEP'=1,'INSIG2'=1,'ERVFRD-1'=1,'GCM1'=1,'DNMT1'=1,'CDH1'=1,'TEAD3'=1,'TEAD4'=1,'TEAD1'=1,'GATA3'=1,'CDKN1A'=1,'SPATA18'=1,'CROT'=1,'PTCHD4'=1,'DDX60'=1,'LIFR'=1,'MAPKAPK3'=1,'VAC14'=1)


##for gmat color_tfdev

#low.q <- c('PAPPA'=0, 'FLT1'=0.2, 'LEP'=0.2, 'PSG8'=0, 'PTCHD4'=0, 'SH3TC2'=0.2)
#high.q <- c('PAPPA'=1, 'FLT1'=1, 'LEP'=1, 'PSG8'=1, 'PTCHD4'=1, 'SH3TC2'=1)

#low.q <- c('PAPPA'=0, 'FLT1'=0.2, 'LEP'=0.2, 'PSG8'=0, 'ERVFRD-1'=0, 'SH3TC2'=0.2)
#high.q <- c('PAPPA'=1, 'FLT1'=1, 'LEP'=1, 'PSG8'=1, 'ERVFRD-1'=1, 'SH3TC2'=1)

low.q <- c('DNMT1' = 0,'ERVFRD-1'= 0, 'PSG8'= 0, 'CGA'= 0, 'LEP'= 0.2,'CSHL1'= 0.3 , 'PAPPA'= 0,'FLT1'= 0.3)
high.q <- c('DNMT1' = 1,'ERVFRD-1'= 1, 'PSG8'= 1, 'CGA'= 1, 'LEP'= 1,'CSHL1'= 1 , 'PAPPA'= 1,'FLT1'= 1)




marker.gmat.df.cutoff <- marker.gmat.df[,1:3]
for(gene in marker.genes){
    #cat('gene ',gene)
    feature.value <- marker.gmat.df[,gene]
    #cutoff <- quantile(feature.value,probs = c(0,1) )
    cutoff <- quantile(feature.value,probs = c(low.q[gene],high.q[gene])) 
    cutoff.low <- cutoff[1]
    cutoff.high <- cutoff[2]
    feature.value[feature.value < cutoff.low] <- cutoff.low
    feature.value[feature.value > cutoff.high] <- cutoff.high
    marker.gmat.df.cutoff <- cbind(marker.gmat.df.cutoff,feature.value)
}
colnames(marker.gmat.df.cutoff)[4:ncol(marker.gmat.df.cutoff)] <- marker.genes


######
#options(repr.plot.height=5.8,repr.plot.width=5)
res.marker <- list()
for(gene in marker.genes){
    p <- ggplot(marker.gmat.df.cutoff,aes_string(x="UMAP_1",y="UMAP_2",col= paste0("`",gene,"`") )) +
    #p <- ggplot(marker.gmat.df,aes_string(x="UMAP_1",y="UMAP_2",col= paste0("`",gene,"`") )) +
      geom_point(size = .2,show.legend = TRUE,alpha= .5 ) +
      #scale_colour_manual(values = c('red','navy'))  +
      #cale_colour_gradientn(colors = viridis(6,option = 'C'))  +
      ##scale_colour_gradientn(colors = color_gradient_my )  +
      ##scale_colour_gradientn(colors = viridis(256, option = "D")) + #good
      scale_colour_gradientn(colors = color_tfdev )  +
      #theme_classic() +
      #theme(legend.position = 'top',axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14, face = "bold"),axis.title.y = element_text(size = 14, face = "bold")  ) +
      theme(
            #legend.position = 'none',
            legend.position = 'right',
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
      ggtitle(gene) +
      xlim(min(cluster.df.add$UMAP_1 )*1.5,1.5*max(cluster.df.add$UMAP_1 )) +
      #ggtitle(paste(sample, "nFeature_RNA",  sep=" ") ) +
      #guides(col = guide_legend(override.aes = list(size = 8))) +  ##no effect ??
      #ylim(0,1) + xlim(0,6) +  ##will overwrite scale_x/y_continuous
      #labs(x = "Fragments per barcode", y = "proportion of promoter reads")
      labs(x = "UMAP1", y = "UMAP2")
      #labs(x = "log10(pass_filters_reads)", y = "proportion of enhancer reads")
   # ggsave(filename = paste('pdfs/marker.gene_activity_score.plot.',gene,'.pdf',sep=''),height=5,width=5,useDingbats=FALSE)
    ggsave(filename = paste('pdfs/marker.gene_activity_score.plot.',gene,'.with_legend.pdf',sep=''),height=5,width=6,useDingbats=FALSE)
    res.marker[[gene]] <- p
    #print(p)
}

# ##arrange plot by patchwork
# options(repr.plot.height=15,repr.plot.width=15)
# res.marker[['DNMT1']] + res.marker[['CDH1']] + res.marker[['ERVFRD-1']] + res.marker[['FLT1']] +
# res.marker[['CSHL1']] + res.marker[['PSG8']] + res.marker[['VIM']] + res.marker[['PECAM1']] +
# res.marker[['CD14']] #+ res.marker[['MKI67']] + res.marker[['PLAC8']] + res.marker[['CD14']] +
# plot_layout(ncol=3,nrow=3)
# #print(res.marker[['FLT1']], vp=viewport(angle=-185))
# ##ggsave(filename = 'pdfs/marker.gene.common.umap.col.grad.my.pdf',height=13.5,width=18,useDingbats=FALSE) #can save grid multiple plots
# ggsave(filename = 'pdfs/marker.gene.common.umap.col.SnapATAC_like.pdf',height=13.5,width=18,useDingbats=FALSE)

# ##arrange plot by patchwork
# options(repr.plot.height=15,repr.plot.width=15)
# res.marker[['FLT1']] + res.marker[['LEP']] + res.marker[['INSIG2']] + res.marker[['CSHL1']] +
# res.marker[['CSH1']] + res.marker[['PAPPA']] + res.marker[['PSG1']] + res.marker[['CGA']] +
# res.marker[['GCM1']] 
# plot_layout(ncol=3,nrow=3)
# #print(res.marker[['FLT1']], vp=viewport(angle=-185))
# ggsave(filename = 'pdfs/marker.gene.hormone.bias.umap.pdf',height=13.5,width=16,useDingbats=FALSE) #can save grid multiple plots

##arrange plot by patchwork
options(repr.plot.height=15,repr.plot.width=15)
#     'PSG1','PSG8','CGA','CGB3',
#      'PAPPA','CSHL1', 'CSH1','CSH2',
#      'FLT1', 'ENG','LEP','INSIG2',
#     'ERVFRD-1','GCM1','DNMT1', 'CDH1',
#      'TEAD3','TEAD4','TEAD1','GATA3',
#     'CDKN1A', 'SPATA18', 'CROT', 'PTCHD4',
#     'DDX60','LIFR','MAPKAPK3','VAC14'

res.marker[['PSG1']] + res.marker[['PSG8']] + res.marker[['CGA']] + res.marker[['CGB3']] +
res.marker[['PAPPA']] + res.marker[['CSHL1']] + res.marker[['CSH1']] + res.marker[['CSH2']] +
res.marker[['FLT1']] + res.marker[['ENG']] + res.marker[['LEP']] + res.marker[['INSIG2']] +
res.marker[['ERVFRD-1']] + res.marker[['GCM1']] + res.marker[['DNMT1']] + res.marker[['CDH1']] 
plot_layout(ncol=4,nrow=4)
#print(res.marker[['FLT1']], vp=viewport(angle=-185))
##ggsave(filename = 'pdfs/marker.gene.hormone.bias.umap.pdf',height=15,width=15,useDingbats=FALSE) #can save grid multiple plots

options(repr.plot.height=12,repr.plot.width=15)
res.marker[['TEAD3']] + res.marker[['TEAD4']] + res.marker[['TEAD1']] + res.marker[['GATA3']] +
res.marker[['CDKN1A']] + res.marker[['SPATA18']] + res.marker[['CROT']] + res.marker[['PTCHD4']] +
res.marker[['DDX60']] + res.marker[['LIFR']] + res.marker[['MAPKAPK3']] + res.marker[['VAC14']] 
plot_layout(ncol=4,nrow=3)



##pick 4 x 3
options(repr.plot.height=12,repr.plot.width=15)
res.marker[['PSG8']] + res.marker[['CGA']] + res.marker[['PAPPA']] + res.marker[['CSHL1']] +
res.marker[['FLT1']] + res.marker[['LEP']] + res.marker[['ERVFRD-1']] + res.marker[['DNMT1']] +
res.marker[['CDKN1A']] + res.marker[['SPATA18']] + res.marker[['MAPKAPK3']] + res.marker[['LIFR']] 
plot_layout(ncol=4,nrow=3)

ggsave(filename = 'pdfs/marker.gene.picked.umap.viridisD256.pdf',height=12,width=15,useDingbats=FALSE)
#ggsave(filename = 'pdfs/marker.gene.picked.umap.my_grad_col.pdf',height=12,width=15,useDingbats=FALSE)



###do for probe gene
pdf(file = "pdfs/marker.gene_activity_score.combined.pdf",useDingbats = FALSE,height = 10, width = 15,compress = FALSE)
#options(repr.plot.height=10,repr.plot.width=15) #fix 5 x 12 for each row with 3 col
options(repr.plot.height=10,repr.plot.width=20)
patchwork::wrap_plots(res.marker, nrow = 2, ncol = 4)
dev.off()




###do gmat DAG with seurat doAllMarkers function (see gmat_DEG) ###



###do DAR identification to start peak centered downstream analysis##

##plot DAR on UMAP , Greenleaf's paper heatmap method?

idy.ls <- readRDS('DARs_exactTest/idy.ls.rds')
covs = Matrix::rowSums(x.after.sp@pmat);

saveRDS(covs,'DARs_exactTest/covs.te.rds')


pdf( "pdfs/PLA-8w-ATAC-DARs.visualize.pdf",height=15,width=10,useDingbats = FALSE)
options(repr.plot.width=10,repr.plot.height=15)
par(mfrow = c(3, 3));
for(cluster_i in levels(x.after.sp@cluster) ){
#for(cluster_i in levels(x.after.sp@cluster)[10:15] ){
	print(cluster_i)
	idy = idy.ls[[cluster_i]];
	vals = Matrix::rowSums(x.after.sp@pmat[,idy]) / covs;
	vals.zscore = (vals - mean(vals)) / sd(vals);
    peaknames = as.data.frame(x.after.sp@peak[idy])
#     write.table(peaknames,file = paste('DARs_exactTest/DARs.cluster',cluster_i,".exactTest.txt",sep=""),
#                sep='\t',quote=FALSE,row.names=FALSE)
    
	plotFeatureSingle(
		obj=x.after.sp,
		feature.value=vals.zscore,
		method="umap", 
		main=paste("Cluster ",cluster_i,sep = ""),
		point.size=0.1, 
		point.shape=19, 
		#down.sample=5000,
		quantiles=c(0.01, 0.99)
 		) 
#       scale_colour_gradientn(colors = color_tfdev   )  +
#       theme(legend.position = 'top',
#         axis.text=element_blank(), 
#         axis.title = element_text(size = 15, face = "bold"),
#         axis.ticks = element_blank(),
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.background = element_rect(color="black", fill = NA,size=1),
#         plot.title = element_text(size = 15, face = "bold"),
#         #complete = TRUE
#         plot.margin = unit(c(1,1,1,1), "lines")
#        )
  }
dev.off()




### use GREAT/ homer to annotate these DAR for (genes and function enrichment)

#output all DAR peak 

lapply(names(idy.ls), function(x){ length(idy.ls[[x]])  }  )
5285
4936
35011
1160
6207
39577
6153
2731
5251
53576
14217
19114
31985
23032
15772


# cluster <- x.sp@cluster #levels(x.sp@cluster)
# names(cluster) = x.sp@barcode
# cluster.mod <- factor(cluster,levels=c('9','6','2','11','7','10','3','4','8','5','1','12'))
# cluster.mod <- cluster.mod[!is.na(cluster.mod)] #6522
# cluster.mod = cluster.mod[order(cluster.mod)] #by cluster levels

order.cl <- c('6','3','9','5','7','1','2','4','8')





###################

