## ----demo_data, eval=TRUE, message=FALSE, warning=FALSE-----------------------
library(vizAPA)

data(scPACds, package='vizAPA')

# summary of the PACdataset
movAPA::summary(scPACds)

# cell meta data
head(scPACds@colData)


## ----coldata------------------------------------------------------------------
colnames(scPACds@colData)
head(scPACds@colData)

## ----vizUMAP0, eval=TRUE------------------------------------------------------
vizUMAP(scPACds, group='celltype', xcol='UMAP_1', ycol='UMAP_2')

## ----save_UMAP, eval=FALSE----------------------------------------------------
#  eoffice::topptx(filename = 'figures.pptx', title="vizUMAP_annotation_RUD",
#                  width = 8, height = 4, append=FALSE)

## ----vizUMAP_RUD, message=FALSE, warning=FALSE--------------------------------
# change the level order of the factor celltype in scPACds
# to make it consistent with sperm differentiation
scPACds@colData$celltype=factor(scPACds@colData$celltype, 
                                levels=c('SC','RS','ES'))

# get 3UTR PACds, only genes with 3'UTR APA can be used for RUD calculation.
# this scPACds has only 3UTR pAs, so this step will remove nothing
scPACds<-scPACds[scPACds@anno$ftr=="3UTR"]
# get RUD
RUD=getAPAindexPACds(scPACds, choose2PA = "PD")

## ----vizStats_dot, eval=TRUE, fig.show="hold"---------------------------------
vizStats(RUD, group='celltype', figType="dot")

## ----gene_convert, message=FALSE, warning=FALSE-------------------------------
head(RUD@anno)

library(Mus.musculus, quietly = TRUE)
orgdb=Mus.musculus
genes=getAnnoGenes(orgdb)
RUD@anno=merge(RUD@anno, genes, by.x='gene', by.y='gene_entrezid', all.x=TRUE)
# there are two entrez ids not in orgdb, so we use entrez id instead symbol
RUD@anno$gene_symbol[is.na(RUD@anno$gene_symbol)]=
  RUD@anno$gene[is.na(RUD@anno$gene_symbol)]
RUD@anno$entrezid=RUD@anno$gene
# set the gene column as gene_symbol
RUD@anno$gene=RUD@anno$gene_symbol
# set rownames of counts slot as gene_symbol
rownames(RUD@counts)=RUD@anno$gene_symbol
# after conversion
head(RUD@anno)

## ----ASRGL1-------------------------------------------------------------------
gene='Asrgl1'
RUD@anno[RUD@anno$gene==gene, ]

## ----vizStats_heatmap, eval=TRUE, fig.show="hold", message=FALSE, warning=FALSE----
## plot a heatmap for summarizing average expression of pAs in this gene
vizStats(RUD, group='celltype', gene=gene, figType="heatmap")

## we can change the text size or remove text in the heatmap (heat.text.size=0)
## or change the number of digits for the averaged value
vizStats(RUD, group='celltype', gene=gene, figType="heatmap", 
         statTheme = list(heat.text.size=6, heat.text.digits=4))


## ----vizStats_heatmap_ggplot2, eval=TRUE, fig.show="hold", message=FALSE, warning=FALSE----
p=vizStats(RUD, group='celltype', gene=gene, figType="heatmap")
p+ggplot2::scale_fill_gradientn(colours = grDevices::heat.colors(50)) + 
  ggplot2::theme_minimal()

## ----vizStats_boxplot, eval=TRUE, fig.show="hold", message=FALSE, warning=FALSE----
vizStats(RUD, group='celltype', gene=gene, figType="box")

## ----vizStats_other, eval=TRUE, fig.show='asis', message=FALSE, warning=FALSE----
# violin plot
vizStats(RUD, group='celltype', gene=gene, figType="violin")

# violin plot with dots, and each dot is one cell
vizStats(RUD, group='celltype', gene=gene, figType="dot")

#eoffice::topptx(filename = 'figures.pptx', title="gene_dot_plot", 
#                width = 4, height = 4, append=TRUE)

## ----vizStats_bubble, eval=TRUE, fig.show='asis', message=FALSE, warning=FALSE----
vizStats(RUD, group='celltype', gene=gene, figType="bubble")

## ----vizStats_heatmap_count, eval=TRUE, fig.show="hold", message=FALSE, warning=FALSE----
# for scPACds, the id type is entrezid
geneid=66514
## plot a heatmap for summarizing average expression of pAs in this gene
vizStats(scPACds, group='celltype', gene=geneid, figType="heatmap")

## ----vizStats_box_count, eval=TRUE, fig.show="hold", message=FALSE, warning=FALSE----
# plot a boxplot to compare the pA usage of this gene in different cell types
vizStats(scPACds, group='celltype', gene=geneid, figType="box")

#eoffice::topptx(filename = 'figures.pptx', title="gene_box_plot", 
#                width = 4, height = 4, append=TRUE)

## ----vizStats_bubble_count, eval=TRUE, fig.show='asis', message=FALSE, warning=FALSE----
vizStats(scPACds, group='celltype', gene=geneid,  figType="bubble")

## ----vizStats_dot_count, eval=TRUE, fig.show="hold"---------------------------
vizStats(scPACds, group='celltype', figType="dot")

## ----vizUMAP_RUD_all----------------------------------------------------------
# show UMAP using RUD scores
vizUMAP(RUD, group='celltype', xcol='UMAP_1', ycol='UMAP_2')

#eoffice::topptx(filename = 'figures.pptx', title="vizUMAP_annotation_RUD", 
#                width = 8, height = 4, append=TRUE)

## ----vizUMAP_agene------------------------------------------------------------
vizUMAP(RUD, 
        group='celltype', xcol='UMAP_1', 
        ycol='UMAP_2', genes=gene)

## ----vizUMAP_agene2-----------------------------------------------------------
vizUMAP(RUD, 
        group='celltype', xcol='UMAP_1', ycol='UMAP_2', 
        annoUMAP=TRUE,
        genes=gene,
        statTheme = list(scale.low.col='green', scale.high.col='red'))


## ----vizStats_agene-----------------------------------------------------------
vizStats(RUD, group='celltype', figType="dot", gene=gene)

## ----vizUMAP_aPA--------------------------------------------------------------
PAids='PA16510'
vizUMAP(scPACds, group='celltype', xcol='UMAP_1', ycol='UMAP_2', 
        PAs=PAids, annoUMAP = FALSE)

## ----getAPAmarkers_count, eval=FALSE------------------------------------------
#  m=getAPAmarkers(scPACds,  group='celltype', everyPair = TRUE)

## ----getAPAmarkers, warning=FALSE, message=FALSE------------------------------
## obtain APA markers by wilcox.test for each pair of cell types
m=getAPAmarkers(RUD,  group='celltype', everyPair = TRUE,  
                min.pct = 0.1, logFC=0.12)

## show marker details
head(m)

## ----bar_markerNum------------------------------------------------------------
countAPAmarkers(m)

#eoffice::topptx(filename = 'figures.pptx', title="marker_num", 
#                width = 8, height = 4, append=TRUE)

## ----markers_list, warning=FALSE, message=FALSE-------------------------------
markers=m$rowid[1:6]

## ----vizAPAmarkers, warning=FALSE, message=FALSE------------------------------
vizAPAMarkers(RUD, group='celltype', 
              markers=markers, 
              figType = 'violin')

#eoffice::topptx(filename = 'figures.pptx', title="markers_violin", 
#                width = 8, height = 6, append=TRUE)

## ----vizAPAmarkers_others, warning=FALSE, message=FALSE-----------------------

# we can plot only two cell types, here SC and RS
vizAPAMarkers(RUD, group='celltype', 
              markers=markers, 
              selGroups = c('SC','RS'),
              figType = 'violin')

# Plot a heatmap for APA markers
vizAPAMarkers(RUD, group='celltype', 
              markers=markers, 
              figType = 'heatmap')

# Plot a bubble plot for markers
vizAPAMarkers(RUD, group='celltype', 
              markers=markers, 
              figType = 'bubble')

# Switch the x/y axis of the bubble plot
vizAPAMarkers(RUD, group='celltype', 
              markers=markers, 
              figType = 'bubble', 
              statTheme=list(xgroup=FALSE))


## ----vizAPAMarkers_umap2, out.height="90%",out.width="90%"--------------------
# Plot the UMAP plot
vizAPAMarkers(RUD, 
              group='celltype', 
              markers=markers,
              figType="umap", 
              umap.x='UMAP_1', umap.y='UMAP_2')

## ----viz_single_marker--------------------------------------------------------
vizStats(RUD, group='celltype', figType="dot", gene=markers[1])

## To plot only the two cell types involved in the marker
vizStats(RUD, group='celltype', figType="dot", gene=markers[1], 
         selGroups=c('SC','RS'))

## plot a heatmap to average RUD values across cell types
vizStats(RUD, group='celltype', figType="heatmap", gene=markers[1])



## ----getAPAmarkers_SC---------------------------------------------------------
m=getAPAmarkers(RUD,  group='celltype', cluster1='SC', only.pos = TRUE)

## -----------------------------------------------------------------------------
head(m)
countAPAmarkers(m, plot=F)

## ----vizAPAmarkers2, warning=FALSE, message=FALSE-----------------------------
vizAPAMarkers(RUD, group='celltype', 
              markers=m$rowid[1:10], 
              figType = 'violin')

## ----readBam------------------------------------------------------------------
#Create the list of BAM files
bam.files=c("dedup_GSM2803334.ES.mini.sorted.bam",
            "dedup_GSM2803334.RS.mini.sorted.bam",
            "dedup_GSM2803334.SC.mini.sorted.bam")
bam.groups=c("ES","RS","SC")
bam.labels=c("ES" ,"RS","SC")
bam.path='./'

# here we set bam.order to make BAM colors consistent with PACdataset
bams<-readBAMFileNames(bam.files=bam.files, 
                       bam.path=bam.path, 
                       bam.labels = bam.labels, 
                       bam.groups = bam.groups,
                       bam.order=c('SC','RS','ES'))
bams

## ----annoHub, warning=FALSE, message=FALSE------------------------------------
annoSource=new("annoHub")
library(TxDb.Mmusculus.UCSC.mm10.knownGene, quietly = TRUE)
txdb=TxDb.Mmusculus.UCSC.mm10.knownGene
annoSource=addAnno(annoSource, txdb)
annoSource

## ----vizTracks, fig.show="hold", message=FALSE, warning=FALSE-----------------
vizTracks(gene=geneid, 
          bams=bams, 
          PACds.list=list(PA=scPACds), 
          PA.show=c("pos"),
          cells=TRUE, cells.group='celltype',
          cells.width=200, cells.annoCoord=NULL,
          cells.method=c('sum'), cells.sort=c('group'),
          annoSource=annoSource,
          PA.columns="coord", PA.width=10,
          space5=0, space3=100,
          )


## ----tracks, warning=FALSE, message=FALSE-------------------------------------
# gene model track
gmtrack=getTrackGeneModel(gene=geneid, annoSource=annoSource, title=gene)


# PA track
pactrack=getTrackPACds(scPACds, gene=geneid, 
                      PA.show=c("pos"), 
                      PA.columns="coord", PA.width=10,
                      title='PA') 
# we can customize the track using ggplot2's grammer
# pactrack$PA=pactrack$PA+ggplot2::theme_minimal()

# cell track
celltrack<-getTrackCells(scPACds, group='celltype', gene=geneid,
                         sortMethod='sum', sortCells='group', log=FALSE,
                         PA.columns="coord", PA.width=50, title='sum.cells')

# bam track
genomicRegion=getGenesRange(gene=geneid, annoSource, rt='list')
# add 1000 bp to both ends of the region
genomicRegion$start=genomicRegion$start-1000
genomicRegion$end=genomicRegion$end+1000
genomicRegion=paste0(unlist(genomicRegion), collapse=':')
bamtrack<-getTrackBams(bams, genomicRegion=genomicRegion)

# combine tracks
tks<-c(gmtrack, pactrack, celltrack, bamtrack)
names(tks)[1:3]=c(gene, 'pA', 'cells')
ggbio::tracks(tks)

## ----tracks_zoom, warning=FALSE, message=FALSE--------------------------------
# plot legend ‘horizontal’ or ‘vertical’
ggbio::tracks(tks, heights=c(0.7, 0.7, 1.2, 1.2, 1.2, 1.2),
              xlim=c(9111000, 9113500))+ 
  ggplot2::theme(legend.direction = 'vertical')

#eoffice::topptx(filename = 'figures.pptx', title="tracks", 
#                width = 8, height = 8, append=TRUE)

## -----------------------------------------------------------------------------
sessionInfo()

