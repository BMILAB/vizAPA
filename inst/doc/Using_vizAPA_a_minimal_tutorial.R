## ----demo_data, eval=TRUE, message=FALSE, warning=FALSE-----------------------
library(vizAPA)

data(scPACds, package='vizAPA')

# summary of the PACdataset
movAPA::summary(scPACds)

# cell meta data
head(scPACds@colData)

## ----vizUMAP0, eval=TRUE------------------------------------------------------
vizUMAP(scPACds, group='celltype', xcol='UMAP_1', ycol='UMAP_2')

## ----vizStats_boxplot, eval=TRUE, fig.show="hold", message=FALSE, warning=FALSE----
## first, we check the gene column in the anno slot
head(scPACds@anno[, c('chr','strand','coord','gene','ftr')])

## the gene ID is entrez id, so we use 252868 instead of Odf4
gene='252868'

# plot a boxplot to compare the pA usage of this gene in different cell types
vizStats(scPACds, group='celltype', gene=gene, PAs=NULL, figType="box")

## ----vizStats_other, eval=TRUE, fig.show='asis', message=FALSE, warning=FALSE----
# violin plot
vizStats(scPACds, group='celltype', gene=gene, PAs=NULL, figType="violin")

# violin plot with dots
vizStats(scPACds, group='celltype', gene=gene, PAs=NULL, figType="dot")

# bubble plot
vizStats(scPACds, group='celltype', gene=gene, PAs=NULL, figType="bubble")

## ----vizStats_violin, eval=TRUE, fig.show="hold"------------------------------
vizStats(scPACds, group='celltype', figType="violin")

## -----------------------------------------------------------------------------
gene='252868'
# First, we check the coordinate labels of the 2D-embedding. 
# For this data, the labels are UMAP_1 and UMAP_2.
colnames(scPACds@colData)
head(scPACds@colData)

## ----vizUMAP_all--------------------------------------------------------------
# Plot the UMAP plot showing cell clusters and another UMAP plot overlaying 
# with the mean expression value of pAs in each cell.
vizUMAP(scPACds, group='celltype', xcol='UMAP_1', ycol='UMAP_2')

## ----vizUMAP_agene------------------------------------------------------------
## Here we plot the Odf4 gene, 
## overlaying the mean value of all pAs in this gene
vizUMAP(scPACds, 
        group='celltype', xcol='UMAP_1', 
        ycol='UMAP_2', genes=gene)

## ----vizUMAP_aPA--------------------------------------------------------------
# here we only overlay one pA in Odf4 gene
# here only show two cell types by specifying selGroup
PAids='PA13514'
vizUMAP(scPACds, group='celltype', xcol='UMAP_1', ycol='UMAP_2', selGroups=c('SC','RS'), PAs=PAids)

## ----getAPAindex, message=FALSE-----------------------------------------------
# First, calculate the RUD index for each gene. 
# Only genes with 3'UTR APA can be used for RUD calculation.
iPACds=getAPAindexPACds(scPACds, choose2PA='PD')
head(iPACds[, 1:10])

## ----vizAPAmarkers, warning=FALSE, message=FALSE------------------------------
## obtain APA markers by wilcox.test for each pair of cell types
m=getAPAmarkers(iPACds,  group='celltype', everyPair = TRUE)

## show marker numbers
table(m$cluster1, m$cluster2)

## show marker details
head(m)

# Visualize the top 6 APA markers, showing all the three cell types
vizAPAMarkers(iPACds, group='celltype', 
              markers=m$rowid[1:6], 
              figType = 'violin')

# Plot a heatmap for APA markers
vizAPAMarkers(iPACds, group='celltype', 
              markers=m$rowid[1:6], 
              figType = 'heatmap')

# Plot a bubble plot for markers
vizAPAMarkers(iPACds, group='celltype', 
              markers=m$rowid[1:6], 
              figType = 'bubble')

# Switch the x/y axis of the bubble plot
vizAPAMarkers(iPACds, group='celltype', 
              markers=m$rowid[1:6], 
              figType = 'bubble', 
              statTheme=list(xgroup=FALSE))


## ----vizAPAMarkers_umap2, out.height="90%",out.width="90%"--------------------
# Plot the UMAP plot
vizAPAMarkers(iPACds, 
              group='celltype', 
              markers=m$rowid[1:6], 
              figType="umap", 
              umap.x='UMAP_1', umap.y='UMAP_2')

## ----getAPAmarkers_ES---------------------------------------------------------
# Detect markers between ES and all other cells.
m=getAPAmarkers(scPACds,  group='celltype', cluster1='ES')

## -----------------------------------------------------------------------------
head(m)
table(m$cluster1, m$cluster2)

## ----readBam------------------------------------------------------------------
#Create the list of BAM files
bam.files=c("dedup_GSM2803334.ES.mini.sorted.bam",
            "dedup_GSM2803334.RS.mini.sorted.bam",
            "dedup_GSM2803334.SC.mini.sorted.bam")
bam.groups=c("ES","RS","SC")
bam.labels=c("ES" ,"RS","SC")
bam.path='./'

bams<-readBAMFileNames(bam.files=bam.files, 
                       bam.path=bam.path, 
                       bam.labels = bam.labels, 
                       bam.groups = bam.groups)

bams

## ----annoHub, warning=FALSE, message=FALSE------------------------------------
annoSource=new("annoHub")
library(TxDb.Mmusculus.UCSC.mm10.knownGene, quietly = TRUE)
txdb=TxDb.Mmusculus.UCSC.mm10.knownGene
annoSource=addAnno(annoSource, txdb)
annoSource

## ----vizTracks, fig.show="hold", message=FALSE, warning=FALSE-----------------
gene=27058
vizTracks(gene=gene, 
          bams=bams, 
          PACds.list=list(PA=scPACds), 
          PA.show=c("pos"),
          annoSource=annoSource,
          PA.columns="coord", PA.width=10,
          space5=1000, space3=1000)


## -----------------------------------------------------------------------------
sessionInfo()

