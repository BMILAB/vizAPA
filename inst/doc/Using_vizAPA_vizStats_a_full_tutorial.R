## ----readPACds, eval=TRUE, message=FALSE, warning=FALSE-----------------------
library(vizAPA)

data(scPACds, package='vizAPA')

# summary of the PACdataset
movAPA::summary(scPACds)

# cell meta data
head(scPACds@colData)

## ----odf4_gene, fig.show="hold", message=FALSE, warning=FALSE-----------------
gene=252868

## show the pAs in this gene
scPACds@anno[scPACds@anno$gene==gene, c(1:6, 10:12)]

## show expression level of each pAs
Matrix::rowSums(scPACds@counts[scPACds@anno$gene==gene, ])

## show the gene model and pAs in a track plot
library(TxDb.Mmusculus.UCSC.mm10.knownGene, quietly = TRUE)
txdb=TxDb.Mmusculus.UCSC.mm10.knownGene
annoSource=new("annoHub")
annoSource=addAnno(annoSource, txdb)

vizTracks(gene=gene, 
          PACds.list=list(pA=scPACds), PA.show=c("pos"),
          annoSource=annoSource,
          PA.columns="coord", PA.width=10,
          space5=1000, space3=1000)


## ----vizStats_violin1, fig.show="hold", message=FALSE, warning=FALSE----------
vizStats(scPACds, group='celltype', gene=gene, PAs=NULL, figType="violin")

## ----vizStats_violin2, fig.show="hold", message=FALSE, warning=FALSE----------
vizStats(scPACds, group='celltype', 
         gene=gene, PAs=NULL, 
         figType="violin",
         log=TRUE)

## ----vizStats_otherplots, fig.show='asis', message=FALSE, warning=FALSE-------

# boxplot
vizStats(scPACds, group='celltype', gene=gene, PAs=NULL, figType="box")

# violin plot with dots
vizStats(scPACds, group='celltype', gene=gene, PAs=NULL, figType="dot")

# bubble plot
vizStats(scPACds, group='celltype', gene=gene, PAs=NULL, 
         figType="bubble")

## ----vizStats_dot_PAs, fig.show="hold"----------------------------------------
# For example, here we show two PAs in another gene
PAids=c('PA6085', 'PA6083')
vizStats(scPACds, group='celltype', PAs=PAids, figType="dot")

## ----vizStats_violin_allPAs, fig.show="hold"----------------------------------
vizStats(scPACds, group='celltype', figType="violin")

## ----vizStats_statTheme, fig.show="hold"--------------------------------------
# display default statTHEME
# setStatTheme(NULL)
  
# To use another color palette 
vizStats(scPACds, group='celltype', gene=gene, figType="violin", 
         statTheme=list(group.cols=c(RColorBrewer::brewer.pal(8, "Set1"))))

## ----vizStats_selGroups, fig.show="hold"--------------------------------------
# change the order to RS>SC>ES
vizStats(scPACds, group='celltype', selGroups=c('RS','SC','ES')) 

## -----------------------------------------------------------------------------
colnames(scPACds@colData)

## ----vizUMAP1-----------------------------------------------------------------
vizUMAP(scPACds, group='celltype', xcol='UMAP_1', ycol='UMAP_2')

# Plot only the overlaying UMAP
vizUMAP(scPACds, group='celltype', annoUMAP=FALSE, xcol='UMAP_1', ycol='UMAP_2')

## ----vizUMAP_genes------------------------------------------------------------
vizUMAP(scPACds, group='celltype', xcol='UMAP_1', ycol='UMAP_2', genes=gene)

## ----vizUMAP_PAs--------------------------------------------------------------
# here we only show two cell types by specifying selGroup
vizUMAP(scPACds, group='celltype', xcol='UMAP_1', ycol='UMAP_2', 
        selGroups=c('SC','RS'), PAs=PAids)

## ----reduceDim, message=FALSE, warning=FALSE----------------------------------
## get embeddings using the pA count matrix with normalization
scPACds=reduceDim(scPACds, dims=1:10, dimLabel='umap_norm', norm=TRUE)

## without normalization
scPACds=reduceDim(scPACds, dims=1:10, dimLabel='umap_raw', norm=FALSE)


## ----reduceDim_res------------------------------------------------------------
## There are new columns adding to the colData slot.
head(scPACds@colData)

## ----vizUMAP_newdim, message=FALSE, warning=FALSE-----------------------------
vizUMAP(scPACds, group='celltype', xcol='umap_raw1', ycol='umap_raw2')

vizUMAP(scPACds, group='celltype', xcol='umap_norm1', ycol='umap_norm2')

## ----cal_RUD, warning=FALSE, message=FALSE, results='hide'--------------------
# First, calculate the RUD index for each gene. 
# Only genes with 3'UTR APA can be used for RUD calculation.
iPACds=getAPAindexPACds(scPACds, choose2PA='PD')

## ----getAPAmarkers, warning=FALSE, message=FALSE------------------------------
m=getAPAmarkers(iPACds,  group='celltype', everyPair = TRUE)
table(m$cluster1, m$cluster2)
head(m)

## ----vizAPAMarkers1, warning=FALSE, message=FALSE-----------------------------
# Visualize the top 6 APA markers, showing all the three cell types
vizAPAMarkers(iPACds, group='celltype', markers=m$rowid[1:6], figType = 'violin')

## ----vizAPAMarkers2, warning=FALSE, message=FALSE-----------------------------
# Actually the top 6 markers are between ES and SC, 
# so we can only show these two cell types
vizAPAMarkers(iPACds, group='celltype', markers=m$rowid[1:6], 
              selGroups=c('ES','SC'), figType = 'violin')

# To show genes as columns and cell types as rows
vizAPAMarkers(iPACds, group='celltype', markers=m$rowid[1:6], 
              figType = 'violin', statTheme=list(xgroup=FALSE))

# Plot a heatmap for APA markers
vizAPAMarkers(iPACds, group='celltype', markers=m$rowid[1:6], 
              figType = 'heatmap')

# Plot a bubble plot for markers
vizAPAMarkers(iPACds, group='celltype', markers=m$rowid[1:6], 
              figType = 'bubble')

# Switch the x/y axis of the bubble plot
vizAPAMarkers(iPACds, group='celltype', markers=m$rowid[1:6], 
              figType = 'bubble', statTheme=list(xgroup=FALSE))


## ----vizAPAMarkers_umap-------------------------------------------------------
# Plot the UMAP plot
vizAPAMarkers(iPACds, group='celltype', markers=m$rowid[1:6], 
              figType="umap", umap.x='UMAP_1', umap.y='UMAP_2')

# Do not show the UMAP embeddings 
vizAPAMarkers(iPACds, group='celltype', markers=m$rowid[1:6], 
              figType="umap", umap.x='UMAP_1', umap.y='UMAP_2', annoUMAP=F)

## ----vizAPAMarkers_sRUD-------------------------------------------------------
# first, get smartRUD APA index by movAPA
pd=movAPA::get3UTRAPApd(pacds=scPACds, 
                        minDist=50, maxDist=5000, 
                        minRatio=0.05, fixDistal=FALSE, 
                        addCols='pd') 

srud=movAPA::movAPAindex(pd, method="smartRUD", sRUD.oweight=FALSE)
head(srud[, 1:10])

# convert to PACds
iPACds=APAindex2PACds(srud, colData=pd@colData)

# get APA markers
m=getAPAmarkers(iPACds,  group='celltype', everyPair = TRUE)

# visualize top APA markers
vizAPAMarkers(iPACds, group='celltype', 
              markers=m$rowid[1:6], 
              figType = 'violin')

## ----getAPAmarkers_cnt, message=FALSE, warning=FALSE--------------------------
# The PACds is a PA-count dataset, 
# in which each row is the read count for each pA. 
# Here we detect differential expression for each pA as markers.
m=getAPAmarkers(scPACds,  group='celltype', everyPair = TRUE)
head(m)

## ----vizAPAMarkers_cnt--------------------------------------------------------
# Show the read count difference of the top 10 markers
vizAPAMarkers(scPACds, group='celltype', 
              markers=m$rowid[1:10], 
              figType="violin")

# Switch the x-y axis to show pA in columns and cell type in rows.
vizAPAMarkers(scPACds, group='celltype', 
              markers=m$rowid[1:10], 
              figType="violin", 
              statTheme=list(xgroup=FALSE))

## ----getAPAmarkers_each, message=FALSE, warning=FALSE-------------------------
# Detect markers between ES and all other cells.
m=getAPAmarkers(scPACds,  group='celltype', cluster1='ES')
table(m$cluster1, m$cluster2)

## ----message=FALSE, warning=FALSE---------------------------------------------
m=getAPAmarkers(scPACds,  group='celltype', cluster1='ES', cluster2='RS')
table(m$cluster1, m$cluster2)

## ----essage=FALSE,warning=FALSE-----------------------------------------------
m=getAPAmarkers(scPACds,  group='celltype', 
                cluster1='ES', cluster2=NULL, 
                everyPair = TRUE)
table(m$cluster1, m$cluster2)

## -----------------------------------------------------------------------------
sessionInfo()

