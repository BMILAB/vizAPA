## ----readPACds, eval=TRUE, message=FALSE, warning=FALSE-----------------------
library(vizAPA)

data(scPACds, package='vizAPA')

# summary of the PACdataset
movAPA::summary(scPACds)

# cell meta data
head(scPACds@colData)

## ----notrun, eval=FALSE-------------------------------------------------------
#  # First, use readPACds or createPACdataset to create the PACdataset object
#  # from a pA list, pA-cell counts table, cell meta data
#  PACds=movAPA::readPACds(pacFile,
#                  colDataFile)
#  
#  movAPA::summary(PACds)
#  
#  # If there are internal priming artifacts, should be removed first
#  library(BSgenome.Mmusculus.UCSC.mm10)
#  bsgenome=BSgenome.Mmusculus.UCSC.mm10
#  PACds=movAPA::removePACdsIP(PACds, bsgenome)
#  PACds=PACds$real
#  
#  # Then annotate the pA list with genome annotation from TxDb or a GFF/GTF file
#  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
#  txdb=TxDb.Mmusculus.UCSC.mm10.knownGene
#  
#  PACds=movAPA::annotatePAC(PACds, aGFF=txdb)
#  
#  # Then extend 3'UTR by some length, e.g., 1000 bp
#  # which can recuite potential pAs in downstream intergenic regions of 3'UTR
#  PACds=movAPA::ext3UTRPACds(PACds, 1000)
#  
#  # For single cell data, we'd better only use 3'UTR pAs for analysis
#  # since pAs in other genomic regions may be artifacts
#  PACds=movAPA::get3UTRAPAds(PACds)
#  
#  # Finally, a full annotated 3'UTR PACdataset can be obtained
#  # which could be used in vizAPA for visualization
#  movAPA::summary(PACds)

## ----readBAM------------------------------------------------------------------
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

## ----annoHub------------------------------------------------------------------
annoSource=new("annoHub")

## ----addAnno_gff_notrun, eval=FALSE-------------------------------------------
#  gff <- useGff(gff ="gencode.vM32.annotation.gff3")
#  annoSource=addAnno(annoSource, gff)

## ----addAnno_orgdb, results='hide', warning=FALSE, message=FALSE--------------
library(Mus.musculus, quietly = TRUE)
orgdb=Mus.musculus
annoSource=addAnno(annoSource, orgdb)
annoSource=setDefaultAnno(annoSource, 'orgdb')

## ----addAnno_txdb, warning=FALSE, message=FALSE-------------------------------
library(TxDb.Mmusculus.UCSC.mm10.knownGene, quietly = TRUE)
txdb=TxDb.Mmusculus.UCSC.mm10.knownGene
annoSource=addAnno(annoSource, txdb)

## ----addAnno_ensdb, warning=FALSE, message=FALSE------------------------------
library(EnsDb.Mmusculus.v79, quietly = TRUE)
ensdb=EnsDb.Mmusculus.v79
annoSource=addAnno(annoSource, ensdb)

## ----addAnno_biomart, warning=FALSE, message=FALSE----------------------------
library(biomaRt)
bm = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
annoSource=addAnno(annoSource, bm)

## ----addAnno_genes------------------------------------------------------------
genes=getAnnoGenes(orgdb)
annoSource=addAnno(annoSource, genes)

## ----annoHub_order------------------------------------------------------------
## re-order annoSource to set the searching priority
names(annoSource@annos)
annoSource@annos=annoSource@annos[c('orgdb','txdb','ensdb','genes','biomart')]
annoSource

## ----check_chrs, warning=FALSE, message=FALSE---------------------------------
## The chr names are different: orgdb is chr1, ensdb is 1.
isChrConsistent(annoSource['orgdb'], annoSource['ensdb'], exact=TRUE) 

isChrConsistent(annoSource['orgdb'], annoSource['ensdb'], exact=FALSE) 

## The chr names are total the same between txdb and orgdb, both are like chr1
isChrConsistent(annoSource['txdb'], annoSource['orgdb'], exact=TRUE) 

## ----getChrs, warning=FALSE, message=FALSE------------------------------------
# chr names are the same in PACds and in BAM file
isChrConsistent(scPACds, bams, exact=TRUE)

# Here for annoSource, the default anno was used. 
# It seems that some chromosomes do not have exacytly the same name
isChrConsistent(scPACds, bams, annoSource, exact=TRUE) 
# However, the main chr names are the same
isChrConsistent(scPACds, bams, annoSource, exact=FALSE) 

# We then get full chromosome names of each object, 
# and it would be fine that main chromosomes have the same name.
getChrs(scPACds)
getChrs(bams)

## ----getChrs_annos, eval=FALSE------------------------------------------------
#  getChrs(annoSource, which='orgdb') # chr1
#  getChrs(annoSource, which='ensdb') # 1
#  getChrs(annoSource, which='biomart') # 1
#  getChrs(annoSource, which='genes') # chr1
#  getChrs(annoSource, which='txdb') # chr1

## ----chr_mapping--------------------------------------------------------------
chrMappings=data.frame(cn1=c(1:19, 'X','Y'))
chrMappings$cn2=paste0('chr', chrMappings$cn1)
annoSource@chrMappings=chrMappings
head(chrMappings)

## ----select_APA---------------------------------------------------------------
gene=66514

## show this gene in PACds
scPACds@anno[scPACds@anno$gene==gene, c(1:6, 10:12)]

## the gene is represented as Entrez ID in PACds, we can show its symbol name 
## using 'genes' table in the annos slot
annoSource@annos$genes[annoSource@annos$genes$gene_entrezid==gene, ]

## ----vizTracks1, fig.show="hold", message=FALSE, warning=FALSE----------------
vizTracks(gene=gene, 
          bams=bams, PACds.list=list(pA=scPACds), PA.show=c("pos"),
          annoSource=annoSource,
          PA.columns="coord", PA.width=10,
          space5=1000, space3=1000)


## ----vizTracks1_cells, message=FALSE, warning=FALSE---------------------------
vizTracks(gene=gene, 
          bams=bams, PACds.list=list(pA=scPACds), PA.show=c("pos"),
          cells=TRUE, cells.group='celltype',
          cells.method=c('sum'), cells.sort=c('group'),
          cells.width=100,
          annoSource=annoSource,
          PA.columns="coord", PA.width=10, logPA=TRUE,
          space5=1000, space3=1000)


## ----vizTheme, eval=FALSE-----------------------------------------------------
#  # Show default vizTheme parameters
#  setVizTheme(NULL)

## ----vizTracks1_colors, message=FALSE, warning=FALSE--------------------------
# first show the order of the bams
bams

# set colors of the cells track as vizTheme's bams.col
cells.group.cols=c(RColorBrewer::brewer.pal(3, "Set2"))
names(cells.group.cols)=c('ES','RS','SC')
vizTheme=list(cells.group.cols=cells.group.cols)

vizTracks(gene=gene, 
          bams=bams, PACds.list=list(pA=scPACds), PA.show=c("pos"),
          cells=TRUE, cells.group='celltype',
          cells.method=c('sum'), cells.sort=c('group'),
          cells.width=100,
          annoSource=annoSource,
          PA.columns="coord", PA.width=10, logPA=TRUE,
          space5=1000, space3=1000,
          vizTheme=vizTheme)

## ----get_gene_gr, message=FALSE, warning=FALSE--------------------------------
# here for demonstration, we get the genomic region of the gene
genesGR=getGenesRange(gene, annoSource, rt='str')
genesGR

## ----vizTracks2, warning=FALSE, message=FALSE---------------------------------
vizTracks(genomicRegion=genesGR, 
          bams=bams, PACds.list=list(pA=scPACds), PA.show=c("pos"),
          cells=TRUE, cells.group='celltype',
          cells.method=c('sum'), cells.sort=c('group'),
          cells.width=100,
          annoSource=annoSource,
          PA.columns="coord", PA.width=10, logPA=TRUE,
          space5=1000, space3=1000,
          vizTheme=vizTheme)

## ----vizTracks_PA_region, fig.show="hold", message=FALSE, warning=FALSE-------
# show the region of the pA
scPACds@anno[1:5, 1:5]

# plot pA regions by setting PA.columns
vizTracks(genomicRegion=genesGR, 
          bams=bams, PACds.list=list(pA=scPACds), PA.show=c("pos"),
          cells=FALSE,
          annoSource=annoSource,
          PA.columns="start:end", PA.width=10, logPA=TRUE,
          space5=1000, space3=1000,
          vizTheme=vizTheme)

## ----vizTracks_PA_value, fig.show="hold", message=FALSE, warning=FALSE--------
## Here we add a tot_tagnum to the PACds@anno, 
## and show the total expression of pAs
scPACds@anno$tot_tagnum=Matrix::rowSums(scPACds@counts)

## We can do log-transformation to show the tot_tagnum clearly.
scPACds@anno$tot_tagnum=log2(scPACds@anno$tot_tagnum+1)

vizTracks(genomicRegion=genesGR, 
          bams=bams, PACds.list=list(pA=scPACds), 
          PA.show=c("pos", "tot_tagnum"),
          annoSource=annoSource,
          PA.columns="coord", PA.width=10, logPA=TRUE,
          space5=1000, space3=1000,
          vizTheme=vizTheme)

## ----vizTracks_multi_PA, fig.show='asis', message=FALSE, warning=FALSE--------

# here we just use a replicated PACds for demonstration
# we can save the tracks and then plot
tks=vizTracks(genomicRegion=genesGR, 
          bams=bams, PACds.list=list(pA1=scPACds, PA2=scPACds), 
          PA.show=c("pos", "tot_tagnum"),
          annoSource=annoSource,
          PA.columns="coord", PA.width=10, logPA=TRUE,
          space5=1000, space3=1000,
          vizTheme=vizTheme, res='list')

# plot tracks
ggbio::tracks(tks)

# we can change the title of the tracks
names(tks)[2:5]=c('pos1', 'tag1', 'pos2', 'tag2')
ggbio::tracks(tks)

## ----vizTracks_PA_colors, fig.show="hold", message=FALSE, warning=FALSE-------
vizTheme=list(PA.col='green', PA.fill='green')
vizTracks(gene=gene,
          bams=bams, PACds.list=list(pA1=scPACds), PA.show=c("pos"),
          annoSource=annoSource,
          PA.columns="coord", PA.width=10,
          space5=1000, space3=1000,
          vizTheme=vizTheme)

## ----vizTracks_cells_1, fig.show='asis', message=FALSE, warning=FALSE---------
# use diffusion map to sort cells within each group
vizTracks(genomicRegion=genesGR, 
          PACds.list=scPACds, PA.show=c("pos"),
          annoSource=annoSource,
          cells=TRUE, cells.group='celltype', 
          cells.method='diff', cells.sort='group', 
          cells.width=100,
          PA.columns="coord", PA.width=10)

## ----vizTracks_cells_2, fig.show='asis', message=FALSE, warning=FALSE---------
vizTracks(genomicRegion=genesGR, 
          PACds.list=scPACds, PA.show=c("pos"),
          annoSource=annoSource,
          cells=TRUE, cells.group='celltype', 
          cells.method='sum', cells.sort='group', 
          cells.width=100,
          PA.columns="coord", PA.width=10, logPA=TRUE)

## ----vizTracks_cells_3, fig.show='asis', message=FALSE, warning=FALSE---------
vizTheme=list(cells.annoBar.pos='right', 
              cells.scale.low='white', cells.scale.high='purple', cells.scale.mid = "grey",
              cells.group.cols=c(RColorBrewer::brewer.pal(8, "Set1")))

vizTracks(genomicRegion=genesGR, 
          PACds.list=scPACds, PA.show=c("pos"),
          annoSource=annoSource,
          cells=TRUE, cells.group='celltype', 
          cells.method='sum', cells.sort='group', 
          cells.width=100,
          PA.columns="coord", PA.width=10, logPA=TRUE, vizTheme=vizTheme)

## ----vizTracks_cells_4, fig.show='asis', message=FALSE, warning=FALSE---------
vizTracks(genomicRegion=genesGR, 
          PACds.list=scPACds, PA.show=c("pos"),
          annoSource=annoSource,
          cells=TRUE, cells.group='celltype', 
          cells.method='sum', cells.sort='group', 
          cells.width=100, cells.annoCoord=9110000,
          PA.columns="coord", PA.width=10, logPA=TRUE)

## ----vizTracks_cells_5, fig.show='asis', message=FALSE, warning=FALSE---------
vizTracks(genomicRegion=genesGR, 
          PACds.list=scPACds, PA.show=c("pos"),
          annoSource=annoSource,
          cells=TRUE, cells.group='celltype', 
          cells.method='sum', cells.sort='none', 
          cells.width=100, 
          PA.columns="coord", PA.width=10, logPA=TRUE)

## ----vizTracks_cells_6, fig.show='asis', message=FALSE, warning=FALSE---------
vizTracks(genomicRegion=genesGR, 
          PACds.list=scPACds, PA.show=c("pos"),
          annoSource=annoSource,
          cells=TRUE, cells.group='celltype', 
          cells.method='sum', cells.sort='all', 
          cells.width=100, 
          PA.columns="coord", PA.width=10, logPA=TRUE)


## ----vizTracks_BAM1, fig.show="hold", message=FALSE, warning=FALSE------------
# set how to show BAM coverage
vizTheme=list(bam.covMerge='avg', 
              bam.showMerge=TRUE, 
              bam.showSingle=TRUE, 
              bam.covMerge.col='yellow', 
              bam.covMerge.fill='yellow')

# we can also set the colors for individual bam coverages
vizTheme$bams.col=RColorBrewer::brewer.pal(8, "Set1")
vizTheme$bams.fill=RColorBrewer::brewer.pal(8, "Set1")

vizTracks(gene=gene,
          bams=bams, 
          PACds.list=list(pA1=scPACds), PA.show=c("pos"),
          annoSource=annoSource,
          PA.columns="coord", PA.width=10,
          space5=1000, space3=1000,
          vizTheme=vizTheme)

## ----vizTracks_BAM2, fig.show="hold", message=FALSE, warning=FALSE------------
vizTheme$bam.collapse=TRUE

vizTracks(gene=gene,
          bams=bams, 
          PACds.list=list(pA1=scPACds), PA.show=c("pos"),
          annoSource=annoSource,
          PA.columns="coord", PA.width=10,
          space5=1000, space3=1000,
          vizTheme=vizTheme)

## ----check_chrs2, message=FALSE, warning=FALSE--------------------------------
## yes, they are consistent
isChrConsistent(annoSource['txdb'], bams, scPACds, exact=FALSE) 
head(getChrs(annoSource['txdb']))

## set the default anno as txdb
annoSource=setDefaultAnno(annoSource, 'txdb')

## ----vizTracks_gene1, message=FALSE, warning=FALSE----------------------------
vizTracks(gene=gene, 
          bams=bams, 
          PACds.list=list(pA=scPACds), PA.show=c("pos"),
          annoSource=annoSource,
          PA.columns="coord", PA.width=10,
          space5=1000, space3=1000)

## ----vizTracks_gene2, message=FALSE, warning=FALSE----------------------------
vizTracks(gene=gene, 
          bams=bams, 
          PACds.list=list(pA=scPACds), PA.show=c("pos"),
          annoSource=new('annoHub'),
          PA.columns="coord", PA.width=10,
          space5=1000, space3=1000)


## ----get_gene_grs, eval=FALSE-------------------------------------------------
#  # first check whether the gene is annotated in different annos.
#  getGenesRange(genes=gene, annoSource['txdb'], rt='gr')
#  getGenesRange(genes=gene, annoSource['ensdb'], rt='gr')
#  getGenesRange(genes=gene, annoSource['orgdb'], rt='gr')
#  getGenesRange(genes=gene, annoSource['biomart'], rt='gr')
#  getGenesRange(genes=gene, annoSource['genes'], rt='gr')

## ----get_indiv_tracks_gm, message=FALSE, warning=FALSE------------------------
# use the orgdb anno
annoSource=setDefaultAnno(annoSource, 'orgdb')
gm.tk1=getTrackGeneModel(gene=gene, annoSource=annoSource, title='orgdb')

annoSource=setDefaultAnno(annoSource, 'ensdb') 
gm.tk2=getTrackGeneModel(gene=gene, annoSource=annoSource, title='ensdb')

annoSource=setDefaultAnno(annoSource, 'txdb') 
gm.tk3=getTrackGeneModel(gene=gene, annoSource=annoSource, title='txdb') 

annoSource=setDefaultAnno(annoSource, 'genes') 
gm.tk4=getTrackGeneModel(gene=gene, annoSource=annoSource, title='genes') 

annoSource=setDefaultAnno(annoSource, 'biomart') 
gm.tk5=getTrackGeneModel(gene=gene, annoSource=annoSource, title='biomart') 


## ----plot_indiv_tracks_gm, message=FALSE, warning=FALSE-----------------------
# show these different annotations for the same gene
tks=c(gm.tk1, gm.tk2, gm.tk3, gm.tk4, gm.tk5)
ggbio::tracks(tks)

## ----get_PA_tracks, message=FALSE, warning=FALSE------------------------------
pac.tk1=getTrackPACds(scPACds, gene=gene, 
                      PA.show=c("pos"), 
                      PA.columns="coord", PA.width=10,
                      title='PApos', vizTheme = list(PA.shape=17)) 

pac.tk2=getTrackPACds(scPACds, gene=gene,
                      PA.show=c("tot_tagnum"), log=TRUE,
                      PA.columns="start:end", title='totTag')

sc.tk1=getTrackCells(scPACds,gene=gene,
                     group='celltype', 
                     sortMethod='sum', sortCells='group',log=TRUE,
                     PA.columns="coord", PA.width=100, title='diffs')

genomicRegion=getGenesRange(genes=gene, annoObj=annoSource['orgdb'], rt='str')
bam.tk1=getTrackBams(bams, genomicRegion=genomicRegion) 

## show gene model, pAs, cells, and bams
tks=c(gm.tk1, pac.tk1, pac.tk2, sc.tk1 ,bam.tk1)
ggbio::tracks(tks)


## ----sessionInfo--------------------------------------------------------------
sessionInfo()

