## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(message = FALSE,warning = FALSE)

## ----results='hide'-----------------------------------------------------------
library(vizAPA)
library(movAPA)

## ----eval=FALSE---------------------------------------------------------------
#  filepath <- "./PBMC/sierra/"
#  peakfile <- paste0(filepath,"pbmc2_peak_annotations.txt")
#  countfile <-paste0(filepath,"pbmc2_count")
#  

## ----eval=FALSE---------------------------------------------------------------
#  peak_annotations <- read.delim(peakfile)
#  head(peak_annotations[1:3,])
#  
#  #Read in peak data saved in MEX format. Files can be in a gzipped (.gz) format.
#  count <- Sierra::ReadPeakCounts(data.dir=countfile,mm.file="matrix.mtx.gz",barcodes.file ="barcodes.tsv.gz")
#  count <- as.data.frame(count)
#  head(count[1:3,1:4])
#  
#  #Remove the numbers such as -1 in barcode
#  colnames(count)<-gsub("-1","",colnames(count))

## ----eval=FALSE---------------------------------------------------------------
#  #Keep the same peaks in count as in peak.annotation
#  infor <- intersect(row.names(count),row.names(peak_annotations))
#  peak_annotations<-peak_annotations[infor,]
#  count<-count[infor,]
#  
#  peak_annotations$coord <- 0
#  #When strand="+", the end position of the identified peak is used as the coordinate of the pA site
#  peak_annotations[peak_annotations$strand == "+",]$coord <- peak_annotations[peak_annotations$strand == "+",]$end
#  #When strand="-",the starting position of the identified peak is used as the coordinate of the pA site
#  peak_annotations[peak_annotations$strand == "-",]$coord <- peak_annotations[peak_annotations$strand == "-",]$start
#  
#  anno<- peak_annotations[,c("seqnames","start","end","strand","coord","gene_id")]
#  colnames(anno) <- c("chr","start","end","strand","coord","gene_id")
#  head(anno[1:3,1:6])
#  
#  #Remove duplicate peak data based on the annotation information of peak
#  anno <- unique(anno)
#  
#  #Preserve PA sites with gene name
#  anno<- anno[!is.na(anno$gene_id),]
#  anno$gene_id <- gsub("[.].*","",anno$gene_id)
#  count<-count[rownames(anno),]

## ----eval=FALSE---------------------------------------------------------------
#  #Retention of cells with annotated cell type information
#  meta<-read.table(paste0(filepath,"sierra_cell.meta.txt"),sep = "\t",header = TRUE)
#  row.names(meta)<-meta$barcode
#  dim(count)
#  length(intersect(meta$barcode,colnames(count)))
#  count<- count[,intersect(meta$barcode,colnames(count))]
#  #combine the result of the annotated peak and peak count
#  pafile<-cbind(anno,count)
#  head(pafile[,1:8])
#  
#  coldata<-data.frame(barcode=colnames(pafile)[7:ncol(pafile)],row.names =colnames(pafile)[7:ncol(pafile)] )
#  scPACds <- readPACds(pacFile = pafile, colDataFile = coldata,noIntergenic = FALSE)
#  

## ----eval=FALSE---------------------------------------------------------------
#  library(BSgenome.Hsapiens.NCBI.GRCh38)
#  bsgenome<-BSgenome.Hsapiens.NCBI.GRCh38
#  #check the consistency of chromosome names between PACdataset and BSgenome
#  isChrConsistent(scPACds,bsgenome)
#  
#  scPACds =removePACdsIP(scPACds , bsgenome, returnBoth=FALSE, up=-140, dn=10, conA=6, sepA=NA,chrCheck = FALSE)
#  #merge pA sites within 24nt of each other
#  scPACds<-mergePACds_v0(scPACds, d=24)
#  
#  dim(scPACds@counts)

## ----eval=FALSE---------------------------------------------------------------
#  loc<- match(scPACds@colData$group,meta$barcode)
#  scPACds@colData$celltype<-meta$celltype[loc]

## -----------------------------------------------------------------------------
annoSource=new("annoHub")
gff<-useGff("./PBMC/Homo_sapiens.GRCh38.109.gtf.gz")
annoSource=addAnno(annoSource, gff)

## ----eval=FALSE---------------------------------------------------------------
#  #check the consistency of chromosome names between PACdataset and genome annotations
#  isChrConsistent(scPACds,gff)
#  scPACds<-annotatePAC(scPACds,aGFF=gff)
#  
#  # extend 3UTR by 1000bp
#  scPACds=ext3UTRPACds(scPACds, 1000)
#  save(scPACds,file =paste0(filepath,"pbmc2_sierra_APA.rda"))

## -----------------------------------------------------------------------------
#Load annotated PAC data
load("./PBMC/sierra/pbmc2_sierra_APA.rda")
summary(scPACds)

## ----eval=FALSE---------------------------------------------------------------
#  fafiles = faFromPACds(scPACds, bsgenome, what = "updn", fapre = "pbmc2.sierra.200nt.GRCh38",
#                         up = -100, dn = 100, byGrp = "ftr", chrCheck = FALSE)
#  
#  plotATCGforFAfile("pbmc2.sierra.200nt.GRCh38.3UTR.fa", ofreq = FALSE, opdf = TRUE, refPos = 101,
#                    filepre = "")

## ----out.height="70%",out.width="60%",echo=FALSE------------------------------
knitr::include_graphics("C:/vizAPA/vizAPA data/PBMC/sierra/pbmc2.sierra.200nt.GRCh38.3UTR.png")

## ----results='hide',out.height="90%",out.width="90%"--------------------------
## get embeddings using the pA count matrix with normalization
scPACds=reduceDim(scPACds, nfeature=2000, dims=1:10, dimLabel='umap_norm', norm=TRUE)
vizUMAP(scPACds, group='celltype', xcol='umap_norm1', ycol='umap_norm2' )


## -----------------------------------------------------------------------------
APAds=get3UTRAPAds(scPACds)
top10=APAds@anno[order(Matrix::rowSums(APAds@counts), decreasing = TRUE)[1:10],]
gene=top10$gene[1]
gene

## ----results='hide'-----------------------------------------------------------

statTheme=setStatTheme(list(scale.high.col="purple",scale.low.col="grey"),check = TRUE)
vizUMAP(scPACds, group='celltype', xcol='umap_norm1', ycol='umap_norm2' ,genes = gene,statTheme = statTheme)

## ----results='hide'-----------------------------------------------------------
iPACds=getAPAindexPACds(scPACds, choose2PA='PD')
markers=getAPAmarkers(iPACds,  group='celltype',cluster1 = "CD4+Tcell" ,logFC = 0.25)

## -----------------------------------------------------------------------------
table(markers$cluster1, markers$cluster2)

## ----results='hide'-----------------------------------------------------------
markers=getAPAmarkers(iPACds,  group='celltype',cluster1 = "Dendritic cell" ,logFC = 0.25)


## -----------------------------------------------------------------------------
table(markers$cluster1, markers$cluster2)

## -----------------------------------------------------------------------------
vizAPAMarkers(iPACds, group='celltype', markers=markers$rowid[1:5], figType = 'bubble', statTheme=list(xgroup=FALSE))

## -----------------------------------------------------------------------------
vizAPAMarkers(iPACds, group='celltype', markers=markers$rowid[1:5], figType = 'violin')

