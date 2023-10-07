## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(message = FALSE, warning = FALSE)

## ----results='hide'-----------------------------------------------------------
library(vizAPA)
library(movAPA)

## ----eval=FALSE---------------------------------------------------------------
#  scPACds <- readPACds(pacFile = expma, colDataFile = coldata, noIntergenic = FALSE)

## ----eval=FALSE---------------------------------------------------------------
#  library(BSgenome.Hsapiens.NCBI.GRCh38)
#  bsgenome<-BSgenome.Hsapiens.NCBI.GRCh38
#  #check the consistency of chromosome names between PACdataset and BSgenome
#  isChrConsistent(scPACds,bsgenome)
#  #chromosome names in scPACds@anno$chr should be the same as in bsgenome$seqnames
#  scPACds@anno$chr<-gsub("chr","",scPACds@anno$chr)
#  scPACds =removePACdsIP(scPACds , bsgenome, returnBoth=FALSE, up=-140, dn=10, conA=6, sepA=NA,chrCheck = FALSE)
#  #merge pA sites within 24nt of each other
#  scPACds=mergePACds_v0(scPACds, d=24)
#  

## ----eval=FALSE---------------------------------------------------------------
#  filepath<- "./PBMC/scAPAtrap/"
#  
#  meta<-read.table(paste0(filepath,"scAPAtrap_cell.meta.txt"),sep = "\t",header = TRUE)
#  row.names(meta)<-meta$barcode
#  
#  #Check that cells with cell type annotations match those in scPACds@counts
#  dim(scPACds@counts)
#  length(intersect(meta$barcode,colnames(scPACds@counts)))
#  
#  count<- scPACds@counts
#  row.names(count)<- row.names(scPACds@anno)
#  count<- count[,intersect(colnames(count),meta$barcode)]
#  scPACds@counts<- count
#  scPACds@colData<- subset(scPACds@colData,scPACds@colData$group %in% colnames(count))
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
#  save(scPACds,file =paste0(filepath,"pbmc2_scAPAtrap_APA.rda"))

## -----------------------------------------------------------------------------
load("./PBMC/scAPAtrap/pbmc2_scAPAtrap_APA.rda")
summary(scPACds)

## ----eval=FALSE---------------------------------------------------------------
#  #paramter for what=updn/dn, specifying the upstream/downstream  region from the PAC.
#  #When up=-100, dn=100, then the output sequence is 201nt
#  fafiles = faFromPACds(scPACds, bsgenome, what = "updn", fapre = "pbmc2.scAPAtrap.200nt.GRCh38",
#                         up = -100, dn = 100, byGrp = "ftr", chrCheck = FALSE)
#  
#  #Plots single nucleotide profile using the sequences of 3'UTR
#  plotATCGforFAfile("pbmc2.scAPAtrap.200nt.GRCh38.3UTR.fa", ofreq = FALSE, opdf = TRUE, refPos = 101,
#                    filepre = "")

## ----out.height="70%",out.width="60%",echo=FALSE------------------------------
knitr::include_graphics("./PBMC/scAPAtrap/pbmc2.scAPAtrap.200nt.GRCh38.3UTR.png")

## -----------------------------------------------------------------------------
vizStats(scPACds, group='celltype', figType="dot", log=TRUE)+ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))

## -----------------------------------------------------------------------------
gene="ENST00000307839"
vizStats(scPACds, group='celltype',gene=gene, figType="box")+ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))

## -----------------------------------------------------------------------------
bam.files<-c("dedup_pbmc2.forward.sorted_Bcell.mini.sorted.bam",
             "dedup_pbmc2.forward.sorted_CD4+Tcell.mini.sorted.bam",
             "dedup_pbmc2.forward.sorted_CD14+monocyte.mini.sorted.bam",
             "dedup_pbmc2.forward.sorted_CD16+monocyte.mini.sorted.bam",
             "dedup_pbmc2.forward.sorted_CytotoxicTcell.mini.sorted.bam",
             "dedup_pbmc2.forward.sorted_Dendriticcell.mini.sorted.bam",
             "dedup_pbmc2.forward.sorted_Megakaryocyte.mini.sorted.bam",
             "dedup_pbmc2.forward.sorted_Naturalkillercell.mini.sorted.bam",
             "dedup_pbmc2.forward.sorted_Plasmacytoiddendriticcell.mini.sorted.bam")
bam.path<-"./PBMC/pbmc2.bam/"
bam.labels<-c("Bcell","CD4+Tcell","CD14+monocyte","CD16+monocyte",
              "CytotoxicTcell","Dendriticcell","Megakaryocyte","Naturalkillercell","Plasmacytoiddendriticcell")
bam.groups<-c("Bcell","CD4+Tcell","CD14+monocyte","CD16+monocyte",
              "CytotoxicTcell","Dendriticcell","Megakaryocyte","Naturalkillercell","Plasmacytoiddendriticcell")
bams<-readBAMFileNames(bam.files=bam.files, bam.path=bam.path, bam.labels = bam.labels, bam.groups = bam.groups)

## -----------------------------------------------------------------------------
APAds=get3UTRAPAds(scPACds)
top10=APAds@anno[order(Matrix::rowSums(APAds@counts), decreasing = TRUE)[1:10],]
head(top10[1:5,c(1:6,13)])
gene=top10$gene[8]
gene


## -----------------------------------------------------------------------------
#check the consistency of chromosome names among PACdataset, BAM files, and genome annotations.
isChrConsistent( scPACds,bams, annoSource["gff"],exact=FALSE)

# We can also get full chromosome names of each object.
getChrs(annoSource, which='gff')
getChrs(bams)


## ----results='hide'-----------------------------------------------------------
##Here we add a tot_tagnum to the PACds@anno and do log-transformation to show the tot_tagnum clearly.
scPACds@anno$tot_tagnum=Matrix::rowSums(scPACds@counts)
scPACds@anno$tot_tagnum=log2(scPACds@anno$tot_tagnum+1)

vizTracks(gene =gene, 
          PACds.list=list(pA=scPACds), 
          PA.show=c("pos", "tot_tagnum"),
          annoSource=annoSource,
          PA.columns="coord", PA.width=10, logPA=TRUE,
          space5=1000, space3=1000)

## ----results='hide',out.height="90%",out.width="90%"--------------------------

#Due to the large number of cell types, we can choose a few of them to plot the bam coverage. 
vizTheme=list(PA.col="#8d9fbb",PA.fill="#8d9fbb",bams.col=c("#8DD3C7","#FFFFB3","#BEBADA"),bams.fill=c("#8DD3C7","#FFFFB3","#BEBADA"))
vizTracks(gene=gene, 
           bams=bams[c(1,8,9),], PACds.list=list(pA=scPACds), PA.show=c("pos"),
           annoSource=annoSource,
           PA.columns="coord", PA.width=10,
           space5=1000, space3=1000,vizTheme =vizTheme )

vizTheme=list(PA.col="#8d9fbb",PA.fill="#8d9fbb",bams.col=c("#FB8072","#80B1D3","#FDB462"),bams.fill=c("#FB8072","#80B1D3","#FDB462"))
vizTracks(gene=gene, 
           bams=bams[c(2,3,4),], PACds.list=list(pA=scPACds), PA.show=c("pos"),
           annoSource=annoSource,
           PA.columns="coord", PA.width=10,
           space5=1000, space3=1000,vizTheme =vizTheme )

vizTheme=list(PA.col="#8d9fbb",PA.fill="#8d9fbb",bams.col=c("#B3DE69","#FCCDE5","#D9D9D9"),bams.fill=c("#B3DE69","#FCCDE5","#D9D9D9"))
vizTracks(gene=gene, 
           bams=bams[c(5,7,6),], PACds.list=list(pA=scPACds), PA.show=c("pos"),
           annoSource=annoSource,
           PA.columns="coord", PA.width=10,
           space5=1000, space3=1000,vizTheme =vizTheme )


## ----results='hide',out.height="90%",out.width="90%"--------------------------

vizTracks(gene=gene, 
          PACds.list=list(pA=scPACds), PA.show=c("pos"),
          cells=TRUE, cells.group='celltype',
          cells.method=c('sum'), cells.sort=c('group'),
          cells.width=100,
          annoSource=annoSource,
          PA.columns="coord", PA.width=10, logPA=TRUE,
          space5=1000, space3=1000,
          )  +ggplot2::guides(color=ggplot2::guide_legend(ncol = 2,
                            byrow = F,
                            reverse = T))


