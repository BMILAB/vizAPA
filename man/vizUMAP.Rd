% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vizAPA_funclib.R
\name{vizUMAP}
\alias{vizUMAP}
\title{Plot a UMAP figure for single-cell PACdataset}
\usage{
vizUMAP(
  PACds,
  group,
  xcol,
  ycol,
  annoUMAP = TRUE,
  selGroups = NULL,
  genes = NULL,
  PAs = NULL,
  log = FALSE,
  statTheme = NULL
)
}
\arguments{
\item{PACds}{A dataset of the PACdataset class.}

\item{group}{The column name recording cell categories in PACds@colData.}

\item{xcol}{the column name in PACds@colData denoting the X-axis coordinates of the 2D embedding.}

\item{ycol}{Same as xcol but for the Y-axis.  Both xcol and ycol are required to be in colData. See reduceDim() for getting 2D embeddings.}

\item{annoUMAP}{Whether to plot the UMAP only with cell annotation but no expression overlay.}

\item{selGroups}{A vector to specify the order and categories to plot, otherwise plot all categories.}

\item{genes}{A vector of gene names present in PACds@anno for filtering PACds.}

\item{PAs}{A vector of pA ids present in rownames of PACds@counts.
If multiple genes or PAs are given, the average is calculated before drawing.
If both gene and PAs are not provided, the average of all counts is obtained.}

\item{log}{TRUE to do log2 scale for the data, which is only applicable when PACds is count type.}

\item{statTheme}{A statTHEME-like list to specify the theme of the plot.}
}
\value{
One (annoUMAP=FALSE) or two (annoUMAP=TRUE) UMAP plots.
}
\description{
vizUMAP plots a UMAP plot where each point is a cell which is positioned based on the cell embeddings determined by the reduction technique.
}
\examples{
\dontrun{
# First, we check the coordinate labels of the 2D-embedding.
# For this data, the labels are UMAP_1 and UMAP_2.
colnames(PACds@colData)
# Plot the UMAP plot showing cell clusters and another UMAP plot
# overlaying with the mean expression value of PAs in each cell.
vizUMAP(PACds, group='celltype', xcol='UMAP_1', ycol='UMAP_2')
# Plot only the overlaying UMAP
vizUMAP(PACds, group='celltype', annoUMAP=FALSE,
       xcol='UMAP_1', ycol='UMAP_2')
# Providing a gene id or a list of genes in the gene column of the PACds,
# we can plot a UMAP overlaying with the mean expression value of the gene(s).
vizUMAP(PACds, group='celltype', xcol='UMAP_1', ycol='UMAP_2', genes=gene)
vizUMAP(PACds, group='celltype', xcol='UMAP_1', ycol='UMAP_2',
       selGroups=c('SC','RS'), PAs=PAids)
## get embeddings using the pA count matrix with normalization
PACds=reduceDim(PACds, nfeature=2000, dims=1:10,
               dimLabel='umap_norm', norm=TRUE)
## without normalization
PACds=reduceDim(PACds, nfeature=2000, dims=1:10,
               dimLabel='umap_raw', norm=FALSE)
## There are new columns adding to the colData slot.
head(PACds@colData)
vizUMAP(PACds, group='celltype', xcol='umap_raw1', ycol='umap_raw2')
vizUMAP(PACds, group='celltype',
       xcol='umap_norm1', ycol='umap_norm2')
vizUMAP(PACds, group='celltype', xcol='UMAP_1', ycol='UMAP_2')
}
}
\seealso{
Other vizUMAP functions: 
\code{\link{reduceDim}()}
}
\concept{vizUMAP functions}
