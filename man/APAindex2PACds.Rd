% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vizAPA_funclib.R
\name{APAindex2PACds}
\alias{APAindex2PACds}
\title{Convert an APA index table to a PACdataset}
\usage{
APAindex2PACds(APAindex, colData = NULL)
}
\arguments{
\item{APAindex}{A table with rownames denoting genes. Columns are samples, and values are the APA index values. If row names is NULL, then will be set as G1..Gn.}

\item{colData}{The sample information for the APAindex's columns. If is NULL, then will be automatically set as group1 for all samples.}
}
\value{
A PACdataset with counts slot the APAindex table. However, this PACdataset is not the real PACdataset, which does not have chr/strand/coord columns.
It is just used in vizAPA for vizAPAmarkers.
}
\description{
APAindex2PACds converts an APA index table to a PACdataset.
The APA index can be obtained by movAPA::movAPAindex, which is a table with each row denoting a gene and each column a sample (e.g., cell).
The element in the APA index table is the APA index for each APA gene, e.g., RUD, denoting the relative usage of distal pAs.
}
\examples{
#Calculate the smartRUD index for each gene.
#Only genes with 3'UTR APA can be used for RUD calculation.
data(scPACds)
pd=movAPA::get3UTRAPApd(pacds=scPACds, minDist=50, maxDist=5000,
                        minRatio=0.05, fixDistal=FALSE, addCols='pd')
srud=movAPA::movAPAindex(pd, method="smartRUD", sRUD.oweight=FALSE)
PACds=APAindex2PACds(srud, colData=pd@colData)
movAPA::summary(PACds)
}
\seealso{
Other vizMarkers functions: 
\code{\link{getAPAindexPACds}()},
\code{\link{getAPAmarkers}()},
\code{\link{vizAPAMarkers}()}
}
\concept{vizMarkers functions}
