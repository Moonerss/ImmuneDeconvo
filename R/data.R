#' @name common_genes
#' @title 10,412 common genes.
#' @description This data set contains the information on 10,412 genes that are common among six platforms:
#'  \itemize{
#'    \item Affymetrix HG-U133Plus2.0
#'    \item Affymetrix HT-HG-U133A
#'    \item Affymetrix Human X3P
#'    \item Agilent 4x44K (G4112F)
#'    \item Agilent G4502A
#'    \item Illumina HiSeq RNA sequence
#'  }
#'  All gene information (according to NCBI) was current as of this package’s release date.
#' @format
#'  All gene information (according to NCBI) was current as of this package’s release date.
#'  \describe{
#'    \item{EntrezID}{integer scalar specifying Entrez Gene ID}
#'    \item{GeneSymbol}{character string specifying official HUGO gene symbol}
#'    \item{Synonyms}{character string specifying alternative gene symbols, each delimited by vertical bar}
#'    \item{GeneName}{character string specifying official HUGO gene name}
#'    \item{Chromosome}{character string specifying chromosomal location}
#'  }
#'
#' @author Erjie Zhao <2055469819@qq.com>
#' @usage
#'   data(common_genes)
"common_genes"

#' @name PurityDataAffy
#' @title Affymetrix data
#' @description
#'  This data set contains stromal, immune, and ESTIMATE scores in all Affymetrix expression data (n=995) that was used to develop the formula for predicting tumor purity based on raw estimate score.
#'  This data set also contains tumor purity based on ABSOLUTE algorithm, predicted tumor purity, and 95 ESTIMATE algorithm.
#' @usage data(PurityDataAffy)
#' @format
#' The object PurityDataAffy is a data.frame with components:
#' \describe{
#'    \item{tumor.purity}{numeric scalar specifying tumor purity calculated by ABSOLUTE algorithm}
#'    \item{StromalScore}{numeric scalar specifying the presence of stromal cells in tumor tissue}
#'    \item{ImmuneScore}{numeric scalar specifying the level of infiltrating immune cells in tumor tissue}
#'    \item{ESTIMATEScore}{numeric scalar specifying tumor cellularity}
#'    \item{fit}{numeric scalar specifying estimated tumor purity based on ESTIMATE algorithm}
#'    \item{lwr.p}{numeric scalar specifying 5 percent confidence interval}
#'    \item{upr.p}{numeric scalar specifying 95 percent confidence interval}
#' }
#' @author Erjie Zhao <2055469819@qq.com>
"PurityDataAffy"

#' @title two signatures for estimate
#' @description This data set contains two gene signatures (stromal and immune signatures). The stromal signature
#'   including 141 stroma-specific genes is designed to capture the presence of stroma in tumor tissue.
#'   The immune signature consisting of 141 immune cell-specific genes represents the infiltration of
#'   immune cells in tumor tissue.
#' @usage data(SI_geneset)
#' @format
#' The object SI_geneset is a data.frame with components:
#' \describe{
#'    \item{StromalSignature}{character string specifying 141 genes in stromal signature}
#'    \item{ImmuneSignature}{character string specifying 141 genes in immune signature}
#' }
#' @author Erjie Zhao <2055469819@qq.com>
"SI_geneset"


#' @title gene expression feature of cibersort analysis
#' @description gene expression matrix of 22 type immune cells,
#' it download from CIBERSORT (https://cibersort.stanford.edu/runcibersort.php).
#' @usage data(lm22)
#' @format
#' The object lm22 is a matrix with row is gene and column is  immune cell type.
#' @author Erjie Zhao <2055469819@qq.com>
"lm22"
