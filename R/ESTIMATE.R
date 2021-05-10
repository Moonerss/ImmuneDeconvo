#' @name estimateScore
#' @title Calculation of stromal, immune, and ESTIMATE scores
#' @description This function computes stromal, immune, and ESTIMATE scores per sample using gene expression data.
#' @param merged_expr the expression of mergeed gene get from \code{filterCommonGenes}.
#' @param platform character string indicating platform type. Defaults to "affymetrix".
#' @details This method is based on single sample gene set enrichment analysis (ssGSEA) algorithm. This function computes stromal, immune, and ESTIMATE scores using gene-level expression data. For
#'   Affymetrix platform data, tumor purity are derived from ESTIMATE scores by applying non-linear
#'   squares methods to TCGA Affymetrix expression data (n=995).
#' @return
#'   Returns \emph{data.frame} with components:
#'   \itemize{
#'     \item StromalScorenumeric scalar specifying the presence of stromal cells in tumor tissue
#'     \item ImmuneScorenumeric scalar specifying the level of infiltrating immune cells in tumor tissue
#'     \item ESTIMATEScorenumeric scalar specifying tumor cellularity
#'     \item TumorPuritynumeric scalar specifying ESTIMATE-based tumor purity with value in range[0,1]
#'   }
#' @author Erjie Zhao <2055469819@qq.com>
#' @export
#' @examples
#' \dontrun{
#'    file <- system.file("extdata", "sample_input.txt", package="ESTIMATE")
#'    expression <- read.table(file, header = TRUE, stringsAsFactors = FALSE)
#'    merged_expr <- filterCommonGenes(expr = data.matrix(expression), id = "GeneSymbol")
#'    scores <- estimateScore(merged_expr, platform = "affymetrix")
#' }
estimateScore <- function(merged_expr, platform=c("affymetrix", "agilent", "illumina")) {
  ## Check arguments
  stopifnot(is.data.frame(merged_expr))
  platform <- match.arg(platform)

  ## input dataset
  ds <- merged_expr[,-1]
  rownames(ds) <- merged_expr[,1]

  descs <- ds[, 1]
  ds <- ds[-1]
  row.names <- row.names(ds)
  names <- names(ds)
  dataset <- list(ds=ds,
                  row.names=row.names,
                  descs=descs,
                  names=names)
  m <- data.matrix(dataset$ds)
  gene.names <- dataset$row.names
  sample.names <- dataset$names
  Ns <- length(m[1, ]) # Number of genes
  Ng <- length(m[, 1]) # Number of samples
  # temp <- strsplit(input.ds, split="/")
  # s <- length(temp[[1]])
  # input.file.name <- temp[[1]][s]
  # temp <- strsplit(input.file.name, split=".gct")
  # input.file.prefix <-  temp[[1]][1]

  ## Sample rank normalization
  for (j in 1:Ns) {
    m[, j] <- rank(m[, j], ties.method="average")
  }
  m <- 10000*m/Ng

  ## SI_geneset
  gs <- as.matrix(SI_geneset[, -1],dimnames=NULL)
  N.gs <- 2
  gs.names <- row.names(SI_geneset)

  ## Loop over gene sets
  score.matrix <- matrix(0, nrow=N.gs, ncol=Ns)
  for (gs.i in 1:N.gs) {
    gene.set <- gs[gs.i,]
    gene.overlap <- intersect(gene.set, gene.names)
    print(paste(gs.i, "gene set:", gs.names[gs.i], " overlap=", length(gene.overlap)))
    if (length(gene.overlap) == 0) {
      score.matrix[gs.i, ] <- rep(NA, Ns)
      next
    } else {
      ES.vector <- vector(length=Ns)

      ## Enrichment score
      for (S.index in 1:Ns) {
        gene.list <- order(m[, S.index], decreasing=TRUE)
        gene.set2 <- match(gene.overlap, gene.names)
        correl.vector <- m[gene.list, S.index]

        TAG <- sign(match(gene.list, gene.set2, nomatch=0))    # 1 (TAG) & 0 (no.TAG)
        no.TAG <- 1 - TAG
        N <- length(gene.list)
        Nh <- length(gene.set2)
        Nm <-  N - Nh
        correl.vector <- abs(correl.vector)^0.25
        sum.correl  <- sum(correl.vector[TAG == 1])
        P0 <- no.TAG / Nm
        F0 <- cumsum(P0)
        Pn <- TAG * correl.vector / sum.correl
        Fn <- cumsum(Pn)
        RES <- Fn - F0
        max.ES <- max(RES)
        min.ES <- min(RES)
        if (max.ES > - min.ES) {
          arg.ES <- which.max(RES)
        } else {
          arg.ES <- which.min(RES)
        }
        ES <- sum(RES)
        EnrichmentScore <- list(ES=ES,
                                arg.ES=arg.ES,
                                RES=RES,
                                indicator=TAG)
        ES.vector[S.index] <- EnrichmentScore$ES
      }
      score.matrix[gs.i, ] <- ES.vector
    }
  }

  score.data <- data.frame(score.matrix)
  names(score.data) <- sample.names
  row.names(score.data) <- gs.names
  estimate.score <- apply(score.data, 2, sum)

  if (platform != "affymetrix"){
    score.data <- rbind(score.data, estimate.score)
    rownames(score.data) <- c("StromalScore",
                              "ImmuneScore",
                              "ESTIMATEScore")
  } else {
    ##---------------------------------------------------------------------
    ## Calculate ESTIMATE-based tumor purity (Affymetrix-specific)
    convert_row_estimate_score_to_tumor_purity <- function(x) {
      stopifnot(is.numeric(x))
      cos(0.6049872018 + 0.0001467884*x)
    }

    est.new <- NULL
    for (i in 1:length(estimate.score)) {
      est_i <- convert_row_estimate_score_to_tumor_purity(estimate.score[i])
      est.new <- rbind(est.new, est_i)
      if (est_i >= 0) {
        next
      } else {
        message(paste(sample.names[i],": out of bounds", sep=""))
      }
    }
    colnames(est.new) <- c("TumorPurity")
    estimate.t1 <- cbind(estimate.score, est.new)
    x.bad.tumor.purities <- estimate.t1[, "TumorPurity"] < 0
    estimate.t1[x.bad.tumor.purities, "TumorPurity"] <- NA
    score.data <- rbind(score.data, t(estimate.t1))
    rownames(score.data) <- c("StromalScore",
                              "ImmuneScore",
                              "ESTIMATEScore",
                              "TumorPurity")
  }
  return(score.data)
}

#' @name filterCommonGenes
#' @title Intersect input data with 10,412 common genes
#' @description This function unifies different number of genes per platform against 10,412 common genes.
#' @param expr a matrix including gene-level expression data, rownames is genes, colnames is samples.
#' @param id character string indicating which gene identifier to use when matching.
#' @details The number of genes in expression data is different for each platform and this difference influences the computational results of stromal and immune scores. To compare stromal, immune and
#'   ESTIMATE scores across different platforms or calculate ESTIMATE-based tumor purity using
#'   Affymetrix expression data, users need to unify the gene identifiers of the input data against the
#'   common genes.
#' @return return a \emph{data.frame} of the result.
#' @author Erjie Zhao <2055469819@qq.com>
#' @export
#' @examples
#'   file <- system.file("extdata", "sample_input.txt", package="ESTIMATE")
#'   expression <- read.table(file, header = TRUE, stringsAsFactors = FALSE)
#'   merged_expr <- filterCommonGenes(expr = data.matrix(expression), id = "GeneSymbol")
filterCommonGenes <- function(expr, id=c("GeneSymbol", "EntrezID")) {
  ## Check arguments
  stopifnot(is.matrix(expr) && is.numeric(expr))
  id <- match.arg(id)

  expr <- as.data.frame(expr)
  merged.df <- merge(common_genes, expr, by.x=id, by.y="row.names")
  rownames(merged.df) <- merged.df$GeneSymbol
  merged.df <- merged.df[, -1:-ncol(common_genes)]
  print(sprintf("Merged dataset includes %d genes (%d mismatched).",
                nrow(merged.df),
                nrow(common_genes) - nrow(merged.df)))

  ## get the object
  res <- data.frame(NAME = rownames(merged.df), Description=rownames(merged.df), merged.df)
  rownames(res) <- NULL
  res <- res[order(res$NAME),]
  return(res)
}

#' @name plotPurity
#' @title Plot tumor purity
#' @description Plot tumor purity based on ESTIMATE score.
#' @param scores A data frame object get from \code{estimateScore}.
#' @param samples ector of character strings specifying sample names to be plotted. Defaults to "all_samples", which creates plots for all input samples
#' @param platform character string indicating platform type. Defaults to "affymetrix"
#' @details This function produces scatterplots for each requested sample; it plots tumor purity against ESTIMATE score. At present, only the Affymetrix platform is supported.
#' @return return an ggplot object
#' @import ggplot2 patchwork
#' @export
#' @author Erjie Zhao <2055469819@qq.com>
#' @examples
#'  \dontrun{
#'    file <- system.file("extdata", "sample_input.txt", package="ESTIMATE")
#'    expression <- read.table(file, header = TRUE, stringsAsFactors = FALSE)
#'    merged_expr <- filterCommonGenes(expr = expression, id = "GeneSymbol")
#'    scores <- estimateScore(merged_expr, platform = "affymetrix")
#'    plot_obj <- plotPurity(scores)
#'  }
plotPurity <- function(scores, samples="all_samples",
                       platform=c("affymetrix", "agilent", "illumina")) {
  options(warn = -1)
  ## Check arguments
  platform <- match.arg(platform)

  if (platform != "affymetrix"){
    stop("not implemented")
  }

  ## Begin processing

  ##-------------------------------------------------------------------------
  convert_row_estimate_score_to_tumor_purity <- function(x) {
    stopifnot(is.numeric(x))
    cos(0.6049872018 + 0.0001467884*x)
  }

  ## Read ESTIMATE data file
  estimate.df <- as.data.frame(t(scores))
  samplenames <- rownames(estimate.df)
  Affy.model <- PurityDataAffy
  pred.p <- Affy.model[, 5:7]
  est <- estimate.df[, 3]
  est.new <- estimate.df[, 4]

  ## ESTIMATE based tumor purity in scatterplot with prediction interval
  message("Plotting tumor purity based on ESTIMATE score")

  max.af <- max(Affy.model$ESTIMATEScore)
  min.af <- min(Affy.model$ESTIMATEScore)

  if (samples[1] == "all_samples"){
    Num.S <- nrow(estimate.df)
  } else {
    Num.S <- as.numeric(length(samples))
  }

  ## create plot object
  plot.obj <- list()
  for (i in 1:Num.S) {
    if(samples[1] =="all_samples"){
      samplename <- samplenames[i]
    } else {
      samplename <- samples[i]
    }

    geMin <- est[i] >= min.af
    leMax <- est[i] <= max.af
    withinMinMax <- geMin && leMax

    xlim <- if (!withinMinMax) {
      ## Expands plot boundary
      adjustment <- 500    # Arbitrary
      if (geMin) {
        from <- min.af
        to   <- est[i] + adjustment
      } else {
        from <- est[i] - adjustment
        to   <- max.af
      }
      c(from, to)
    } else {
      NULL
    }

    plot.obj[[i]] <- ggplot(Affy.model) +
      aes(x = Affy.model$ESTIMATEScore, y = Affy.model$tumor.purity) +
      geom_point(size = 1, colour = "lightgrey", fill = "white", shape = 21) +
      annotate("point", x = est[i], y = est.new[i], size = 3, colour = "black") +
      geom_line(mapping = aes(x = Affy.model$ESTIMATEScore, y = Affy.model$fit),
                colour = "darkgrey", linetype = 1, size = 1) +
      geom_line(mapping = aes(x = Affy.model$ESTIMATEScore, y = Affy.model$lwr.p),
                colour = "darkgrey", linetype = 2, size = 1) +
      geom_line(mapping = aes(x = Affy.model$ESTIMATEScore, y = Affy.model$upr.p),
                colour = "darkgrey", linetype = 2, size = 1) +
      geom_vline(xintercept = est[i], colour = "black", linetype = 2, size = 1) +
      geom_hline(yintercept = est.new[i], colour = "black", linetype = 2, size = 1) +
      labs(title = samplename, x = "ESTIMATE score", y = "Tumor purity") +
      ylim(0, 1) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5),
            panel.grid = element_blank())

    if (!withinMinMax) {
      plot.obj[[i]] <- plot.obj[[i]] +
        stat_function(fun = convert_row_estimate_score_to_tumor_purity,
                      n = 10000, colour = "grey", linetype = 1, size = 1) +
        xlim(from, to)
    }
  }

  ## merge all ggplot object
  res_obj <- patchwork::wrap_plots(plot.obj)

  return(res_obj)
}
