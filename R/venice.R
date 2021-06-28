  #' VeniceAllMarkers
#'
#' This function helps to find marker genes for all clusters in a Seurat v3 object.
#' We adapted some utilities functions in Seurat v3 package (Butler et al., Nature Biotechnology 2018) to handle the annotation of an
#' object, extract the counts data from an object, and build the cluster tree.
#'
#' @param object An Seurat v3 object which contain at least one identification of clusters
#' @param assay name of the assay to be used
#' @param slot name of the slot to be used
#' @param pvalue cut off p.value
#' @param node input a tree of clusters instead of a cluster name
#' @param verbose display the progress bar or not
#' @param only.pos filter the positive markers only
#' @export
VeniceAllMarkers <- function(
  object,
  assay = "RNA",
  slot = 'data',
  pvalue = 0.05,
  node = NULL,
  verbose = TRUE,
  only.pos = FALSE,
  logfc.threshold = 0.25,
  nperm = 0
) {
  data.use <-  Seurat:::GetAssayData(object = object[[assay]], slot = slot)
  MapVals <- function(vec, from, to) {
    vec2 <- setNames(object = to, nm = from)[as.character(x = vec)]
    vec2[is.na(x = vec2)] <- vec[is.na(x = vec2)]
    return(unname(obj = vec2))
  }
  if (is.null(x = node)) {
    idents.all <- sort(x = unique(x = Seurat:::Idents(object = object)))
  } else {
    tree <- Seurat:::Tool(object = object, slot = 'BuildClusterTree')
    if (is.null(x = tree)) {
      stop("Please run 'BuildClusterTree' before finding markers on nodes")
    }
    descendants <- Seurat:::DFT(tree = tree, node = node, include.children = TRUE)
    all.children <- sort(x = tree$edge[, 2][!tree$edge[, 2] %in% tree$edge[, 1]])
    descendants <- MapVals(
      vec = descendants,
      from = all.children,
      to = tree$tip.label
    )
    drop.children <- setdiff(x = tree$tip.label, y = descendants)
    keep.children <- setdiff(x = tree$tip.label, y = drop.children)
    orig.nodes <- c(
      node,
      as.numeric(x = setdiff(x = descendants, y = keep.children))
    )
    tree <- ape::drop.tip(phy = tree, tip = drop.children)
    new.nodes <- unique(x = tree$edge[, 1, drop = TRUE])
    idents.all <- (tree$Nnode + 2):max(tree$edge)
  }
  genes.de <- list()
  messages <- list()
  for (i in 1:length(x = idents.all)) {
    if (verbose) {
      message("Calculating cluster ", idents.all[i])
    }
    genes.de[[i]] <- tryCatch(
      expr = {
        Signac::VeniceMarker(
          S4_mtx  = data.use,
          cluster = GetCellsIndices(data.use, object, idents.all[i],  NULL),
          perm = nperm,
          verbose = verbose
        )
      },
      error = function(cond) {
        return(cond$message)
      }
    )
    if (class(x = genes.de[[i]]) == "character") {
      messages[[i]] <- genes.de[[i]]
      genes.de[[i]] <- NULL
    }
  }
  gde.all <- data.frame()
  for (i in 1:length(x = idents.all)) {
    if (is.null(x = unlist(x = genes.de[i]))) {
      next
    }
    gde <- genes.de[[i]]
    if (nrow(x = gde) > 0) {
      gde <- gde[order(gde[["Log10.adjusted.p.value"]], -gde[["Log2.fold.change"]]), ]
      gde <- subset(x = gde, subset = Log10.adjusted.p.value < log10(pvalue))
      gde$cluster <- idents.all[i]
      gde.all <- rbind(gde.all, gde)
      }
  }
  if ((only.pos) && nrow(x = gde.all) > 0) {
    gde.all <- subset(x = gde.all, subset = gde.all[["Log2.fold.change"]] > 0)
  }
  rownames(x = gde.all) <- make.unique(names = as.character(x = gde.all[["Gene.Name"]]))
  if (nrow(x = gde.all) == 0) {
    warning("No DE genes identified", call. = FALSE, immediate. = TRUE)
  }
  if (length(x = messages) > 0) {
    warning("The following tests were not performed: ", call. = FALSE, immediate. = TRUE)
    for (i in 1:length(x = messages)) {
      if (!is.null(x = messages[[i]])) {
        warning("When testing ", idents.all[i], " versus all:\n\t", messages[[i]], call. = FALSE, immediate. = TRUE)
      }
    }
  }
  if (!is.null(x = node)) {
    gde.all$cluster <- MapVals(
      vec = gde.all$cluster,
      from = new.nodes,
      to = orig.nodes
    )
  }
  gde.all <- gde.all[, c("Gene.Name", "Log2.fold.change", "Log10.adjusted.p.value",
                         "Log10.p.value", "cluster",
                         "Dissimilarity", "Bin.count",
                         "Up.Down.score")]
  gde.all <- subset(gde.all, abs(Log2.fold.change) > logfc.threshold)
  return(gde.all)
}

#' This function helps to find marker genes for two clusters in a Seurat v3 object.
#' We adapted some utilities functions in Seurat v3 package (Butler et al., Nature Biotechnology 2018) to handle the annotation of an
#' object, extract the counts data from an object, and build the cluster tree.
#' Some parameters are derived from the FindMarkers function of Seurat package v3
#' https://github.com/satijalab/seurat/blob/master/R/differential_expression.R
#' @param ident.1 Identity class to define markers for; pass an object of class
#' \code{phylo} or 'clustertree' to find markers for a node in a cluster tree;
#' passing 'clustertree' requires \code{\link{BuildClusterTree}} to have been run
#' @param ident.2 A second identity class for comparison; if \code{NULL},
#' use all other cells for comparison; if an object of class \code{phylo} or
#' 'clustertree' is passed to \code{ident.1}, must pass a node to find markers for
#' @param group.by Regroup cells into a different identity class prior to performing differential expression (see example)
#' @param subset.ident Subset a particular identity class prior to regrouping. Only relevant if group.by is set (see example)
#' @param assay Assay to use in differential expression testing
#' @param slot Slot to pull data from; note that if \code{test.use} is "negbinom", "poisson", or "DESeq2",
#' \code{slot} will be set to "counts"
#'
#' @importFrom methods is
#'
#' @rdname VeniceFindMarkers
#' @export
#'
VeniceFindMarkers <- function(
  object,
  ident.1 = NULL,
  ident.2 = NULL,
  group.by = NULL,
  subset.ident = NULL,
  assay = NULL,
  slot = 'data',
  reduction = NULL,
  features = NULL,
  logfc.threshold = 0.25,
  verbose = TRUE,
  pvalue = 0.05,
  only.pos = FALSE,
  nperm = 0
) {
  if (!is.null(x = group.by)) {
    if (!is.null(x = subset.ident)) {
      object <- subset(x = object, idents = subset.ident)
    }
    Idents(object = object) <- group.by
  }
  if (!is.null(x = assay) && !is.null(x = reduction)) {
    stop("Please only specify either assay or reduction.")
  }
  assay <- if (is.null(assay)) DefaultAssay(object = object) else assay
  data.use <-  GetAssayData(object = object[[assay]], slot = slot)

  counts <- switch(
    EXPR = slot,
    'scale.data' = Seurat:::GetAssayData(object = object[[assay]], slot = "counts"),
    numeric()
  )

  indices <- GetCellsIndices(data.use, object, ident.1, ident.2)
  de.results <- Signac::VeniceMarker(
    S4_mtx  = data.use,
    cluster = indices,
    perm = nperm,
    verbose = verbose
  )
  de.results <- de.results[, c("Gene.Name", "Log2.fold.change", "Log10.adjusted.p.value",
                         "Log10.p.value", "Dissimilarity", "Bin.count", "Up.Down.score")]
  de.results <- subset(de.results, abs(Log2.fold.change) > logfc.threshold)
  return(de.results)
}

#' This function helps to extract the cell indices given two indents
#' We adapted some utilities functions in Seurat v3 package (Butler et al., Nature Biotechnology 2018) to handle the annotation of an
#' object
#' This function are edited from the FindMarkers function of Seurat package v3
#' https://github.com/satijalab/seurat/blob/master/R/differential_expression.R
#'
#' @param data.use The sparse matrix object that contains the counts data
#' @param object The Seurat v3 object
#' @param ident.1 Identity class to define markers for; pass an object of class
#' \code{phylo} or 'clustertree' to find markers for a node in a cluster tree;
#' passing 'clustertree' requires \code{\link{BuildClusterTree}} to have been run
#' @param ident.2 A second identity class for comparison; if \code{NULL},
#' use all other cells for comparison; if an object of class \code{phylo} or
#' 'clustertree' is passed to \code{ident.1}, must pass a node to find markers for
GetCellsIndices <- function(data.use, object, ident.1, ident.2) {
  if (is.null(x = ident.1)) {
    stop("Please provide ident.1")
  } else if ((length(x = ident.1) == 1 && ident.1[1] == 'clustertree') || is(object = ident.1, class2 = 'phylo')) {
    if (is.null(x = ident.2)) {
      stop("Please pass a node to 'ident.2' to run FindMarkers on a tree")
    }
    tree <- if (is(object = ident.1, class2 = 'phylo')) {
      ident.1
    } else {
      Seurat:::Tool(object = object, slot = 'BuildClusterTree')
    }
    if (is.null(x = tree)) {
      stop("Please run 'BuildClusterTree' or pass an object of class 'phylo' as 'ident.1'")
    }
    ident.1 <- tree$tip.label[Seurat:::GetLeftDescendants(tree = tree, node = ident.2)]
    ident.2 <- tree$tip.label[Seurat:::GetRightDescendants(tree = tree, node = ident.2)]
  }
  if (length(x = as.vector(x = ident.1)) > 1 &&
      any(as.character(x = ident.1) %in% colnames(x = data.use))) {
    bad.cells <- colnames(x = data.use)[which(x = !as.character(x = ident.1) %in% colnames(x = data.use))]
    if (length(x = bad.cells) > 0) {
      stop(paste0("The following cell names provided to ident.1 are not present in the object: ", paste(bad.cells, collapse = ", ")))
    }
  } else {
    ident.1 <- Seurat:::WhichCells(object = object, idents = ident.1)
  }
  # if NULL for ident.2, use all other cells
  if (length(x = as.vector(x = ident.2)) > 1 &&
      any(as.character(x = ident.2) %in% colnames(x = data.use))) {
    bad.cells <- colnames(x = data.use)[which(!as.character(x = ident.2) %in% colnames(x = data.use))]
    if (length(x = bad.cells) > 0) {
      stop(paste0("The following cell names provided to ident.2 are not present in the object: ", paste(bad.cells, collapse = ", ")))
    }
  } else {
    if (is.null(x = ident.2)) {
      ident.2 <- setdiff(x = colnames(x = data.use), y = ident.1)
    } else {
      ident.2 <- Seurat:::WhichCells(object = object, idents = ident.2)
    }
  }
  ### Label the cells before finding marker genes: 1: Ident.1,  0: Ident2, >1:Other
  clusters <- rep(2, ncol(data.use))
  clusters[match(ident.1, colnames(x = data.use))] <- 1
  clusters[match(ident.2, colnames(x = data.use))] <- 0
  if (sum(is.na(clusters)) > 0) {
    stop("Couldn't match some cell names with the colname of data.use")
  }
  return(clusters)
}
