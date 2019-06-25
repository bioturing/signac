#' VeniceAllMarkers
#'
#' This function helps to find marker genes for all clusters in a Seurat v3 object.
#' We adapted some utilities functions in Seurat v3 package (Butler et al., Nature Biotechnology 2018) to handle the annotation of an
#' object.
#'
#' @param object An Seurat v3 object which contain at least one identification of clusters
#' @param ident.1 Name of the first indentification
#' @param ident.2 NULL if finding all markers, name of the second indentification otherwise
#' @param assay name of the assay to be used
#' @param slot name of the slot to be used
#' @param pvalue cut off p.value
#' @param node input a tree of clusters instead of a cluster name
#' @param verbose display the progress bar or not
#' @param only.pos filter the positive markers only
#'
#' @return
#' @export
#'
#' @examples
VeniceAllMarkers <- function(
  object,
  ident.1 = NULL,
  ident.2 = NULL,
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
          cluster = GetCellsIndices(data.use, object, idents.all[i], ident.2),
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
      any(as.character(x = ide5nt.2) %in% colnames(x = data.use))) {
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
  ### Label the cells before finding marker genes: 1: Ident.1,  2: Ident2, 3:Other
  clusters <- rep(0, ncol(data.use))
  clusters[match(ident.1, colnames(x = data.use))] <- 1
  clusters[match(ident.2, colnames(x = data.use))] <- 2
  if (sum(is.na(clusters)) > 0) {
    stop("Couldn't match some cell names with the colname of data.use")
  }
  return(clusters)
}
