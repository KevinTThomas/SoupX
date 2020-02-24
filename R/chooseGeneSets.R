#' Selects sets of genes from a preselected \code{cadidateList} to estimate the contamination fraction of a soupChannel object. 
#'
#' The purpose of this function is to identify gene sets and specific genes in the candidate list that are also good candidates for estimating soup contamination.
#'
#' The curated sets can be compared and passed to \code{\link{calculateContaminatingFraction}} to estimate the Soup contamination in each channel.
#'
#' @export
#' @param sc A SoupChannel object.
#' @param inferredGenes A gene table passed from \code{\link{inferNonExpressedGenes}}. If \code{NULL}, the function is run on the SoupChannel object. 
#' @param candidateList A list of pre-selected candidate gene sets to cross-reference with the inferred non-expressed genes.
#' @param ... Passed to \code{\link{inferNonExpressedGenes}}.
#' @return A dataframe with each of the tested gene sets and specific genes. They are listed in order by the estimated global contamination fraction and the summed extremity score for all the genes in the set.
#' @importFrom utils head
#' @importFrom BiocGenerics rownames, unlist, sapply, lapply, mapply
#' @importFrom dplyr intersect
#' @importFrom rlang set_names
chooseGeneSets = function (sc, inferredGenes = NULL, candidateList, top_n = 20,...) {
  if(is.null(inferredGenes)) {
    print("Inferring non-expressed genes.")
    inferredGenes = inferNonExpressedGenes(sc,...)
  }
  topGenes = inferredGenes %>% head(top_n)
  NEgns = topGenes %>%
    row.names() %>%
    intersect(y = unlist(candidateList)) %>%
    (function(y) lapply(candidateList, function (x) intersect(x, y))) %>%
    .[lapply(.,length)>0]
  flagged_sets = intersect(names(NEgns), names(candidateList))
  
  if (is_empty(flagged_sets)) {
    NEgns = rownames(topGenes) %>% list(.) %>% `names<-`(rownames(topGenes))
    candidateList = NEgns
    flagged_sets = intersect(names(NEgns), names(candidateList))
  }
  
  # nullMat = sapply(flagged_sets,
  #                  function (x) estimateNonExpressingCells(sc, candidateList[x]),
  #                  simplify = FALSE
  # ) %>%
  #   lapply(function (x) rownames_to_column(data.frame(x), var = "Barcode")) %>%
  #   (function(x) {if(length(x)>1) {reduce(x, merge)} else {x[[1]]}}) %>%
  #   `rownames<-`(.$Barcode) %>%
  #   .[, which(colnames(.) != "Barcode"), drop=FALSE] %>%
  #   data.matrix(rownames.force = TRUE)
  # 
  # globalRhos = character(length = length(flagged_sets)) %>% set_names(nm = flagged_sets)
  # globalRhosLow = character(length = length(flagged_sets)) %>% set_names(nm = flagged_sets)
  # globalRhosHigh = character(length = length(flagged_sets)) %>% set_names(nm = flagged_sets)
  # for (e in flagged_sets) {
  #   if(sum(nullMat[,e])>0){
  #     tmp = tryCatch(
  #       {
  #       message(paste0("Calculating contamination fraction of ", e, " gene set."))
  #       suppressMessages(calculateContaminationFraction(sc,candidateList[e],nullMat[,e,drop=FALSE], verbose = FALSE))
  #       },
  #       error = function (cond) {
  #         message("Can't calculate contamination fraction:")
  #         message(cond)
  #         sct = sc
  #         sct$metaData$rho = character(length = dim(sct$toc)[2])
  #         sct$metaData$rhoLow = character(length = dim(sct$toc)[2])
  #         sct$metaData$rhoHigh = character(length = dim(sct$toc)[2])
  #         return(sct)
  #       }
  #     )
  #     globalRhos[e] = tmp$metaData$rho[1]
  #     globalRhosLow[e] = tmp$metaData$rhoLow[1]
  #     globalRhosHigh[e] = tmp$metaData$rhoHigh[1]
  #   }else{
  #     globalRhos[e] = NA
  #     globalRhosLow[e] = NA
  #     globalRhosHigh[e] = NA
  #   }
  # }
  
  tab = data.frame(genes = sapply(NEgns, function(x) str_flatten(x, collapse = ", ") %>% str_remove_all(., "hg38_")),
                   pct.set = mapply(function (x,y) length(x)/length(y), x = NEgns, y = candidateList[flagged_sets]),
                   sum.extremity = sapply(NEgns,  function (x) sum(topGenes$extremity[rownames(topGenes) %in% x])),
                   frac.is.Useful = sapply(NEgns, function (x) {sum(topGenes$isUseful[rownames(topGenes) %in% x])/length(x)}),
                   # est.global.contam = globalRhos,
                   # est.global.lo = globalRhosLow,
                   # est.global.hi = globalRhosHigh,
                   check.rows = TRUE
  )
  # tab$est.global.contam = as.numeric(as.character(tab$est.global.contam))
  # tab$est.global.lo = as.numeric(as.character(tab$est.global.lo))
  # tab$est.global.hi = as.numeric(as.character(tab$est.global.hi))
  tab = tab[order(tab$sum.extremity, decreasing = TRUE),]
  
  return(tab)
}