#  This function will return a list of differentially expressed genes
#' the function outputs a vector with gene symbols
#' @description For a pre-specified contrast this function will return the top differentially expressed genes
#' @param BigTab is a big tab object with GENENAME as the column with name of the gene
#' @param conts this is the label of the contrast that is prefixed with lgFCH, FCH, pvals , fdr or status
#' @examples

extract_degs <-function( BigTab, conts) {

  # vector of contrasts
  conts = paste0('_',conts)

  # create a list objects with differentially expressed genes for each contrast
  list.tab.degs = lapply(conts, function(var) {
    out = BigTab %>%
      dplyr::select(GENENAME , contains('Status')) %>%
      dplyr::select(GENENAME , status = contains(var)) %>%
      dplyr::filter(status != 0) %>% dplyr::select(GENENAME) %>% pull()
    return(out)
  })

    deg.genes = unique( do.call('c', list.tab.degs) )

    return(deg.genes)

  }
