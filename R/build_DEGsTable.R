#  This function will build a table with the number of differentially expressed genes
#' the function outputs a table with the frequency of down and up regulated genes within each contrast
#' @description For a pre-specified contrast this function will build a frequency table for the number of DEGs
#' @param BigTab is a big tab object with GENENAME as the column with name of the gene
#' @param FCH is the threshold for fold-changes
#' @param pcut is a threshold for the p-value
#' @param cont.list is a list of contrasts
#' @param method is the method for type-I error control (none, fdr, bonferroni)
#' @examples


build_DEGsTable <-function(BigTab, FCH, pcut, cont.list , method = 'none') {

    cont.list = paste0('_', cont.list)


    list.tab.degs = lapply(cont.list, function(var) {


      out = BigTab %>% dplyr::select(contains('Status')) %>%
        dplyr::select(contains(var)) %>%
        pull() %>%
        tabyl() %>%
        dplyr::rename('Regulation' = '.')

      if (dim(out)[1] == 1) {
        out = rbind.data.frame(out, cbind.data.frame(
          Regulation = c('Down', 'Up'),
          n = c(0, 0),
          percent = c(0, 0)
        ))
      } else{
        out = out
      }

      out = out %>%
        mutate(Regulation = plyr::mapvalues(
          Regulation,
          from = c('-1', '0', '1'),
          to = c('Down', 'None', 'Up')
        )) %>%
        mutate(contrast = var) %>%
        dplyr::select(-percent) %>%
        reshape2::melt(., id.vars = c('Regulation', 'contrast')) %>%
        dplyr::filter(Regulation != 'None') %>%
        dplyr::select(-variable)
      return(out)
    })

    tab.degs = do.call(rbind, list.tab.degs) %>%
      reshape2::dcast(., contrast ~ Regulation, value.var = 'value') %>%
      mutate(FCH = FCH)

    if (method == 'fdr') {
      tab.degs = tab.degs %>% mutate(FDR = pcut)
    }
    if (method == 'none') {
      tab.degs = tab.degs %>% mutate(PVAL = pcut)
    }

    tab.degs = tab.degs %>%
      mutate(contrast = factor(contrast, levels = cont.list)) %>% arrange(contrast) %>%
      mutate(contrast = gsub('_', '', contrast))


    return(tab.degs)


    }
