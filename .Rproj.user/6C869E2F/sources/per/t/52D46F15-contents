#-----------------------------------------------
# correlation function
#-----------------------------------------------
runCor = function(varx,vary,dat = newdata){

  tempdata = dat %>%
    dplyr::select( SUBJID,
                   Y = all_of(varx),
                   X = all_of(vary))%>%
    distinct()%>%
    mutate(Y = numerize(Y)) %>%
    mutate(X = numerize(X)) %>%
    na.omit()

  cor.spearman = with(tempdata,cor.test(Y,X,method='spearman'))
  cor.pearson  = with(tempdata,cor.test(Y,X,method='pearson'))

  outx =  data.frame(
    X = varx,
    Y = vary,
    cor.spearman = cor.spearman$estimate,
    p.value.spearman = cor.spearman$p.value,
    cor.pearson = cor.pearson$estimate,
    p.value.pearson = cor.pearson$p.value,
    nobs = dim(tempdata)[1]
  ) %>% remove_rownames()


  return(outx)

}
