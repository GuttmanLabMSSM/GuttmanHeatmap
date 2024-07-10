# Generate ebfit output for mixed effects model.

runEbfit = function(tb){
  coefs = reshape2::dcast(tb,biomarker~contrasts,value.var = 'estimate')%>%
    column_to_rownames('biomarker')%>% as.matrix()

  pvals = reshape2::dcast(tb,biomarker~contrasts,value.var = 'p.value')%>%
    column_to_rownames('biomarker')%>% as.matrix()

  return(list(coef = coefs , p.value = pvals))
}
