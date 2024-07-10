extractImprovement = function(PRE.LS,PRE.NL,POST.LS,genes,coefs){

  d = coefs %>%
    data.frame()  %>%
    rownames_to_column('GENENAME') %>%
    dplyr::filter(GENENAME %in% genes) %>%
    dplyr::select(GENENAME,bl.ls = all_of(PRE.LS) , bl.nl = all_of(PRE.NL) , post.ls = all_of(POST.LS))%>%
    mutate(Regulation = bl.ls - bl.nl) %>%
    mutate(Reg.Type = ifelse(Regulation<0,'Down','Up'))%>%
    mutate(Modulation = post.ls - bl.ls ) %>%
    mutate(Mod.Type = ifelse(Modulation<0,'Down','Up'))%>%
    mutate(Normalization = post.ls-bl.nl)%>%
    mutate(post.expr = post.ls)%>%
    dplyr::filter(GENENAME %in% genes) %>%
    mutate(Improvement =  -100*(Modulation/Regulation))

  regData = d%>%
    dplyr::select(GENENAME,Improvement,Regulation,Reg.Type,Modulation,Mod.Type,Normalization,post.expr) %>%
    mutate(PRE.LS = PRE.LS, PRE.NL=PRE.NL, POST.LS = POST.LS)


  return(regData)
}
