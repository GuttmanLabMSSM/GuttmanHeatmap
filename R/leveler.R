leveler<-function(fact){
  lev<-c()
  for (i in unique(fact)){
    lev<-c(lev,list(which(fact==i)))
  }
  names(lev)<-unique(fact)
  return(lev)}
