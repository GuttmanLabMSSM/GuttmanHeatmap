
doTab<-function(fdrsx,coefsx){
  require(stringr)
  require(weights)

  p.cuts.stars=c(0.001,0.01, 0.05, 0.1,1)
  p.symbols.stars= c('***',"**", "*", "+","  ")

  coefsx[coefsx==0]<-1e-10

  fdrsx[] <- sapply(fdrsx,starmaker, p.levels =p.cuts.stars, symbols = p.symbols.stars)

  transformfch <- function(lgfch) {
    fch <- sign(lgfch) * (2^abs(lgfch))
    return(fch)
  }

  padding<-function(x,n){
    m=ceiling(log10(abs(x)+.001))
    x<-as.character(x)
    initspace=n-m
    if (substr(x,1,1)!='-'){
      initspace<-initspace+1
    }
    if (!str_detect(x,'\\.')){
      x<-paste(x,'.',sep='')
    }
    while (nchar(x)<n-initspace+4){
      x<-paste(x,'0',sep='')
    }
    initstr=paste0("<span style='color:white'>",str_c(rep('\\-',initspace),collapse=''))
    x<-paste(initstr,"</span>",x,sep='')
    return(x)
  }



  padstars<-function(x){
    while (nchar(x)<3){x<-paste(x,'')}
    return(x)
  }

  coefsx[]<-sapply(coefsx,transformfch)
  coefsx<-round(coefsx,2)

  ss<-rownames(coefsx)
  maxss <- max(sapply(ss, nchar))
  adspaceText <- function(x, n) {paste0(x,"<span style='font-size:12pt;color:black'>",' ',str_c(rep("\\-",(n-nchar(x))),collapse = ''),"</span>","\\:")}

  ss <- sapply(ss, adspaceText, maxss+1)

  maxdigits=ceiling(log10(max(abs(coefsx))))
  Tab <- data.frame(Symbol = ss)
  for (i in c(1:ncol(coefsx))) {
    Tab[colnames(coefsx)[i]] <-
      paste(sapply(coefsx[, i],padding,maxdigits),
            fdrsx[,i], sep = "")
  }
  rownames(Tab)= ss
  colnames(Tab)= c('Symbol', colnames(coefsx))

  return(Tab)
}
