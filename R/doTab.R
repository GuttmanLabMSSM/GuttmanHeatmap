doTab<-function(fdrsx,coefsx){
  # loading libraries that are needed
  require(stringr)
  require(weights)

  #----------------------------------------------------------
  # functions to be used inside this function
  #----------------------------------------------------------

  # a function to convert log-fold changes into fold-changes
  transformfch <- function(lgfch) {
    fch <- sign(lgfch) * (2^abs(lgfch))
    return(fch)
  }

  # a padding function just to fill some gaps
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

  # a function to pad stars
  padstars<-function(x){
    while (nchar(x)<3){x<-paste(x,'')}
    return(x)
  }

  # a function to pad spaces
  adspaceText <- function(x, n) {
    paste0(x,"<span style='font-size:12pt;color:black'>",' ',str_c(rep("\\-",(n-nchar(x))),collapse = ''),"</span>","\\:")
  }

  # creating the map between the p-values thresholds and stars
  p.cuts.stars=c(0.001,0.01, 0.05, 0.1,1)
  p.symbols.stars= c('***',
                     paste0('**',"<span style='font-size:12pt;color:white'>",'*',"</span>",collapse = ''),
                     paste0('*',"<span style='font-size:12pt;color:white'>",'*',"</span>","<span style='font-size:12pt;color:white'>",'*',"</span>",collapse = ''),
                     paste0('+',"<span style='font-size:12pt;color:white'>",'*',"</span>","<span style='font-size:12pt;color:white'>",'*',"</span>",collapse = ''),
                     paste0("<span style='font-size:12pt;color:white'>",'*',"</span>","<span style='font-size:12pt;color:white'>",'*',"</span>","<span style='font-size:12pt;color:white'>",'*',"</span>",collapse = '')
  )
  #"*  ", "+  ","   ")

  # to avoid problems with 0 values
  coefsx[coefsx==0]<-1e-10

  # replace p-values with stars
  fdrsx[] <- sapply(fdrsx,starmaker, p.levels =p.cuts.stars, symbols = p.symbols.stars)


  coefsx[]<-sapply(coefsx,transformfch)
  coefsx<-round(coefsx,2)

  ss<-rownames(coefsx)
  maxss <- max(sapply(ss, nchar))


  ss <- sapply(ss, adspaceText, maxss+1)

  maxdigits=ceiling(log10(max(abs(coefsx))))

  # a for loop to create a table with all coefficients and p-values
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
