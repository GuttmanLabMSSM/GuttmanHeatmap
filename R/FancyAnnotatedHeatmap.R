
FancyAnnotatedHeatmap <-function(
    grafname, # grafname
    matx, # matrix of expression
    annot.col, # annotation data
    colsplitx=NULL, # splitting column
    column_ha = NULL,
    BigTab = NULL,
    wdt=NULL, # width
    hgt=NULL, # height
    setname='', # description of the gene set
    column_names_max_height = unit(5, "cm"),
    BigTabcols=NULL,
    cfx = NULL, # lgfch of contrasts of interest
    fdx  = NULL, # fdr or pvalue for the contrasts of interest
    ShowColumnNames = F,
    fnt_size_title = 8,
    scalerows=TRUE,
    rampvalues = 2,
    ...
    #tabNames = NULL # example of 3 contrasts,...

){
  if (!is.null(BigTabcols)){
    allBigTabcols<-c()
    for (i in BigTabcols){
      allBigTabcols<-append(allBigTabcols,colnames(BigTab)[str_detect(colnames(BigTab),i)])
    }
    BigTab<-BigTab[allBigTabcols]
  }

  ncontrasts = dim(cfx)[2]

  # number of genes
  n.genes = dim(matx)[1]

  if (is.null(wdt)){wdt=dim(matx)[2]*2+dim(cfx)[2]/10}
  if (is.null(hgt)){hgt=max(8,n.genes/4)}

  mat = matx[,rownames(annot.col)]

  scale_rows<-function (x)
  {
    m = apply(x, 1, mean, na.rm = T)
    s = apply(x, 1, sd , na.rm = T)
    return((x - m)/s)
  }

  if(scalerows == TRUE){
  mat=scale_rows(mat)}else if(scalerows == FALSE)
  {
    mat = mat
  }

if(!is.null(cfx)){

  coefsx = cfx[rownames(mat),]
  fdrsx = fdx[rownames(mat),]
  require(stringr)
  require(weights)

  # transform lgfch into fch
  transformfch <- function(lgfch) {
    fch <- sign(lgfch) * (2^abs(lgfch))
    return(fch)
  }

  # padding function
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


  p.cuts.stars=c(0.001,0.01, 0.05, 0.1,1)
  p.symbols.stars= c(
    "<span style='font-size:8pt;color:black'>***</span>",
    paste0("<span style='font-size:8pt;color:black'>**</span>","<span style='font-size:8pt;color:white'>*</span>"),
    paste0("<span style='font-size:8pt;color:black'>*</span>","<span style='font-size:8pt;color:white'>**</span>"),
    paste0("<span style='font-size:8pt;color:black'><sup>+</sup></span>","<span style='font-size:8pt;color:white'>**</span>"),
    paste0("<span style='font-size:8pt;color:black'></span>","<span style='font-size:8pt;color:white'>***</span>")
  )

  coefsx[coefsx==0]<-1e-10

  fdrsx[] <- sapply(fdrsx,starmaker, p.levels =p.cuts.stars, symbols = p.symbols.stars)



  coefsx[]<-sapply(coefsx,transformfch)
  coefsx<-round(coefsx,2)

  ss<-rownames(coefsx)
  maxss <- max(sapply(ss, nchar))
  adspaceText <- function(x, n) {
    paste0(
      x,
      "<span style='font-size:8pt;color:black'>",
      ' ',
      str_c(rep("\\-",(n-nchar(x))),collapse = ''),
      "</span>","\\:"
    )
    }

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

}else{
  Tab = cbind.data.frame(Gene = rownames(matx))
  }

  require(circlize)
  require(ComplexHeatmap)
  col_fun = colorRamp2(c(-1*rampvalues,0,rampvalues), c("blue",'white', "red"),space= 'sRGB')
  scalename = 'Z-score'
  if(scalerows==FALSE){

    col_fun = colorRamp2(c(min(mat)+1,mean(mat),max(mat)-1), c("blue",'white', "red"),space= 'sRGB')
    scalename =  'log2-expression'
  }
  row_dend = as.dendrogram(hclust(dist(mat)))

  ht = ComplexHeatmap::Heatmap(mat,
                               name = scalename,
                               heatmap_legend_param = list(
                                 legend_height = unit(2.5, "cm"),
                                 legend_width  = unit(0.5, "cm"),
                                 legend_gp = gpar(fontsize = 8),
                                 title_gp = gpar(fontsize = 8, fontface = "bold")
                               ),
                               # row_km = 3,
                               column_split = colsplitx ,
                               column_title_gp = gpar(fontsize =fnt_size_title),
                               col = col_fun,
                               top_annotation = column_ha ,
                               cluster_rows = T,
                               cluster_columns = F,
                               show_row_names = F,
                               show_column_names = ShowColumnNames,
                               column_order = (rownames(annot.col)),
                               use_raster = TRUE)


  ha_names = rowAnnotation(
    Gene = anno_text(
      gt_render(Tab[,1], align_widths = TRUE),
      location = 0,
      just = 'left',
      gp = gpar(
        fill ='white',
        border= 'white',
        col = "black",
        fontface = 2,
        fontsize = 8,
        fontfamily = "mono"
        )
      )
    )

  htlist=ht+ha_names

  if(!is.null(cfx)){
  for (i in 1:ncontrasts){

    thisha=rowAnnotation(
      x=anno_text(
        gt_render(
          Tab[,i+1],
          align_widths = FALSE,
          gp = gpar(box_col='white',fontface=2,fontsize=8,fontfamily='mono')
          ),
        just = "left",
        show_name = FALSE
        )
      )
    names(thisha)<-paste('Tab',i,sep='')
    htlist=htlist+thisha
  }
}

  #---------------------------------------------
  # export heatmap to a pdf file
  #---------------------------------------------
  pdf(file = grafname,width=wdt,height=hgt)

  #----------------------------------------------
  # draw the heatmap
  #----------------------------------------------
  draw(
    htlist,
    padding = unit(c(8,8,8,8), "mm"), # top, right, bottom, left
    ht_gap =  unit(0.6, "mm"),
    column_title = setname,
    column_title_side = 'bottom'
    )
  #----------------------------------------------
  # decorate the annotations
  #----------------------------------------------
  if(!is.null(cfx)){
  for (i in 1:ncontrasts){
    nm=paste('Tab',i,sep='')
    decorate_annotation(nm, {
      grid.text(
        gsub('\\.','\n',substring(colnames(Tab)[i+1],7)),
        y = unit(1, "npc") + unit(2, "mm"),
        x = unit(0, "npc") + unit(2, "mm"),
        just = "bottom",
        hjust =0,
        rot=0,
        gp = gpar(
          box_col = "white",
          box_lwd =1,
          fontface=2,
          fontsize=8,
          fontfamily='mono'
          )
        )
    })
  }
  }


  dev.off()

}
