
doAnnotatedHeatmap <-function(
    grafname, # grafname
    matx, # matrix of expression
    annot.col, # annotation data
    colsplitx=NULL, # splitting column
    column_ha = NULL,
    BigTab = NULL,
    wdt=0, # width
    hgt=0, # height
    setname='', # description of the gene set
    column_names_max_height = unit(8, "cm"),
    BigTabcols=NULL,
    cfx = NULL, # lgfch of contrasts of interest
    fdx  = NULL, # fdr or pvalue for the contrasts of interest
    ShowColumnNames = F,
    fnt_size_title = 20,
    scalerows=TRUE,
    rampvalues = 1.5,
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

  #cfxfunc<-function(x){substr(x,1,5)=='lgFCH'}
  #fdxfunc<-function(x){substr(x,1,5)=='pvals'}
  #cfx=BigTab[,colnames(BigTab)[unlist(lapply(colnames(BigTab),cfxfunc))]]
  #fdx=BigTab[,colnames(BigTab)[unlist(lapply(colnames(BigTab),fdxfunc))]]

  ncontrasts = dim(cfx)[2]

  # number of genes
  n.genes = dim(matx)[1]

  if (wdt==0){wdt=dim(matx)[2]*2+dim(BigTab)[2]/10}
  if (hgt==0){hgt=max(8,n.genes/4)}

  mat = mat[,rownames(annot.col)]

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
  Tab = doTab(fdrsx = fdx , coefsx = cfx)
}else{Tab = cbind.data.frame(Gene = rownames(matx))}

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


  ha_names = rowAnnotation(Gene = anno_text(gt_render(Tab[,1], align_widths = TRUE),location=0,just='left' ,
                                            gp = gpar(fill ='white', border= 'white', col = "black", fontface = 2,  fontsize = 12,
                                                      fontfamily = "mono")))
  htlist=ht+ha_names

  if(!is.null(cfx)){
  for (i in 1:ncontrasts){

    thisha=rowAnnotation( x=anno_text( gt_render((Tab[,i+1]), align_widths = FALSE),
                                       gp = gpar(box_col = "white", box_lwd =1, fontface=2, fontsize = 12,
                                                 fontfamily = "mono"),
                                       just = "left",show_name = FALSE))
    names(thisha)<-paste('Tab',i,sep='')
    htlist=htlist+thisha
  }
}

  # write the file to a pdf file
  pdf(file = grafname,width=wdt,height=hgt)
  draw(htlist,padding = unit(c(8, 8, 12, 12), "mm") ,  column_title = setname ,column_title_side = 'bottom')

  if(!is.null(cfx)){
  for (i in 1:ncontrasts){
    nm=paste('Tab',i,sep='')
    decorate_annotation(nm, {
      grid.text(gsub('\\.','\n',substring(colnames(Tab)[i+1],7)), y = unit(1, "npc") + unit(2, "mm"), just = "bottom",hjust =0,rot=0,
                gp = gpar(box_col = "white", box_lwd =1, fontface=2,fontsize=10,fontfamily='mono'))
    })
  }
  }


  dev.off()

}
