
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
    column_names_max_height = unit(8, "cm"),
    BigTabcols=NULL,
    cfx = NULL, # lgfch of contrasts of interest
    fdx  = NULL, # fdr or pvalue for the contrasts of interest
    ShowColumnNames = F,
    fnt_size_title = 20,
    fnt_size_decoration=12,
    scalerows=TRUE,
    rampvalues = 1.5,
    bottom_annotation = NULL,
    border = TRUE,  # Add borders around cells
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.rect(x, y, width, height, gp = gpar(fill = "white", col = NA))
    },
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

 # the number of contrasts in the contrast matrix
  ncontrasts = dim(cfx)[2]

  # number of genes
  n.genes = dim(matx)[1]

  # automatic setting of the column dimendions
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
    Tab = doTab(fdrsx = fdx , coefsx = cfx)
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
                               use_raster = TRUE,
                               bottom_annotation = bottom_annotation,
                               border = border,
                               cell_fun = cell_fun)

  # create a column with gene names
  ha_names = rowAnnotation(
    Gene = anno_text(gt_render(Tab[,1], align_widths = TRUE),
                     location=0,just='left',
                     gp = gpar(
                       fill ='white',
                       border= 'white',
                       col = "black",
                       fontface = 2,
                       fontsize = 12,
                       fontfamily = "mono"
                     )
    )
  )

  # this first iteration just generate the heatmap with the genenames
  htlist=ht+ha_names

  if(!is.null(cfx)){
    for (i in 1:ncontrasts){

       column_box = rowAnnotation(
        x = anno_empty(),
        width = unit(0.01, "cm")
       )



      thisha=rowAnnotation(
        x=anno_text(
          gt_render(Tab[,i+1],
                    align_widths = FALSE,
                    gp = gpar(
                      border='white',
                      box_col='white',
                      box_lwd =1,
                      fontface=2,
                      fontfamily='mono'
                    )
          ),
          just = 'left',
          show_name = FALSE

        )#, x=anno_empty(border=TRUE)
      )


      #
      #                     )),
      # gp = gpar(box_col = "white", box_lwd =1, fontface=2, fontsize = 12,fontfamily = "mono"),
      #                     just = "left",show_name = FALSE))




      names(thisha)<-paste('Tab',i,sep='')

        htlist=htlist+column_box+thisha

    }
  }#htlist = htlist+column_box
  htlist = htlist+column_box




  # write the file to a pdf file
  pdf(file = grafname,width=wdt,height=hgt)

  draw(htlist,
       #padding = unit(c(8, 8, 12, 12), "mm"),
       padding = unit(c(8, 8, 12, 12), "mm"),

       column_title = setname,
       column_title_side = 'bottom',
       heatmap_legend_side = "bottom",
       gap=unit(0.25,'mm')
  )




  if(!is.null(cfx)){
    for (i in 1:ncontrasts){

      nm=paste('Tab',i,sep='')



      # draw a table around the column
      decorate_annotation(nm, {

        AnnotationText =  gsub('\\.','\n',substring(colnames(Tab)[i+1],7))

        grid.text(
          AnnotationText,
          x = unit(0, "npc") + unit(2, "mm"),
          y = unit(1, "npc") + unit(2, "mm"),
          just = "bottom",
          hjust = 0,
          rot=0,
          gp = gpar(
            box_col = "black",
            box_lwd =1,
            fontface=2,
            fontsize=fnt_size_decoration,
            fontfamily='mono'
          )
        )

        grid.lines(c(0,0),c(1,1.08),
                   gp=gpar(col = "black",lwd=1.1)
         )

        grid.lines(c(1,1),c(1,1.08),
                   gp=gpar(col = "black",lwd=1.1)
        )

        grid.lines(c(0,1),c(1.08,1.08),
                   gp=gpar(col = "black",lwd=1.1)
        )

         grid.lines(c(0,1),c(1.001,1.001),
                   gp=gpar(col = "black",lwd=1)
         )

         grid.lines(c(0,1),c(-0.001,-0.001),
                    gp=gpar(col = "black",lwd=1)
         )
        #   gp=gpar(col = "black")
        # )



      })

    }
  }
  dev.off()

}
