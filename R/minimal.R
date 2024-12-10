minimal<-function(){
require(tidyverse)
  mat = cbind.data.frame(
    G1 = rnorm(3),
    G2 = rnorm(3),
    G3 = rnorm(3),
    G4 = rnorm(3),
    G5 = rnorm(3)
  )

  mat= as.matrix(t(mat))

  colnames(mat) = c('A','B','C')

  annot = cbind.data.frame(
    Group = c('A','B','C')
  )

  row.names(annot) = c('A','B','C')


  lgfch.matrix = cbind.data.frame(
    AvsB = runif(5,-2,2),
    AvsC = runif(5,-2,2),
    BvsC = runif(5,-2,2)
  ) %>% as.matrix()

  fdr.matrix = cbind.data.frame(
    AvsB = runif(5,0,1),
    AvsC = runif(5,0,1),
    BvsC = runif(5,0,1)
  ) %>% as.matrix()

  rownames(lgfch.matrix) = c('G1','G2','G3','G4','G5')
  rownames(fdr.matrix) = rownames(lgfch.matrix)

return(list(mat = mat, lgfch = lgfch.matrix, fdr = fdr.matrix, annot = annot))

}


