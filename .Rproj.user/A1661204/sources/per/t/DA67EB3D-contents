#  This function will normalize Ct values from TLDA or RT-PCR
#' the function outputs a matrix of normalized Ct values according to a house keeping gene
#' @description For a matrix of TLDA or RT-PCR expression it will normalize the Ct values for running linear models
#' @param mat is a matrix with sample ids as rownames and genes as column names
#' @param hk is the name of the house keeping gene
#' @examples


deltaCtNormalize<-function(mat,hk='RPLP0'){


  for(i in 1:dim(mat)[1]){
    for(j in 1:dim(mat)[2]){
      mat[i,j] = ifelse(mat[i,j]==40,NA,mat[i,j])
    }
  }

  mat = mat[,hk]-mat


  # replace NA by 20 % of the minimum unlogged values

  for(i in 1:dim(mat)[2]){

    if ( sum( is.na(mat[,i]) )>0 ){

      mat[which( is.na(mat[,i])==T )  , i] <- log2(0.2*(2^min(mat[,i],na.rm=T)))

    }

  }

  return(mat)

}
