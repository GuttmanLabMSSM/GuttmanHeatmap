  #  This function will generate zscores for any gene set and other contexts rather than gene expression
  #' The function outputs a data frame in long format with rownames as sample ids, geneset names and zscores
  #' @description A function that transform values into Z scores
  #' @param mat is a matrix of expressions with genes as rownames and samples ids as columns
  #' @param geneset is a vector with gene symbols
  #' @param geneset.name is a string with the name of the gene set
  #' @examples

generate.Zscores = function(mat,geneset,geneset.name){

  # subset the immune genes set
  sub.mat = mat[which(rownames(mat) %in% geneset),]

  # evaluate the zscores for the immune genes set
  z.sub.mat = (sub.mat - apply(sub.mat, 1, mean)) / apply(sub.mat, 1, sd)

  # evaluate the zscores for the entire gene set
  k = dim(z.sub.mat)[1]
  z.scores.mat = apply(z.sub.mat, 2, function(x) sum(x) / sqrt(k))

  # create a data frame in long format with the zscores evaluated for each sample
  z.scores.data = data.frame(geneset.name = geneset.name, z.scores = z.scores.mat)

  return(z.scores.data)
}



