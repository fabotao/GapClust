#################
## MIT license ##
#################
# Copyright <YEAR> <COPYRIGHT HOLDER>
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
# associated documentation files (the "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the
# following conditions:
# The above copyright notice and this permission notice shall be included in all copies or substantial
# portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT
# LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
######################################################################################################

##################
## Main program ##
##################
#' Detecting rare cells from expression profiles by single cell RNA-seq.
#'
#' GapClust takes advantage of the gap between minor cluster and neighbouring
#' abundant cluster to let rare cells within minor cluster stand out through
#' delicately designed statistics. Meanwhile, GapClust does not
#' struggle to search for rare cell informative genes like most of the competitors,
#' but learns the cluster size as well as rare cells using simple arithmetic calculation.
#'
#' @param data A numeric matrix-like object of counts which has been normalized, where rows are genes and columns are cells.
#' @param k The upper limit of minor cluster size. k can be adjusted accordingly
#'
#' @return An object with rare cell indices
#'
#' @aliases GapClust
#'

GapClust <- function(data, k=200){

  ## Fano genes for clustering
  pbmc <- CreateSeuratObject(counts = data)
  pbmc <- FindVariableFeatures(object = pbmc, selection.method='vst', nfeatures=dim(data)[1], verbose = F)
  vst <- (pbmc@assays$RNA@meta.features$vst.variance.standardized)
  den <- density(vst)
  features.vst <- dimnames(data)[[1]][vst > find_elbow(den$x[which.max(den$y):length(den$x)], den$y[which.max(den$y):length(den$y)])]
  tmp <- data[dimnames(data)[[1]] %in% (features.vst),]
  tmp <- log2(as.matrix(tmp)+1)
  N = dim(tmp)[2]
  pca <- irlba(t(tmp), nv=min(c(50, dim(tmp)-1))) # More robust no error, contrast to calcul.pca
  rm(list=c('pbmc', 'data', 'tmp'))
  
  pca$pca <-t(pca$d*t(pca$u))
  knn.res <- Neighbour(pca$pca, pca$pca, k=k)

  distance.diff <- (knn.res$distances[, -1, drop = FALSE] - knn.res$distances[, -ncol(knn.res$distances), drop = FALSE])
  diff.left <- distance.diff[, -1, drop = FALSE] - distance.diff[, -ncol(distance.diff), drop = FALSE]
  diff.both <- diff.left[, -ncol(diff.left), drop=FALSE] - diff.left[, -1, drop=FALSE]
  diff.both[,1] <- diff.both[,1] + distance.diff[,1]  # Very important due to distance variation to the first neighbor.

  v1.k <- matrix(NA, N, k-3)
  rm(list=c('pca', 'distance.diff', 'diff.left'))
  gc()
  skew <- c()
  top.values.ave <- c()
  for(j in 1:dim(diff.both)[2]){
    v <- diff.both[,j]
    v1 <- v
    for(m in 1:length(v)){
      v1[m] <- (v[m] + v[knn.res$indices[m,2]])/2
    }
    v1.k[, j] <- (v1)
    v2 <- v1[order(v1, decreasing = T)[(j+2):length(v1)]]
    v2[is.na(v2)] <- 0
    top.values <- v1[knn.res$indices[which.max(v1),1:(j+1)]]
    v2 <- c(v2[v2 <= (quantile(v2, 0.75)+1.5*IQR(v2)) & v2 >= (quantile(v2, 0.25)-1.5*IQR(v2))], rep(sum(top.values[top.values>0])/length(top.values), (2)))
    skew <- c(skew, skewness(v2))
    top.values.ave <- c(top.values.ave, mean(top.values))
  }

  ids <- which(skew > 2)
  if(length(ids) < 1){
    print('No rare cell cluster identified')
    return(NA)
  }
  col.mat <- matrix(0, length(ids), N)
  for(i in 1:length(ids)){
    top.cell <- which.max(v1.k[,(ids[i])])
    col.mat[i, knn.res$indices[top.cell,1:(ids[i]+1)]] <- skew[ids[i]] * top.values.ave[ids[i]]
  }

  id.max <- apply(col.mat, 2, which.max)
  max.val <- apply(col.mat, 2, max)
  id.max[max.val==0] <- 0
  cnt <- table(id.max)
  cnt <- cnt[names(cnt)!='0']
  id.max.match <- cnt[which(cnt == (ids[as.integer(names(cnt))] + 1))] - 1

  cls <- rep(0, N)
  for(id.match in id.max.match){
    cls[id.max==(id.match)] <- which(id.max.match %in% id.match)
  }

  rare.cells <- list()
  for(id.match in id.max.match){
    rare.cells[[as.character(id.match)]] <- knn.res$indices[which.max(v1.k[,id.match]), 1:(id.match+1)]
  }
  results <- list(skewness=skew, rare_cell_indices=rare.cells, rare_score=v1.k)
  return(results)
}


