require(Matrix)
require(igraph)
require(gmodels)
require(ggplot2)
require(sva)
require(RANN)
require(reshape)
require(preprocessCore)


rnaObj <- setClass("rnaObj", slots =
                     c(data = "data.frame", tissue = "vector", study = "vector", counts = "data.frame",
                       pca.scores = "data.frame", pca.load = "data.frame", percentage = "vector"))


setGeneric("initialize", function(object, min.expression = 0.01, keep.tissue = c(), log = FALSE, normalize = FALSE, center = TRUE, scale = TRUE) 
  standardGeneric("initialize"))
setMethod("initialize", "rnaObj",
          function(object, min.expression = 0.01, keep.tissue = c(), log = FALSE, normalize = FALSE, center = TRUE, scale = TRUE) {
            print("Initializing S4 object")
            tmp = object@data
            # print(dim(tmp))
            
            # remove unwanted tissues
            if (length(keep.tissue) != 0) {
              print("removing unwanted tissues")
              cols.use = unlist(lapply(keep.tissue, function(x) colnames(data)[grep(x, colnames(data))]))
              tmp = tmp[, cols.use]
              rm(cols.use)
            }
            
            # gene filtering
            # remove genes with row means less than min expression level
            print("removing genes below mean expression level")
            # tmp = tmp[-which(rowMeans(tmp[,]) <= min.expression), ]
            print(dim(tmp))
            
            # if log then take log of data
            if (log) {
              print("applying log2 transformation")
              tmp = log2(tmp + 1)
            }
            
            # if normalize is true then do quantile normalization
            if (normalize) {
              print("applying quantile normalization per gene")
              # tmp = t(scale(t(tmp),  center = center, scale = scale))
              tmp.rownames = rownames(tmp)
              tmp.colnames = colnames(tmp)
              tmp = normalize.quantiles(as.matrix(tmp))
              tmp = data.frame(tmp, row.names = tmp.rownames)
              colnames(tmp) = tmp.colnames
              
              rm(tmp.rownames, tmp.colnames)
            }
            
            object@data = data.frame(tmp)
            rm(tmp)
            
            print("calculating tissue counts")
            object@study = factor(unlist(lapply(colnames(object@data), function(x) unlist(strsplit(x, "\\.")[1][1])[1])))
            object@tissue = factor(unlist(lapply(colnames(object@data), function(x) unlist(strsplit(x, "\\.")[1][1])[2])))
            
            study.table = table(object@study)
            tissue.table = table(object@tissue)
            
            total.gtex = study.table[1]
            total.tcga = study.table[2]
            
            gtex.tissues = colnames(object@data)[grep("GTEX", colnames(object@data))]
            gtex.tissues = lapply(gtex.tissues, function(x) unlist(strsplit(x, "\\.")[1][1]))
            
            tcga.tissues = colnames(object@data)[grep("TCGA", colnames(object@data))]
            tcga.tissues = lapply(tcga.tissues, function(x) unlist(strsplit(x, "\\.")[1][1]))
            
            count.gtex = as.vector(table(unlist(lapply(gtex.tissues, function (x) x[2]))))
            count.tcga = as.vector(table(unlist(lapply(tcga.tissues, function (x) x[2]))))
            
            unique.tissues = levels(object@tissue)
            
            count = data.frame(GTEX = count.gtex, TCGA = count.tcga, row.names = unique.tissues)
            total.row = data.frame(GTEX = total.gtex, TCGA = total.tcga, row.names = c("total"))
            count = rbind(count, total.row)
            
            object@counts = count
            rm(count, total.gtex, total.tcga, count.gtex, count.tcga, unique.tissues, total.row, study.table, tissue.table, gtex.tissues, tcga.tissues)
            
            # print(dim(object@data))
            
            return(object)
          })


setGeneric("doPCA", function(object, pcs.store = 100) standardGeneric("doPCA"))
setMethod("doPCA", "rnaObj",
          function(object, pcs.store=100) {
            data.use = object@data
            pca.obj = fast.prcomp(t(data.use), center = FALSE, scale = FALSE)
            perc = 100*(pc.obj$sdev)^2 / sum(pc.obj$sdev^2)

            object@percentage = (100*(pc.obj$sdev)^2 / sum(pc.obj$sdev^2))[1:pcs.store]
            object@pca.scores = data.frame(pca.obj$x[ ,1:pcs.store])
            object@pca.load = data.frame(pca.obj$rotation[ ,1:pcs.store])
            
            rm(perc, pca.obj)

            return(object)
          })


# Loads and merges RNA data from a list of file names and paths
FileMultiMerge = function(file.names, path) {
  tissueTypes = c("bladder", "breast", "prostate", "thyroid")
  path.names = lapply(file.names, function(file.name) {
    paste(path, file.name, sep = "/")
  })
  
  data.list = lapply(path.names, function(file.name) {
    read.table(file = file.name, header = TRUE)
  })
  
  data = Reduce(function(x, y) {
    merge(x, y, by=c('Hugo_Symbol', 'Entrez_Gene_Id'))
  }, data.list)
  
  return(data)
}


# Loads and merges RNA data from all files in given path
PathMultiMerge = function(path) {
  file.names = list.files(path = path, full.names = TRUE)
  data.list = lapply(file.names, function(x) {
    read.table(file = x, header = TRUE)
  })
  
  data = Reduce(function(x, y) {
    merge(x, y, by=c('Hugo_Symbol', 'Entrez_Gene_Id'))
  }, data.list)
  
  return(data)
}