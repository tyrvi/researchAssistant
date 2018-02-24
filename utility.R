require(gmodels)
require(Matrix)

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

# rnaObj <- setClass("rnaObj", slots =
#                      c(data = "data.frame", pca.scores = "data.frame", pca.load = "data.frame",
#                        meta = "data.frame", tissue = "vector", genes = "vector", group = "vector",
#                        percentage = "Matrix"))
# 
# 
# setMethod("doPCA", "rnaObj",
#           function(object, pcs.store=100) {
#             data.use = object@data
#             pca.obj = fast.prcomp(t(data.use), center = FALSE, scale = FALSE)
#             perc = 100*(pc.obj$sdev)^2 / sum(pc.obj$sdev^2)
# 
#             object@percentage = Matrix(100*(pc.obj$sdev)^2 / sum(pc.obj$sdev^2))[1:pcs.store]
#             object@pca.scores = data.frame(pca.obj$x[ ,1:pcs.store])
#             object@pca.load = data.frame(pca.obj$rotation[ ,1:pcs.store])
# 
#             return(object)
#           })