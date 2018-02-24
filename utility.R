require(gmodels)
require(Matrix)

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



MultiMerge = function(path) {
  file.names = list.files(path = path, full.names = TRUE)
  data.list = lapply(file.names, function(x) {
    read.table(file = x, header = TRUE)
  })
  
  Reduce(function(x, y) {
    merge(x, y, by=c('Hugo_Symbol', 'Entrez_Gene_Id'))
  }, data.list)
}