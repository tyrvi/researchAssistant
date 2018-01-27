MultiMerge = function(path) {
  file.names = list.files(path=path, full.names = TRUE)
  data.list = lapply(file.names, function(x) {
    read.table(file = x, header = TRUE)
  })
  
  Reduce(function(x, y) {
    merge(x, y, by=c('Hugo_Symbol', 'Entrez_Gene_Id'))
  }, data.list)
}