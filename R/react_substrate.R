#'Gets information whether metabolite is a product or substrate
#'@param met_gene_info dataframe with gene information obtained from enzyme_gene_info
#'@export
#'@examples
#'rclass_info = react_substrate(met_gene_info)

react_substrate = function(met_gene_info){
rclass_info= met_gene_info
rclass_info$reactant_product=NA
for (i in 1:nrow(met_gene_info)){
  
  id=grep(met_gene_info[["KEGG ID"]][i],strsplit(as.character(met_gene_info[["EQUATION"]][i])," <=> ")[[1]])
  #id2=agrep(rclass_info[["KEGG ID"]][i], strsplit(rclass_info[["RCLASS1"]], "_")[[i]])
  if (id==1){
    rclass_info[["reactant_product"]][i]="Substrate"
  }else if (id==2 ){
    rclass_info[["reactant_product"]][i]="Product"  
  }
  
}
return(rclass_info)
}
