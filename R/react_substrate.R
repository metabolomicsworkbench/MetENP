#'Gets information whether metabolite is a product or substrate
#'@param met_gene_info dataframe with gene information obtained from enzyme_gene_info
#'@export
#'@examples
#'rclass_info = react_substrate(met_gene_info)

react_substrate = function(met_gene_info){
debug = 0; # Mano
rclass_info= met_gene_info
rclass_info$reactant_product=NA
for (i in 1:nrow(met_gene_info)){
  
  id=grep(met_gene_info[["KEGG ID"]][i],strsplit(as.character(met_gene_info[["EQUATION"]][i])," <=> ")[[1]])
  #id2=agrep(rclass_info[["KEGG ID"]][i], strsplit(rclass_info[["RCLASS1"]], "_")[[i]])
  # Mano: 2023/08/15: For one case, length(id) was > 1, so, error; if so ask it to print the KEGG ID and id and i
  if(length(id) == 1){ # Mano: added condition
    if (id==1){
      rclass_info[["reactant_product"]][i]="Substrate"
    }else if (id==2 ){
      rclass_info[["reactant_product"]][i]="Product" 
    }
  } else { # Mano
    if(length(id) > 1){ # Mano: added condition
      rclass_info[["reactant_product"]][i]="Substrate,Product"; # Mano: in both
      if(debug > 0) {
        print("length(id) > 1 detected: printing the row of met_gene_info:");
        print(paste0("i = ", i));
        print(paste0("KEGG ID:", met_gene_info[["KEGG ID"]][i]));
        print(paste0("EQUATION:", met_gene_info[["EQUATION"]][i]));
        print(as.character(met_gene_info[i,]));
      }
    }
  }
}
return(rclass_info)
}
