
#'Gives details of reactions for the metabolites
#'@param class_mapped_data metenrichment score dataframe having metabolites, metabolite class and their enrichment score obtained 
#'from metclassenrichment
#'@importFrom KEGGREST keggGet
#'@importFrom tidyr uncount
#'@export
#'@examples
#'res= rxninfo(df_metenrichment)

rxninfo = function(class_mapped_data)
{
  
  aa = class_mapped_data
  aakegg = aa[["KEGG ID"]]
  aakegg = as.vector(aakegg)
  
  ### pass the argument in list of 10s since keggrest takes 10 inputs
  aakegg = split(aakegg,  ceiling(seq_along(aakegg)/10))
  #aakegg=llply(aakegg, as.list)
  
  ### kegg info
  info = llply(aakegg, function(x)keggGet(x))
  unlist_info <- unlist(info, recursive = F)
  extract_info <- lapply(unlist_info, '[', c("ENTRY",'REACTION'))
  dd=do.call(rbind, extract_info)
  df = as.data.frame(dd)
  names(df)[2] = 'REACTION'
  aa[['rxn']] = df[,'REACTION'][match(aa[, 'KEGG ID'], df[, 'ENTRY'])]
  r = sapply(aa[,'rxn'], paste0, collapse=" ")
  for (i in 1: length(r))
  {
    aa$countrxn[i] = length(strsplit(r[[i]], " ")[[1]])
  }
  aa[,'rxn'] = sapply(aa[,'rxn'], paste0, collapse=" ")
  #aa2 = aa[rep(seq_len(nrow(aa)), aa$countrxn), ]
  aa3=aa %>%
    uncount(weights = countrxn, .id = "n", .remove = F) 
  for (j in 1:nrow(aa3)){
    
    aa3[j, 'rxn']= strsplit(aa3[,'rxn'][j], " ")[[1]][aa3$n[j]]
  }
  
  aa3 = aa3[,setdiff(colnames(aa3), c("countrxn","n"))]

  return(aa3)
}
