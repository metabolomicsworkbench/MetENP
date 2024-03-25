#'Provides pathway or orthology information from reaction data of significant metabolite
#'@param res dataframe having metabolite information along with reaction information as obtained from rxninfo 
#'@param col_name "PATHWAY" or "ORTHOLOGY"
#'@importFrom KEGGREST keggGet
#'@importFrom tidyr uncount
#'@export
#'@examples
#'res_path= pathinfo(res, 'PATHWAY') ###rxn pathway

pathinfo = function(res, col_name)
{
  
  query = as.vector(res$rxn)
  query = unique(query)
  query = query[!is.na(query)]; # remove NAs
  query = unique(query[query != ""]) # exclude empty
  ### pass the argument in list of 10s since keggrest takes 10 inputs
  query_split = split(query,  ceiling(seq_along(query)/10))
 
  
  ### kegg info
  info = llply(query_split, function(x)keggGet(x))
  unlist_info <- unlist(info, recursive = F)
  extract_info <- lapply(unlist_info, '[', c("ENTRY","NAME","PATHWAY","RCLASS","ORTHOLOGY","DEFINITION","EQUATION","ENZYME"))
  dd=do.call(rbind, extract_info)
  df = as.data.frame(dd)
  ## BUG FIX; Added by Sumana Srinivasan 0n 04042022
  ### named all columns to be consistent
  names(df)[1] = 'ENTRY'
  names(df)[2] = 'NAME'
  names(df)[3] = 'PATHWAY'
  names(df)[4] = 'RCLASS'            
  names(df)[5] = 'ORTHOLOGY'
  names(df)[6] = 'EQUATION_more'
  names(df)[7] = 'EQUATION'
  names(df)[8] = 'ENZYME' 
  ## BUG FIX; Added by Sumana Srinivasan 0n 04042022, for rno, the ORTHOLOGY ENZYME column, keggGET returns NA and hence have to be named
  
  ###pathway
  if (col_name == "PATHWAY"){
  path_info = df[,c("ENTRY","NAME","PATHWAY")]
  }else{
    path_info = df[,c("ENTRY","NAME","RCLASS","ORTHOLOGY","EQUATION",'EQUATION_more',"ENZYME")]
  }
  
  r = sapply(path_info[,col_name], paste0, collapse=";")
  path_info[["id"]]=NA
  for (i in 1: length(r))
  {
    path_info$countrxn[i] = length(strsplit(r[[i]], ";")[[1]])
    if (!is.null(names(extract_info[[i]][[col_name]]))){
    path_info[["id"]][i] = do.call(paste, c(as.list(names(extract_info[[i]][[col_name]])), sep = ";"))
    }
  }
  path_info[,col_name] = sapply(path_info[,col_name], paste0, collapse=";")
  #aa2 = aa[rep(seq_len(nrow(aa)), aa$countrxn), ]
  path_info3=path_info %>%
    uncount(weights = countrxn, .id = "n", .remove = F) 
  for (j in 1:nrow(path_info3)){
    
    path_info3[j, col_name]= strsplit(path_info3[,col_name][j], ";")[[1]][path_info3$n[j]]
    path_info3[j, "id"]= strsplit(path_info3[,"id"][j], ";")[[1]][path_info3$n[j]]
  }
  
  path_info3 = path_info3[,setdiff(colnames(path_info3), c("countrxn","n"))]
  
  if (col_name == "PATHWAY"){
    colnames(path_info3)[4] = "pathway_id"
  }
  # else if (col_name == "ORTHOLOGY")
  # {colnames(path_info3)[4] = "orthology_id"}
 
  return(path_info3)
}

##### Test code, Sumana Srinivasan, 05052022 #############
# library(plyr); library(KEGGREST);
# rxn_ids <- c("R00021", "R00086") ## give rxn ids that give NA
# info = llply(rxn_ids, function(x)keggGet(x))
# unlist_info <- unlist(info, recursive = F)
# extract_info <- lapply(unlist_info, '[', c("ENTRY","NAME","PATHWAY","RCLASS","ORTHOLOGY","DEFINITION","EQUATION","ENZYME"))
# dd=do.call(rbind, extract_info)
# df = as.data.frame(dd)
