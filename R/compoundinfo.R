#'Counts total no of compound in a pathway that has been linked to significant metabolites
#' This will help in calculating HG score for metabolite enrichment
#'@param met_path dataframe of pathways mapped to significant metabolites as obtained from met_pathways
#'@param sps species name
#'@importFrom utils data
#'@export
#'@examples
#'compound_count = compoundinfo(met_path, sps)

compoundinfo = function(met_path, sps)
{
  #data("korg", package="pathview")
query = as.vector(met_path$sps_path_id)
query = unique(query)
query_split = split(query,  ceiling(seq_along(query)/10))
info = llply(query_split, function(x)keggGet(x))
unlist_info <- unlist(info, recursive = F)
extract_info <- lapply(unlist_info, '[', c("ENTRY","NAME","CLASS","COMPOUND"))
dd=do.call(rbind, extract_info)
df = as.data.frame(dd)
names(df)[3] = 'CLASS'
names(df)[4] = 'COMPOUND'
compound_info = df
r = sapply(compound_info[,'COMPOUND'], paste0, collapse=";")
for (i in 1: length(r))
{
  compound_info$Total_no._of_comps_in_pathway[i] = length(strsplit(r[[i]], ";")[[1]])
  if (!is.null(names(extract_info[[i]][['COMPOUND']]))){
    compound_info$id[i] = do.call(paste, c(as.list(names(extract_info[[i]][['COMPOUND']])), sep = ";"))
  }else
  {compound_info$id[i] = ""}
}
count_met_path = compound_info[compound_info[,'Total_no._of_comps_in_pathway'] > 0,]
# if(sps %in% c("hsa")){
#   spsname=c("Homo sapiens");
# }
nam = korg[,'scientific.name'][which(korg[,'kegg.code']==sps)]

count_met_path[["NAME"]]=gsub(paste0(" - ",nam,".+"), "", count_met_path[["NAME"]])
count_met_path[["NAME"]]=gsub(paste0(" - ",nam), "", count_met_path[["NAME"]])
count_met_path[["NAME"]]=gsub(paste0(" - ",".+"), "", count_met_path[["NAME"]])
return(count_met_path)
}

