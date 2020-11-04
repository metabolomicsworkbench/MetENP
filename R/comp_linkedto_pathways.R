#'counts total no of compounds that are linked to KEGG pathways
#'@importFrom KEGGREST keggList
#'@importFrom KEGGREST keggGet
#'@export
#'@examples
#'ls_path =comp_linkedto_pathways()

comp_linkedto_pathways <- function()
{
kegg_compounds = keggList("compound")
c_names_df = data.frame()
for (i in 1:length(kegg_compounds))
{
  
  c_names_df[i,"names"] = names(kegg_compounds[i])
}
query_c = gsub("cpd:","",c_names_df$names)
query_c = as.vector(query_c)
query_c = unique(query_c)
query_split = split(query_c,  ceiling(seq_along(query_c)/10))
info = llply(query_split, function(x)keggGet(x))  

ls_path = list()
for (i in 1:length(info))
{
  for (j in 1:length(info[[i]])) ### because we have splited the queries in 10s
  {
    if ("PATHWAY" %in% names(info[[i]][[j]]))
    {
      ls_path = c(ls_path,c(info[[i]][[j]][["ENTRY"]]))
      
    }
  }}
return(ls_path)
}

      

