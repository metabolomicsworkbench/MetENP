#' Gives information on pathways linked to significant metabolites specific to species of choice
#'@param df_metenrichment metenrichment score dataframe having metabolites, metabolite class and their enrichment score obtained 
#'from metclassenrichment
#'@param sps species name 
#'@export
#'@examples
#'met_path = met_pathways(df_metenrichment = metenrichment, 'hsa')


met_pathways<- function(df_metenrichment,sps)
{
res= rxninfo(df_metenrichment)
res_path= pathinfo(res, 'PATHWAY') ###rxn pathway
names(res_path)[c(1,2)] = c("Rxn_id","Rxn_name")
org_pathways_map= mapspspath(sps=sps,rxn_pathway_df= res_path)
met_path = merge(res, org_pathways_map, by.x = 'rxn', by.y="Rxn_id")
return(unique(met_path))
}
