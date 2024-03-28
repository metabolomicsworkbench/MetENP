#' Gives information on pathways linked to significant metabolites specific to species of choice
#'@param df_metenrichment metenrichment score dataframe having metabolites, metabolite class and their enrichment score obtained 
#'from metclassenrichment
#'@param sps species name 
#'@param debug 1 to print some informative lines, 0 to suppress them
#'@export
#'@examples
#'met_path = met_pathways(df_metenrichment = metenrichment, sps = 'hsa', debug = 0)

met_pathways<- function(df_metenrichment,sps, debug = 0)
{
if(debug > 0) { tic("rxninfo:"); }
res= rxninfo(df_metenrichment)
if(debug > 0) { toc(); tic("pathinfo:"); }
if(debug > 0) { print("Res:"); print(res$rxn); }
res_path= pathinfo(res, 'PATHWAY') ###rxn pathway
names(res_path)[c(1,2)] = c("Rxn_id","Rxn_name")
if(debug > 0) { toc(); tic("mapspspath:"); }
org_pathways_map= mapspspath(sps=sps,rxn_pathway_df= res_path)
if(debug > 0) { toc(); tic("merge:"); }
met_path = merge(res, org_pathways_map, by.x = 'rxn', by.y="Rxn_id")
if(debug > 0) { toc(); }
return(unique(met_path))
}

#
# Additional notes:
# If the resulting list of reactions for the metabolites is large, keggGet may throw an error like:
# Error in .getUrl(url, .flatFileParser): Forbidden (HTTP 403).
# Then, may be pass a subset of the dataframe as below; 
#met_path = met_pathways(df_metenrichment = metenrichment[1:80, ], 'hsa', debug = 1)

