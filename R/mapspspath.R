#'Maps species specific pathways
#'This is important because sometimes some species may not have some pathways
#'@param sps species name
#'@param rxn_pathway_df dataframe having pathway information on metabolites as obtained from pathinfo
#'@importFrom KEGGREST keggList
#'@export
#'@examples
#'org_pathways_map= mapspspath(sps=sps,rxn_pathway_df= res_path)

mapspspath = function(sps, rxn_pathway_df) ###rxn_pathway_df is the result from pathinfo.R
{
  org_pathways = data.frame(keggList("pathway", sps))
  colnames(org_pathways)[1] = "pathway"
  org_pathways$path_id  = rownames(org_pathways)
  rownames(org_pathways) = NULL
  org_pathways$path_id=gsub("path:", "", org_pathways$path_id)
  ###map reaction pathwya to sps pathway
  sps_path_match= rxn_pathway_df[gsub("rn", "", rxn_pathway_df$pathway_id) %in% gsub(sps, "", org_pathways$path_id),]
  sps_path_match$sps_path_id = gsub("rn", sps, sps_path_match$pathway_id)

  return (sps_path_match)
}
