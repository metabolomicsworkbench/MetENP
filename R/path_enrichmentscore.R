#'Calculates pathway enrichment score of the pathways linked to significant metabolites
#'@param met_path dataframe of pathways mapped to significant metabolites as obtained from met_pathways
#'@param ls_path list of all the KEGG compounds linked to pathways. This can be obtained by running comp_linkedto_pathways().R or is also available
#'please use ls_path when you want to calculate HG score depending on all the compounds in KEGG linked to a pathway
#'sysdata.rda
#'@param sps species name
#'@param padj Method for p value adjustment, i.e. "fdr"
#'@param sig_metabolite_kegg_id dataframe of significantly altered metabolites for HG calculation
#'@param refmet_class dataframe of metabolite data with refmet class for a study, alternatively,
#'data obtained from rest can also be given here which is yet not associated to refmet class
#'We need this to get the number of all the metabolites detected in the study for calculating HG score.
#'@param kegg_comp_path TRUE or FALSE depending you want to use the ls_path as your background. Default is FALSE, where background set
#'is the total metabolite in a study
#'@importFrom stats p.adjust
#'@importFrom stats phyper
#'@export
#'@examples
#'kegg_es = path_enrichmentscore(met_path,sig_metabolite_kegg_id, ls_path,refmet_class,sps='hsa',padj='fdr')

path_enrichmentscore <- function(met_path,sig_metabolite_kegg_id,ls_path,refmet_class,sps, padj,kegg_comp_path=FALSE)
{
  compound_count = compoundinfo(met_path, sps)
  met_path_selected = met_path[,c('Metabolite', 'PATHWAY')] %>% distinct
  #met_compound_count=merge(met_path_selected, compound_count, by.x = "PATHWAY", by.y = "NAME")
  comp_path_count = unique(compound_count[,c("NAME","Total_no._of_comps_in_pathway")])
  met_compound_count=partial_join(met_path_selected, comp_path_count, by_x = "PATHWAY", pattern_y = "NAME")
  comp_path_count=comp_path_count[order(comp_path_count$NAME),]
  #levels(met_compound_count[['PATHWAY']])=(factor(met_compound_count[['PATHWAY']]))
  freqtable = cbind(data.frame(table(met_compound_count[['PATHWAY']])), comp_path_count)
  colnames(freqtable)[2] = "No.of mets in study"
  colnames(freqtable)[1] = "Pathway name"
  freqtable[3] = NULL

  if (kegg_comp_path==TRUE){
  N=length(ls_path)### met with pathway annotations, noofcompounds_linked_pathway.R
  }else{
    N=nrow(refmet_class)}
  ### met with pathway annotations, noofcompounds_linked_pathway.R
  for(i in 1:nrow(freqtable)){
    k = freqtable[i,3]## count of met inthe pathway
    M = freqtable[i,2]### altered met in the pathway
    L = nrow(sig_metabolites_kegg_id) ## total altered met
    pp<-phyper(M-1,L, N-L,k, lower.tail=FALSE)
    freqtable[i,"pathway_HG p-value"] <- pp
  }

  freqtable$Padjust = p.adjust(freqtable$pathway_HG, padj)
  return(freqtable)
}
