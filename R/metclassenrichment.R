#'Provides enrichment score of the selected metabolic class via hypergeometric score
#'@param df_metclass dataframe having significant metabolite with metabolite class and kegg ids
#'@param refmet_class dataframe having all the metabolites in the study with refmet classification
#'@param metclass Sub class, Main class or Super class
#'@param enrich_stats HG for hypergeometric score ## leaves room for further including of KS stats
#'@param  no number of significant metabolites that should be in a class, default 1, should increase to 3 or more for a
#'@param  debug flag to indicate if run in debug mode
#'better statistically significant score. Default is chosen as 1 so as to get information on all the significant metabolites without any filtering.
#'@export
#'@examples
#'metenrichment = metclassenrichment(df_metclass=sig_metabolites_kegg_id, refmet_class, metclass= "Sub class",enrich_stats="HG",no=1)
#'@importFrom stats phyper

metclassenrichment <- function(df_metclass,refmet_class, metclass,enrich_stats,no, debug = 0)
{

  #debug = 1; # Mano: 2022/06/21: added the line but decided to pass debug as argument with default value of 0

  # path= "https://www.metabolomicsworkbench.org/rest/refmet/all"
  # r <- GET(url = path)
  # r <- httr::content(r, as = "text", encoding = "UTF-8")
  # df <- fromJSON(r)
  # res=lapply(df, function(i) list(unlist(i, recursive = TRUE)))
  # extract_info <- lapply(df, '[', c(names(df[["Row1"]])))
  # dd=do.call(rbind, extract_info)
  # ref=as.data.frame(dd)
  tt <- table(factor(df_metclass[[metclass]]))

  df_metclass <- df_metclass[df_metclass[[metclass]] %in% names(tt[tt >= no]),]

  if (enrich_stats == "HG"){
    # Mano: 2022/06/21: add the column "HG p-value" before the for loop, initialize to NA
    df_metclass["HG p-value"] = NA; # Mano: 2022/06/21

    for(i in 1:nrow(df_metclass)){
    	classname <-as.character(df_metclass[i,metclass])


	#L = length(df_metclass) # Mano: 2022/06/21: This line seems wrong since length(data.frame) returns number of columns # it should be nrow(df_metclass)
	L = nrow(df_metclass);
	#L = nrow(ref[ref[[metclass]] == classname,])
    	N = nrow(refmet_class)
    	k = nrow(refmet_class[refmet_class[[metclass]] == classname,])
    	M = nrow(df_metclass[df_metclass[[metclass]] == classname,])
    	pp<-phyper(M-1,L, N-L,k, lower.tail=FALSE)

    	#print(pp)
    	if(debug>0){print(paste0("Name = ", df_metclass[i,"refmet_name"])); print(paste0(c("M", "L", "N", "k", "pp"), " = ", c(M, L, N, k, pp)));}

    	df_metclass[i,"HG p-value"] <- pp # Mano: 2022/06/21: Without adding the column in advance as Mano did above, this line is also problematic, as it addes a new columns for i = `1, then length(df_metclass) will be 1 more for i>1, that is why some sub_class-es has two HG p-values and the barplot has more than one element there (appears stacked)
    }
  } # Mano: 2022/06/21
  if (dim(na.omit(df_metclass))[1] == 0)
  {cat(blue ("metabolite enrichment could not be calculated. please check the number of metabolites, you may want to use less stringent number\n"))}
  return(df_metclass)
}
