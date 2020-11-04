#'Assigns kegg id to the significant metabolites
#'@param sig_metabolites dataframe having significant metabolite with metabolite class
#'@export
#'@examples
#'sig_metabolites_kegg_id= map_keggid(sig_metabolites)


map_keggid = function(sig_metabolites){

  datalist = list()
  for (i in 1:nrow(sig_metabolites))
  {
    path= paste0('https://www.metabolomicsworkbench.org/databases/refmet/name_to_refmet_kegg.php?metabolite_name=',sig_metabolites[["refmet_name"]][i])

    if (grepl("\\s*", path)){
      path=URLencode(path)
    }
    r <- httr::GET(url = path)
    #r <- content(r, as = "text", encoding = "UTF-8")
    df=read.csv(textConnection(httr::content(r, 'text')), sep="\t", header = FALSE)
    #df_t=as.data.frame(t(df))
    names(df) <- as.character(unlist(df[1,]))

    df <- df[-1, ]
    df=dplyr::rename(df, "refmet_name"="Standardized name")
    #df[['refmet_name']]=df[,"Standardized name"][i]
    datalist[[i]]  = df
  }
  refmet_kegg_match = do.call(rbind, datalist)
  refmet_kegg_match=dplyr::select(refmet_kegg_match,!c('In RefMet'))
  refmet_kegg_match2 =merge(refmet_kegg_match, sig_metabolites, by = 'refmet_name')
  refmet_kegg_match2=unique(refmet_kegg_match2)
  refmet_kegg_match2= dplyr::select(refmet_kegg_match2, c(-Formula, -'Input name', -'Sub class', -'Main class', -'Super class'))
  return(refmet_kegg_match2)

      }
