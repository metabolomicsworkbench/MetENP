#' Associates refmet class to the metabolites
#' @param metdata metabolomics data retrieved or uploaded by the user
#'@export
#'@importFrom dplyr select
#'@importFrom sjmisc is_empty
#'@examples
#'refmet_class= convert_refmet(metdata)

convert_refmet = function(metdata){

  datalist2=list()
  for (i in 1:nrow(metdata)){
  path= paste0("https://www.metabolomicsworkbench.org/rest/refmet/match/",as.character(metdata[["metabolite_name"]][i]))
  path=URLencode(path)
  r <- httr::GET(url = path)
  r <- httr::content(r, as = "text", encoding = "UTF-8")
  df1 <- jsonlite::fromJSON(r)
  if (!sjmisc::is_empty(df1)){
  if (typeof(df1[[names(df1)[1]]]) == 'character'){
    res=lapply(df1, function(i) (unlist(i, recursive = TRUE)))
    #res2=res
  }else{
    #res=t(as.data.frame((lapply(df, unlist))))
    res= bind_rows(df1)}
  }else{res=lapply(df1, function(i) (unlist(i, recursive = TRUE)))}


  datalist2[[i]]  = res
  datalist2[[i]][["metabolite_name"]] = as.character(metdata[["metabolite_name"]][i])
  }
  #extract_info <- lapply(df, '[', c(names(df[["Row1"]])))
  final=bind_rows(datalist2)
  ## Sumana: change 3/31/2022, class_index deprecated in MW, both brances seem to do the same thing
  if ("metabolite_id" %in% colnames(metdata)){
    ##final2=unique(dplyr::select(final,c( -exactmass, -class_index)))
    final2=unique(dplyr::select(final,c( -exactmass)))
  }else{
    ##final2=unique(dplyr::select(final,c( -exactmass, -class_index)))
    final2=unique(dplyr::select(final,c( -exactmass)))
  }
  final2=final2[!duplicated(final2$metabolite_name),]
  #final3=final2[,c("metabolite_name","refmet_name")]
  metdata$metabolite_name=as.character(metdata$metabolite_name)
  metdata=metdata[order(metdata$metabolite_name),]
  data_refmet=merge(metdata, final2, by="metabolite_name")#cbind(metdata, final2)

  ### get all the row of duplicates, just duplicated gives one row
  tes=data_refmet[data_refmet$metabolite_name %in% data_refmet$metabolite_name[duplicated(data_refmet$metabolite_name)],]

  #refmet_data2= subset(data_refmet, !(data_refmet$refmet_name.y %in% tes$refmet_name.y[tes$refmet_name.x !=tes$refmet_name.y]))

  if ("metabolite_id" %in% colnames(metdata)){
  #refmet_data2=refmet_data2[, !duplicated(colnames(data_refmet))]
  refmet_data2=dplyr::select(data_refmet,c( -refmet_name.y))
  refmet_data2= dplyr::rename(refmet_data2, refmet_name = refmet_name.x)
  refmet_data2 =unique(refmet_data2)
  }else{refmet_data2 =unique(data_refmet)}
return(refmet_data2)
}
