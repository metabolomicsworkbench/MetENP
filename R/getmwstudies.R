#'Gets metabolomics data, metadata and metabolite info from Metabolomics Workbench using REST service
#'@param studyid from MW and typedata whether metadata (factors) or metabolic data 
#'@param typedata can be "factors" for metadata, "data" for metabolomics data and 'metabolites' for metabolite info
#'@importFrom tidyr unnest_wider
#'@examples
#' data = getmwstudies('ST001304', 'data')
#'@export

getmwstudies <- function(studyid, typedata)
{

if (typedata == "factors"){
path = paste0('https://www.metabolomicsworkbench.org/rest/study/study_id/', studyid,'/factors')
r <- httr::GET(url = path)
r <- httr::content(r, as = "text", encoding = "UTF-8")

df_met <- jsonlite::fromJSON(r)

extract_info <- lapply(df_met, '[', c(names(df_met[["1"]])))
dd=do.call(rbind, extract_info)
df_met1 = as.data.frame(dd)
  if (any(grepl("\\|", df_met1$factors))){
result <- data.frame(df_met1,do.call(rbind,stringr::str_split(df_met1$factors," \\| ")))
for (i in 1:length(stringr::str_split(df_met1$factors," \\| ")[[1]])){
 col_name = paste0("X",i) 
 name = strsplit(as.character(result[[col_name]]),":")[[1]][1]
# name1= strsplit(result$X1,":")[[1]][1]
# name2 = strsplit(result$X2,":")[[1]][1]
names(result)[names(result) == col_name] <- name
result[[name]]=gsub(paste0(name,":"),"",result[[name]])
#result[[name2]]=gsub(paste0(name2,":"),"",result[[name2]])
}
}else{
  name1= do.call(rbind,stringr::str_split(df_met1$factors,":"))[[1]][1]
  result=df_met1
  result[[name1]]=result$factors
  #result =  rename(df_met1, c('factors'=name1))
  result[[name1]]=gsub(paste0(name1,":"),"",result[[name1]])
}
}else if (typedata == "data"){
  path = paste0('https://www.metabolomicsworkbench.org/rest/study/study_id/', studyid,'/data')
  r <- httr::GET(url = path)
  r <- httr::content(r, as = "text", encoding = "UTF-8")
  
  df_met <- jsonlite::fromJSON(r)
  id=which.max(lapply(df_met, function(x) sum(lengths(x))))
  extract_info <- lapply(df_met, '[', c(names(df_met[[id]])))
  #extract_info <- lapply(df_met, '[', c(names(df_met[["1"]])))
  dd=do.call(rbind, extract_info)
  df_met1 = as.data.frame(dd)
  id2=which(is.na(colnames(df_met1)))
  names(df_met1)[id2] = names(df_met[[id]])[id2]
  data1=tidyr::unnest_wider(df_met1, DATA)
  result = dplyr::select(data1,-study_id, -units)

}else{
  path = paste0('https://www.metabolomicsworkbench.org/rest/study/study_id/', studyid,'/metabolites')
  r <- httr::GET(url = path)
  r <- httr::content(r, as = "text", encoding = "UTF-8")
  
  df_met <- jsonlite::fromJSON(r)
  id=which.max(lapply(df_met, function(x) sum(lengths(x))))
  extract_info <- lapply(df_met, '[', c(names(df_met[[id]])))
  dd=do.call(rbind, extract_info)
  df_met1 = as.data.frame(dd)
  id2=which(is.na(colnames(df_met1)))
  names(df_met1)[id2] = names(df_met[[id]])[id2]
  result = df_met1
}
  return(result)
}
