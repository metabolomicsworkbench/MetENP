#'Separates user data into metadata and metabolomics data with metabolite measurements
#'@param example user data file
#'@param typedata type of data desired from the user data uploaded: metadata or metabolomics data 
#'@param sample_col TRUE or FALSE, TRUE if sample names are in the first column, FALSE if metabolites names are in the 
#' first column.
#'@export
#'@examples
#'separate_data('example.txt', "metadata",FALSE)


separate_data<- function(example,typedata,sample_col){

if (getExtension(example) == "txt"){
example_data = read.csv(example, sep='\t', header=TRUE,check.names = FALSE)
}else if (getExtension(example) == "csv"){
example_data = read.csv(example, sep=',', header=TRUE,check.names = FALSE)
}
if (sample_col == TRUE){
row.names(example_data) = example_data[[1]]
example_data=dplyr::select(example_data, -names(example_data[1]))
example_data2= as.data.frame(t(example_data))
example_data2=cbind(Samples = rownames(example_data2), example_data2)
rownames(example_data2) <- 1:nrow(example_data2)
}else{example_data2=example_data}
if (typedata=="metadata"){
metadata= example_data2[1,1:(ncol(example_data2))]
metadata2=as.data.frame(t(metadata))
metadata2= cbind(local_sample_id = rownames(metadata2 ), metadata2 )
#names(metadata2)[2] = 'combined_factors'
names(metadata2)[2] = as.character(metadata2[[1,2]])
rownames(metadata2)=NULL
metadata2=metadata2[-1,]

if (any(grepl("\\|", metadata2[[2]]))){
names(metadata2)[2] = 'combined_factors'
  result <- data.frame(metadata2,do.call(rbind,stringr::str_split(metadata2$combined_factors," \\| ")))
  for (i in 1:length(stringr::str_split(metadata2$combined_factors," \\| ")[[1]])){
    col_name = paste0("X",i) 
    name = strsplit(as.character(result[[col_name]]),":")[[1]][1]
    # name1= strsplit(result$X1,":")[[1]][1]
    # name2 = strsplit(result$X2,":")[[1]][1]
    names(result)[names(result) == col_name] <- name
    result[[name]]=gsub(paste0(name,":"),"",result[[name]])
    #result[[name2]]=gsub(paste0(name2,":"),"",result[[name2]])
  }
  }else{result=metadata2}
  
  }else if (typedata=="data"){

metabolomics_data=example_data2[-1,]
result= dplyr::rename(metabolomics_data, metabolite_name=Samples)
  }
return(result)}
