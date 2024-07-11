#' Provides log2fold change and p value on metabolomics data using t test.
#'@description Provides log2fold change and p value on metabolomics data using t test. One can then calculate significant metabolites by using the dataframe generated from this function.
#'@param metabolomics_data metabolomics data associated to refmet class
#'@param met_col column with metabolite names
#'@param analysis_type type of analysis ex-GCMS, HILIAC positive ion mode.
#'@param metadata Metadata
#'@param factor1 first independent variable
#'@param factor2 second independent variable
#'@param factor_col column name of the independent variables
#'@param sample_col the column name having samples
#'@param p_adjust Method for p value adjustment, i.e. "fdr"
#'@param normalization method for normalization : "half_of_min": where the NAs are replaced by half of min values in the data
#' "remove_NAs": where Cols with NAs values are removed
#' "50percent": where cols with more than 50% NAs values are removed
#'@importFrom stats p.adjust
#'@importFrom stats aggregate
#'@importFrom stats t.test
#'@importFrom plyr ldply
#'@importFrom dplyr select
#'@importFrom magrittr %>%
#'@export
#'@examples
#'stats_metabolites = significant_met(metabolomics_data=data,met_col="metabolite_name",
#'analysis_type=c('HILIC NEGATIVE ION MODE','HILIC POSITIVE ION MODE'),
#'metadata=metadata, factor1='No shield', factor2='plus shield',
#'factor_col='treatment',sample_col='local_sample_id', p_adjust='fdr',normalization="half_of_min")

significant_met <- function(metabolomics_data,met_col,analysis_type,
                            metadata,normalization, factor1, factor2, factor_col,sample_col, p_adjust)
{
  combinedresult=list()
  for (a in 1:length(analysis_type)){
  metabolomics_data$refmet_name= gsub("\t", "", metabolomics_data$refmet_name)
  metabolomics_data$metabolite_name=as.character(metabolomics_data$metabolite_name)

  analysis_selected = metabolomics_data[metabolomics_data[["analysis_summary"]] == analysis_type[a],]


  analysis_selected=analysis_selected %>% dplyr::distinct(metabolite_name, .keep_all = TRUE)



if ("metabolite_id" %in% colnames(analysis_selected)){
  refmet_id = analysis_selected[,c('metabolite_name','metabolite_id','refmet_name')]
analysis_selected=analysis_selected %>% dplyr::select(-metabolite_id, -refmet_name) ### numeric data only
analysis_selected=analysis_selected %>% dplyr::select(-analysis_id) ### numeric data only
analysis_selected=analysis_selected %>% dplyr::select(-analysis_summary)
### numeric data only
}else{
  refmet_id = analysis_selected[,c('metabolite_name','refmet_name')]
  analysis_selected=analysis_selected %>% dplyr::select(-refmet_name,-analysis_summary)

}
  met_col_name=analysis_selected[[met_col]]

  #analysis_selected[[met_col]] <- NULL
  row.names(analysis_selected) <- met_col_name

#metabolomics_data=metabolomics_data[complete.cases(metabolomics_data), ]

metabolomics_data_transposed <- as.data.frame(t(analysis_selected))
metabolomics_data_transposed <- metabolomics_data_transposed [-1,]
metadata[, factor_col]=as.character(metadata[, factor_col])
Factor = factor(metadata[, factor_col])
metabolomics_data_transposed$Factor<-Factor[match(rownames(metabolomics_data_transposed), metadata[,sample_col])]

metabolomics_data_subset=metabolomics_data_transposed[metabolomics_data_transposed$Factor %in% factor1 | metabolomics_data_transposed$Factor %in% factor2,]
#metabolomics_data_subset<- metabolomics_data_transposed[which(metabolomics_data_transposed$Factor == factor1|metabolomics_data_transposed$Factor == factor2),]

metabolomics_data_subset2 <- metabolomics_data_subset%>% dplyr::select(-Factor) ### numeric data only
for (i in 1:ncol(metabolomics_data_subset2)){

  metabolomics_data_subset2[[i]]=suppressWarnings(as.numeric(paste(metabolomics_data_subset2[[i]])))
}

if (normalization=="half_of_min") ### for this you need to remove the columns which are only in one factor
{
#lowst_per= min(table(Factor))/length(Factor) ## lowest percentage of factors, to be sure we do not leave meatbolites only in one factor
  min_val = min(apply(metabolomics_data_subset2,2,min, na.rm=T))
  # metabolomics_data_subset2=metabolomics_data_subset2 %>%
    # purrr::discard(~sum(is.na(.x))/length(.x)* 100 >=lowst_per)
  metabolomics_data_subset2[is.na(metabolomics_data_subset2)] <- min_val/2
}else if (normalization == "remove_NAs")
{
  metabolomics_data_subset2=metabolomics_data_subset2 %>%
      dplyr::select_if(~ !any(is.na(.)))
}else if (normalization == "50percent")
{
lowst_per= min(table(Factor))/length(Factor) ## lowest percentage of factors, to be sure we do not leave meatbolites only in one factor

  metabolomics_data_subset2=metabolomics_data_subset2 %>%
    purrr::discard(~sum(is.na(.x))/length(.x)* 100 >=lowst_per*100)
  min_val = min(apply(metabolomics_data_subset2,2,min, na.rm=T))
  # metabolomics_data_subset2=metabolomics_data_subset2 %>%
  # purrr::discard(~sum(is.na(.x))/length(.x)* 100 >=lowst_per)
  metabolomics_data_subset2[is.na(metabolomics_data_subset2)] <- min_val/2
}

  metabolomics_data_normalized <- metabolomics_data_subset2


### add factors column again
metabolomics_data_normalized$Factor = metabolomics_data_subset$Factor

cdatamod<- colnames(metabolomics_data_normalized[, which(names(metabolomics_data_normalized) != "Factor")]) ### colnames of just metabolites
fac_col = metabolomics_data_normalized[, "Factor"]

fac_col=factor(metabolomics_data_normalized$Factor, levels = c(factor1, factor2))

if (length(factor1)>1|length(factor2) >1){
  y <- matrix(nrow=ncol(metabolomics_data_normalized)-1,ncol=6)
  num=ncol(metabolomics_data_normalized)-1
  for(i in 1:num){
    d<- metabolomics_data_normalized[metabolomics_data_normalized$Factor %in% factor2,i]
    n<- metabolomics_data_normalized[metabolomics_data_normalized$Factor %in% factor1,i]
    #tt <- t.test(n,d,alternative="two.sided", paired=FALSE) # Mano: original till 2024/07/10
    tt <- t.test(n,d,alternative="two.sided") # Mano: 2024/07/10: removed arg paired due to error
    ratio <- (mean(d)/mean(n))
    log2Fold_change <- log2(ratio)

    j=i
    y[j,1] <- as.character(colnames(metabolomics_data_transposed)[i])
    y[j,2] <- tt$statistic
    y[j,3] <- tt$p.value

    y[j,5]=ratio
    y[j,6]=log2Fold_change
  }
  y[, 4] <- p.adjust(y[, 3], method=p_adjust)
results = as.data.frame(y)
names(results)[1]="Metabolite"
names(results)[2]="t_value"
names(results)[3]="pval"
names(results)[6]="log2Fold_change"
names(results)[5]="Fold_change"
names(results)[4]="padj"
}else{

res <- ldply(
  cdatamod, function(Metabolite) {
        #t_val =t.test(metabolomics_data_normalized[[Metabolite]]~fac_col,alternative="two.sided", paired=FALSE)$statistic # Mano: original till 2024/07/10
        #p_val = t.test(metabolomics_data_normalized[[Metabolite]]~ fac_col,alternative="two.sided", paired=FALSE)$p.value # Mano: original till 2024/07/10
        t_val =t.test(metabolomics_data_normalized[[Metabolite]]~fac_col,alternative="two.sided")$statistic # Mano: 2024/07/10: removed arg paired due to error
        p_val = t.test(metabolomics_data_normalized[[Metabolite]]~ fac_col,alternative="two.sided")$p.value # Mano: 2024/07/10: removed arg paired due to error
        return(data.frame(Metabolite=Metabolite, t_value=t_val, pval = p_val))
    })


res$padj = p.adjust(res$pval, method = p_adjust)

if (typeof(metabolomics_data_subset$Factor)=='list'){
  mean_val=aggregate(metabolomics_data_subset2, list(as.character(metabolomics_data_subset$Factor)), mean, na.rm=TRUE)
}else{

mean_val = aggregate(metabolomics_data_subset2, list(metabolomics_data_subset$Factor), mean, na.rm=TRUE)}
Means_transposed = as.data.frame(t(mean_val))
names(Means_transposed) <- paste(as.character(unlist(Means_transposed[1,])), "mean",sep="_")
Means_transposed = Means_transposed[-1,]
Means_transposed$Fold_change=as.numeric(as.character(Means_transposed[,paste(factor2,"mean", sep="_")])) / as.numeric(as.character(Means_transposed[,paste(factor1,"mean", sep="_")]))
Means_transposed$log2Fold_change=log2(as.numeric(as.character(Means_transposed[,3])))
Means_transposed$Metabolite = row.names(Means_transposed)
#row.names(Means_transposed)=NULL
#Means_transposed2= data.frame(Means_transposed, Means_transposed, row.names = NULL)
results = merge(Means_transposed, res, by = "Metabolite", all.y = TRUE)
results %>% dplyr::select(-(paste(factor1, "mean", sep = "_")),-(paste(factor2, "mean", sep = "_")))
}
#results =results[is.finite(results[,"log2Fold_change"]),]
results=merge(results, refmet_id, by.x = "Metabolite", by.y="metabolite_name")
refmet_class2 = dplyr::select(analysis_selected, c(metabolite_name, super_class, main_class, sub_class, formula))
results2= merge(results, refmet_class2, by.y='metabolite_name', by.x='Metabolite')
combinedresult[[a]]=results2
  }
  combinedresult2=dplyr::bind_rows(combinedresult)
  combinedresult2[,"pval"]= as.numeric(paste(combinedresult2[,"pval"]))
  combinedresult2[,"log2Fold_change"]= as.numeric(paste(combinedresult2[,"log2Fold_change"]))
  combinedresult2=unique(combinedresult2)
return(combinedresult2)
  }

