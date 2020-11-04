#'Performs anova analysis
#'Performs an anova across all response variables. The function takes 2
#'independent variables.
#'@param metabolomics_data metabolomics data
#'@param metadata Metadata
#'@param met_col column with metabolite names
#'@param analysis_type type of analysis ex-GCMS, HILIAC positive ion mode.
#'@param sample_col the column name having samples
#'@param normalization normalization methods any of the three: half_of_min, remove_NAs and 50percent
#' 1) half_of_min: where the NAs are replaced by half of min values in the data
#'2) remove_NAs: where Cols with NAs values are removed
#'3) 50percent: where cols with more than 50% NAs values are removed
#'@param factor_col1 column name of first independent variable
#'@param factor_col2 column name of second independent variable
#'@param p_adjust Method for p value adjustment, i.e. "fdr"
#'@importFrom plyr llply
#'@importFrom dplyr distinct
#'@importFrom dplyr bind_rows
#'@importFrom stats p.adjust
#'@importFrom stats anova
#'@importFrom stats lm
#'@examples
#'anova_result = anova_ana(metabolomics_data=data,met_col="metabolite_name",analysis_type=c('Phospholipids, Chol. esters and Diacylglycerols','Sphingolipids'), metadata=metadata,normalization="50percent", factor_col1="TreatmentGroup",factor_col2="SamplingTimePoint",p_adjust="fdr")
#'@export


anova_ana = function(metabolomics_data,met_col,analysis_type,sample_col,
                     metadata,normalization, factor_col1,factor_col2,p_adjust)
{
  combinedresult=list()
  for (a in 1:length(analysis_type)){
    metabolomics_data$refmet_name= gsub("\t", "", metabolomics_data$refmet_name)
    metabolomics_data$metabolite_name=as.character(metabolomics_data$metabolite_name)

  analysis_selected = metabolomics_data[metabolomics_data[["analysis_summary"]] == analysis_type[a],]
  #refmet_id = analysis_selected[,c('metabolite_name','metabolite_id','refmet_name')]
  analysis_selected=analysis_selected %>% distinct(metabolite_name, .keep_all = TRUE)


  if ("metabolite_id" %in% colnames(analysis_selected)){
    refmet_id = analysis_selected[,c('metabolite_name','metabolite_id','refmet_name')]
  #analysis_selected[[met_col]] <- NULL
  analysis_selected=analysis_selected %>% dplyr::select(-metabolite_id, -refmet_name) ### numeric data only
  analysis_selected=analysis_selected %>% dplyr::select(-analysis_id) ### numeric data only
  analysis_selected=analysis_selected %>% dplyr::select(-analysis_summary) ### numeric data only
  }else{
    refmet_id = analysis_selected[,c('metabolite_name','refmet_name')]
    analysis_selected=analysis_selected %>% dplyr::select(-refmet_name,-analysis_summary)

  }
  row.names(analysis_selected) <- analysis_selected[[met_col]]
  analysis_selected[[met_col]] <- NULL


  metabolomics_data_transposed <- as.data.frame(t(analysis_selected))


  for (i in 1:ncol(metabolomics_data_transposed)){

    metabolomics_data_transposed[[i]]=suppressWarnings(as.numeric(paste(metabolomics_data_transposed[[i]])))
  }
  
var2 = metadata[, factor_col2]
var1 = metadata[, factor_col1]

  if (normalization=="half_of_min")
  {
    min_val = min(apply(metabolomics_data_transposed,2,min, na.rm=T))
    metabolomics_data_transposed[is.na(metabolomics_data_transposed)] <- min_val/2
  }else if (normalization == "remove_NAs")
  {
    metabolomics_data_transposed=metabolomics_data_transposed %>%
      select_if(~ !any(is.na(.)))
  }else if (normalization == "50percent")
  {
    lowst_per= min(min(table(var1))/length(var1),min(table(var1))/length(var1))
    metabolomics_data_transposed=metabolomics_data_transposed %>%
      purrr::discard(~sum(is.na(.x))/length(.x)* 100 >=lowst_per)
  }

  metabolite_col = colnames(metabolomics_data_transposed)
 


metabolomics_data_transposed[,factor_col1]<-var1[match(rownames(metabolomics_data_transposed), metadata[,sample_col])]
metabolomics_data_transposed[,factor_col2]<-var2[match(rownames(metabolomics_data_transposed), metadata[,sample_col])]



results <- llply(
  metabolite_col, function(x) {
    models <- lm(metabolomics_data_transposed[[x]] ~ var1 + var2 + var1 * var2, data = metabolomics_data_transposed)
  #posthoc_test=TukeyHSD(models)
    })
names(results) <- as.character(metabolite_col)
#results <- lapply(results, summary)
results <- lapply(results, anova)
results <- sapply(results, cbind)
results <- t(results)
#results1 <- sapply(results1, cbind)
#results1 <- t(results1)
results_df<- as.data.frame(results[,"Pr(>F)"])

results_df<- as.data.frame(t(results_df))
results_df <- results_df[,1:3]
colnames(results_df)[1] <- factor_col1
colnames(results_df)[1] <- paste(colnames(results_df)[1], "pval", sep = ".")
colnames(results_df)[2] <- factor_col2
colnames(results_df)[2] <- paste(colnames(results_df)[2], "pval", sep = ".")
colnames(results_df)[3] <- "Interaction.pval"
results_df$padj = p.adjust(results_df[,1], method = p_adjust)
colnames(results_df)[4] <- paste(colnames(results_df)[4], factor_col1, sep = ".")
results_df$padj = p.adjust(results_df[,2], method = p_adjust)
colnames(results_df)[5] <- paste(colnames(results_df)[5], factor_col2, sep = ".")
results_df$padj = p.adjust(results_df[,3], method = p_adjust)
colnames(results_df)[6] <- paste(colnames(results_df)[6], "Interaction", sep = ".")
#row.names(results_df) = as.character(metabolite_col)
results_df[,"metabolic_name"]=as.character(metabolite_col)
combinedresult[[a]]=results_df

  }
  combinedresult2=bind_rows(combinedresult)
  combinedresult2=combinedresult2[, c(7, 1, 2, 3,4,5,6)]
return(combinedresult2)
}
