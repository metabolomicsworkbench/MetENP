#'counts and plots significant metabolite
#'@param df_metclass dataframe having significant metabolite with metabolite class as obtained from map_refmet_class.R
#'@param metclass sub_class, main_class or super_class
#'@param plot logical TRUE or FALSE, default true
#'@param thres_logfC log2 fold change threshold
#'@importFrom reshape2 dcast
#'@importFrom reshape2 melt
#'@importFrom dplyr filter
#'@importFrom dplyr mutate
#'@importFrom ggplot2 geom_bar
#'@importFrom ggplot2 element_line
#'@importFrom ggplot2 element_blank
#'@importFrom ggplot2 scale_size_manual
#'@importFrom ggplot2 labs
#'@export
#'@examples
#'count_changes = metcountplot(df_metclass=sig_metabolites_pubchem, metclass='sub_class', plotting=TRUE, thres_logfC = 0.5)


metcountplot <- function(df_metclass, metclass, plotting, thres_logfC)
{
  df_metclass=df_metclass[abs(df_metclass$log2Fold_change)>=thres_logfC,]
  df_metclass=df_metclass %>%
    mutate(Significant_Changes = (log2Fold_change>thres_logfC),
           neg = (log2Fold_change<thres_logfC))
df_metclass=as.data.frame(df_metclass)

  selected_refmet =df_metclass[, c(metclass, 'Significant_Changes','neg')]
  class_group = melt(selected_refmet, id.vars = c(metclass))
  class_group2=class_group[class_group$value == TRUE,]
  class_group2[[metclass]]=as.character(class_group2[[metclass]])
  dcast_grouping= dcast(class_group2, class_group2[,metclass] ~ variable, value.var = "value", fun.aggregate = sum)
  if (!("neg" %in% colnames(dcast_grouping))) {
    count_changes = filter(dcast_grouping, abs(Significant_Changes) > 0)
    names(count_changes)[1] = metclass

  }else{
  dcast_grouping$neg = dcast_grouping$neg* -1
  names(dcast_grouping)[1] = metclass
  neg_dcast_group = subset(dcast_grouping, select = c(metclass,'neg'))
  colnames(neg_dcast_group)[2] <- "Significant_Changes"
  dcast_grouping = subset(dcast_grouping, select = -neg)
  dcast_grouping2 <- bind_rows(dcast_grouping, neg_dcast_group)
  #count_changes = rbind(dcast_grouping, neg_dcast_group)

  count_changes = filter(dcast_grouping2, abs(Significant_Changes) > 0)
}

  if (plotting==TRUE)
  {
    for(i in 1:nrow(count_changes)){
      if (count_changes$Significant_Changes[i] >0){
        count_changes[i, "color"] = "increased metabolites"
      }else
      { count_changes[i, "color"] = "decreased metabolites"}
    }

 p<- ggplot(data=count_changes, aes(x=count_changes[[metclass]], y=Significant_Changes,fill=color)) +
    geom_bar(stat="identity",color="black", width=0.5) + theme_bw()+

    theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=1, size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size=14, face="bold"),
            axis.title.y = element_text(size=14, face="bold") ) +
    theme(panel.border = element_blank(), axis.line = element_line(color = "Black", size=1,
                                                                   lineend = "square"))+
    #theme(legend.position="top")+

    theme(plot.title = element_text(hjust = 0.5)) +
    scale_size_manual(values = c(1,1)) +
    #scale_fill_manual(values = c("blue","red")) +
	xlab(paste0("Metabolite class:", metclass)) +
    scale_color_manual(values = c("black", "black"))+ labs(fill = "sub-class")+
   geom_hline(yintercept = 0)


  return(list(sig_met_count= count_changes,plotimg= p))

  } else{return(count_changes)}
}
