#'counts and plots significant metabolite
#'@param df_metclass dataframe having significant metabolite with metabolite class as obtained from map_refmet_class.R
#'@param metclass sub_class, main_class or super_class
#'@param plotting logical TRUE or FALSE, default true
#'@param thres_logfC log2 fold change threshold
#'@param updown_fillcolor user-defined fill color for the bars for up/increased and down/decreased metabolites
#'@importFrom reshape2 dcast
#'@importFrom reshape2 melt
#'@importFrom dplyr filter
#'@importFrom dplyr mutate
#'@importFrom ggplot2 geom_bar
#'@importFrom ggplot2 element_line
#'@importFrom ggplot2 element_blank
#'@importFrom ggplot2 scale_size_manual
#'@importFrom ggplot2 scale_fill_manual
#'@importFrom ggplot2 labs
#'@export
#'@examples
#'count_changes = metcountplot(df_metclass=sig_metabolites_pubchem, metclass='sub_class', plotting=TRUE, thres_logfC = 0.5, updown_fillcolor=c("red", "green"))


#Original:metcountplot <- function(df_metclass, metclass, plotting, thres_logfC)
#Mano: 2022/06/18: adding arguments to let user decide bar fill color for up (increased) and down (decreased)
metcountplot <- function(df_metclass, metclass, plotting, thres_logfC, updown_fillcolor=c("red", "green"))
{
  df_metclass=df_metclass[abs(df_metclass$log2Fold_change)>=thres_logfC,]
  df_metclass=df_metclass %>%
    mutate(No.of_metabolites = (log2Fold_change>thres_logfC),
           neg = (log2Fold_change<thres_logfC))
df_metclass=as.data.frame(df_metclass)

  selected_refmet =df_metclass[, c(metclass, 'No.of_metabolites','neg')]
  class_group = melt(selected_refmet, id.vars = c(metclass))
  class_group2=class_group[class_group$value == TRUE,]
  class_group2[[metclass]]=as.character(class_group2[[metclass]])
  dcast_grouping= dcast(class_group2, class_group2[,metclass] ~ variable, value.var = "value", fun.aggregate = sum)
  if (!("neg" %in% colnames(dcast_grouping))) {
    count_changes = filter(dcast_grouping, abs(No.of_metabolites) > 0)
    names(count_changes)[1] = metclass

  }else{
  dcast_grouping$neg = dcast_grouping$neg* -1
  names(dcast_grouping)[1] = metclass
  neg_dcast_group = subset(dcast_grouping, select = c(metclass,'neg'))
  colnames(neg_dcast_group)[2] <- "No.of_metabolites"
  dcast_grouping = subset(dcast_grouping, select = -neg)
  dcast_grouping2 <- bind_rows(dcast_grouping, neg_dcast_group)
  #count_changes = rbind(dcast_grouping, neg_dcast_group)

  count_changes = filter(dcast_grouping2, abs(No.of_metabolites) > 0)
}

  if (plotting==TRUE)
  {
    for(i in 1:nrow(count_changes)){
      if (count_changes$No.of_metabolites[i] >0){
        count_changes[i, "color"] = "increased metabolites"
      }else
      { count_changes[i, "color"] = "decreased metabolites"}
    }

    # Mano: 2022/06/18: https://ggplot2.tidyverse.org/reference/scale_manual.html
    # Mano: 2022/06/18: used named colors
    fillcolor = c("increased metabolites" = updown_fillcolor[1], "decreased metabolites" = updown_fillcolor[2]);

 p<- ggplot(data=count_changes, aes(x=count_changes[[metclass]], y=No.of_metabolites,fill=color)) +
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
    geom_hline(yintercept = 0) +

    # Mano: 2022/06/18: had to preped ggplot:: before scale_fill_manual, else was saying function not found
    ggplot2::scale_fill_manual(values = fillcolor) ; # Mano: 2022/06/18: the user can decide the color, with default set to meaningful values


  return(list(sig_met_count= count_changes,plotimg= p))

  } else{return(count_changes)}
}
