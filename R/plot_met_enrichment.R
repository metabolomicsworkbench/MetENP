#'Plots metabolite enrichment score
#'@param df_metenrichment metenrichment score dataframe having metabolites, metabolite class and their enrichment score obtained
#'from metclassenrichment
#'@param metclass = sub_class, main_class or super_class
#'@param enrich_stats HG for hypergeometric score ## leaves room for further including of KS stats
#'@param  no number of significant metabolites that should be in a class, default 1, should increase to 3 or more for a
#'better statistically significant score. Default is chosen as 1 so as to get information on all the significant metabolites without any filtering.
#'@importFrom ggplot2 coord_flip
#'@export
#'@examples
#'plot_met_enrichment(metenrichment, "sub_class","HG", no=1)

plot_met_enrichment <- function(df_metenrichment, metclass,enrich_stats, no)
{
  df_metenrichment[[metclass]]=factor(as.character(df_metenrichment[[metclass]]))
  tt <- table(unlist(df_metenrichment[[metclass]]))
  df2 <- subset(df_metenrichment, sub_class %in% names(tt[tt >= no]))

  # Mano: 2022/06/19: streamline, if use something many times, use an intermediate variable
  enrich_stats_pvalue_colname = paste0(enrich_stats, " p-value");
  ylab_value = paste("-log10",enrich_stats_pvalue_colname); 

  #metclass_stats = df2[,c(metclass,paste0(enrich_stats, " p-value"))] %>% distinct() # Mano commented
  #metclass_stats[[enrich_stats]] = -log10(metclass_stats[[paste0(enrich_stats, " p-value")]]) # Mano commented
  metclass_stats = df2[,c(metclass, enrich_stats_pvalue_colname)] %>% distinct()
  metclass_stats[[enrich_stats]] = -log10(metclass_stats[[enrich_stats_pvalue_colname]])

  ggplot(data=metclass_stats, aes(x=metclass_stats[[metclass]], y=metclass_stats[[enrich_stats]],fill=metclass_stats[[enrich_stats]])) +
    geom_bar(stat="identity",color="black", width=0.5) + theme_bw()+

    theme(axis.text.x = element_text(angle = 30, hjust=1, vjust=1, size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size=14, face="bold"),
          axis.title.y = element_text(size=14, face="bold")) +
    theme(panel.border = element_blank(), axis.line = element_line(color = "Black", size=1,
                                                                   lineend = "square"))+
    #theme(legend.position="top")+

    theme(plot.title = element_text(hjust = 0.5)) +
    scale_size_manual(values = c(1,1)) +
    #scale_fill_manual(values = c("blue","red")) +
    #scale_color_manual(values = c("black", "black")) + labs(fill = "-log10(p value)")+ # Mano commented
    scale_color_manual(values = c("black", "black")) + labs(fill = ylab_value)+
    coord_flip()+
    #ylab(paste("-log10",enrich_stats_pvalue_colname)) + xlab(paste("Class:",metclass)) # Mano commented
    ylab(ylab_value) + xlab(paste("Class:",metclass))
}
