#'Plots Volcano plot
#'@param stats_metabolites dataframe obtained by running sig_metabolites.R
#'@param thres_pval p value threshold
#'@param thres_log2foldchange log2fold change threshold
#'@examples
#'plot_volcano(stats_metabolites, thres_pval= 0.05,thres_log2foldchange = 0.5)
#'@export
#'@importFrom ggplot2 ggplot
#'@importFrom ggplot2 geom_point
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 element_text
#'@importFrom ggrepel geom_text_repel
#'@importFrom magrittr %>%

plot_volcano <- function(stats_metabolites, thres_pval, thres_log2foldchange, text_nooverlap)
{
# threshold_OE <- stats_metabolites$pval < thres_pval & abs(stats_metabolites$log2Fold_change) > thres_log2foldchange
# stats_metabolites$threshold <- threshold_OE
stats_metabolites$log2Fold_change=  as.numeric(paste(stats_metabolites$log2Fold_change))
stats_metabolites$pval = as.numeric(paste(stats_metabolites$pval))


stats_metabolites <- stats_metabolites %>%
  dplyr::mutate(threshold = factor(dplyr::case_when(stats_metabolites$log2Fold_change > thres_log2foldchange & stats_metabolites$pval < thres_pval ~ "cond1",
                                       stats_metabolites$log2Fold_change < -thres_log2foldchange & stats_metabolites$pval < thres_pval ~ "cond2",
                                      abs(stats_metabolites$log2Fold_change) < thres_log2foldchange & stats_metabolites$pval < thres_pval ~ "cond4",
                                                                         TRUE ~ "cond3")))
stats_metabolites_ordered <- stats_metabolites[order(stats_metabolites$pval), ]
stats_metabolites_ordered$Metabolite = as.character(stats_metabolites_ordered$Metabolite)

## Create a column to indicate which met to label

stats_metabolites_ordered=stats_metabolites_ordered %>%
  dplyr::mutate(metlabels = factor(dplyr::case_when(stats_metabolites_ordered$threshold == "cond1" ~ T,
                                      stats_metabolites_ordered$threshold == "cond2" ~ T,
                                      stats_metabolites_ordered$threshold == "cond4" ~ T,
                                      TRUE ~ F)))

p=ggplot(stats_metabolites_ordered) +
  geom_point(aes(x=log2Fold_change, y=-log10(pval),colour=threshold)) +

  theme_bw()+
  ggtitle("Volcano Plot") +
  xlab("log2 fold change") +
  ylab("-log10 p-value") +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5,size=14, face="bold"),
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"),
        axis.text.x=element_text(size = 10),
        axis.text.y = element_text(size = 10))+
  geom_hline(yintercept = -log10(thres_pval),colour="#990000", linetype="dashed")+
#geom_vline(xintercept = -as.numeric(thres_log2foldchange),colour="#990000", linetype="dashed")
geom_vline(xintercept = c(-thres_log2foldchange, thres_log2foldchange),colour="#990000", linetype="dashed")+
  scale_color_manual(name = "Threshold",
                     values = c("cond1" = "red", "cond2" = "green", "cond3" = "black", "cond4"= "orange"))
if (text_nooverlap){
  p+geom_text_repel(aes(x = log2Fold_change, y = -log10(pval), label = ifelse(metlabels == T, Metabolite,"")), size = 3)

}else{p+geom_text(aes(x = log2Fold_change, y = -log10(pval), label = ifelse(metlabels == T, Metabolite,"")), size = 3)}




}

