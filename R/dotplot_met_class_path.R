#'Plots dot plot pf PATHWAY and Metabolite class
#'@param kegg_es dataframe of enrichment score of pathways
#'@param met_path dataframe of pathways mapped to significant metabolites
#'@param classm metabolite class super_class, main_class or sub_class
#'@param xaxis font of x axis
#'@param yaxis font of y axis
#'@importFrom dplyr count
#'@importFrom dplyr rename
#'@importFrom ggplot2 scale_color_continuous
#'@importFrom ggplot2 guide_colorbar
#'@export
#'@examples
#'dotplot_met_class_path (met_path, kegg_es,"sub_class", xaxis, yaxis)

dotplot_met_class_path = function(met_path, kegg_es,classm,xaxis,yaxis){

met_pathdotplot = unique(met_path[,c('Metabolite','PATHWAY','log2Fold_change',classm)])
global_pathways = c('Metabolic pathways','Biosynthesis of secondary metabolites',
                    'Microbial metabolism in diverse environments',
                    'Carbon metabolism',
                    'Oxocarboxylic acid metabolism','Fatty acid metabolism',
                    'Biosynthesis of amino acids','Degradation of aromatic compounds')

##remove global pathways from plotting
met_pathtoplot2=met_pathdotplot[!met_pathdotplot$PATHWAY %in% global_pathways,]
freq_cal=met_pathtoplot2[,c("PATHWAY",classm,"Metabolite")]
freq_cal$PATHWAY = as.character(freq_cal$PATHWAY)
freq_cal[[classm]] = as.character(freq_cal[[classm]])
freq_cal=unique(freq_cal)

## get frequency of metabolites in each sub class and each pathway.
freq= dplyr::rename(count(freq_cal, PATHWAY,freq_cal[[classm]]), Freq = n)
names(freq)[2] =classm
# freq2=freq %>% arrange(Freq, sub_class, PATHWAY)
# freq2$sub_class = factor(freq2$sub_class, levels=unique(freq2$sub_class))
# freq2$PATHWAY = factor(freq2$PATHWAY, levels=unique(freq2$PATHWAY))
pathway_enrichment = merge(freq, kegg_es, by.x="PATHWAY","Pathway name")
significant_pathways = pathway_enrichment
 
if (nrow(significant_pathways)==0)
{significant_pathways = pathway_enrichment}else{significant_pathways =significant_pathways}
ggplot(significant_pathways, aes_(y=~PATHWAY, x=significant_pathways[[classm]])) +

  geom_point(aes(size=Freq, color=-log10(significant_pathways[["pathway_HG p-value"]])))+
  scale_color_continuous(low = "blue", high = "red",guide=guide_colorbar(reverse=TRUE))+
  theme_bw() + labs(color = "-log10 p-value HG")+

  theme(#panel.grid.major = element_blank(),
        axis.text.x=element_text(angle = 90, hjust = 1,size = xaxis),
          axis.text.y = element_text(size = yaxis),
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"))+
  xlab("Metabolite class")

}

