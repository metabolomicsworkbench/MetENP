#'Plots heatmap
#'@param met_path dataframe of pathways mapped to significant metabolites
#'@param shorten_name default TRUE if metabolite names are too big
#'@param refmet_name default false, set to true if you want to plot according to refmet names
#'@export
#'@importFrom ggplot2 geom_tile
#'@importFrom ggplot2 scale_fill_continuous
#'@examples
#'plot_heatmap(met_path)

plot_heatmap = function(met_path, shorten_name=TRUE,refmet_name=FALSE)
{
met_pathtoplot = unique(met_path[,c('Metabolite','refmet_name','PATHWAY','log2Fold_change')])
global_pathways = c('Metabolic pathways','Biosynthesis of secondary metabolites',
                    'Microbial metabolism in diverse environments',
                    'Carbon metabolism',
                    'Oxocarboxylic acid metabolism','Fatty acid metabolism',
                    'Biosynthesis of amino acids','Degradation of aromatic compounds')

##remove global pathways from plotting
met_pathtoplot2=met_pathtoplot[!met_pathtoplot$PATHWAY %in% global_pathways,]

tidy_name <- function(name, n_char) {
  ifelse(nchar(name) > (n_char - 2),
    {substr(name, 1, n_char) %>% paste0(., "..")},
    name)
}

if (shorten_name == TRUE){

met_pathtoplot2$Metabolite = tidy_name(as.character(met_pathtoplot2$Metabolite),25)
}else{met_pathtoplot2$Metabolite=met_pathtoplot2$Metabolite}



if (refmet_name==TRUE){
p=ggplot(met_pathtoplot2, aes_(~refmet_name, ~PATHWAY))
}else{p=ggplot(met_pathtoplot2, aes_(~Metabolite, ~PATHWAY)) }

p + geom_tile(aes_(fill = ~log2Fold_change), color = "white") +
  scale_fill_continuous(low="blue", high="red", name = "fold change")+
  theme_minimal() +
  xlab("Metabolite name") +
  theme(panel.grid.major = element_blank(),
        axis.text.x=element_text(angle = 90, hjust = 1,size = 8),
          axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"))

}
