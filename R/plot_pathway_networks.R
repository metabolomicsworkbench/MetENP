#' Plots 1st degree pathway network of only pathways having significant pathway score
#'@param met_path dataframe of pathways mapped to significant metabolites
#'@param kegg_es dataframe of enrichment score of pathways
#'@param colorEdge default TRUE to color the edges by pathways
#'@importFrom tidygraph as_tbl_graph
#'@importFrom tidygraph activate
#'@importFrom tidygraph pull
#'@importFrom igraph E
#'@importFrom igraph V
#'@importFrom ggraph geom_edge_diagonal
#'@importFrom ggraph geom_edge_diagonal0
#'@importFrom ggraph ggraph
#'@importFrom ggraph geom_node_point
#'@importFrom ggraph geom_node_text
#'@importFrom ggplot2 theme_minimal
#'@importFrom ggplot2 theme_set
#'@importFrom ggplot2 aes_
#'@importFrom ggplot2 scale_colour_gradient2
#'@importFrom ggplot2 scale_shape_identity
#'@export
#'@examples
#'plot_pathway_networks (met_path,kegg_es, TRUE)

plot_pathway_networks = function(met_path, kegg_es, colorEdge)
{
  pathway_enrichment = merge(met_path, kegg_es, by.x="PATHWAY",by.y="Pathway name")

  #significant_pathways = pathway_enrichment[which(pathway_enrichment$pathway_HG < 0.05),]

significant_pathways = pathway_enrichment
met_pathtoplot = unique(significant_pathways[,c('Metabolite','PATHWAY')])
global_pathways = c('Metabolic pathways','Biosynthesis of secondary metabolites',
                    'Microbial metabolism in diverse environments',
                    'Carbon metabolism',
                    'Oxocarboxylic acid metabolism','Fatty acid metabolism',
                    'Biosynthesis of amino acids','Degradation of aromatic compounds')

##remove global pathways from plotting
met_pathtoplot2=met_pathtoplot[!met_pathtoplot$PATHWAY %in% global_pathways,]


met_path_fc = unique(met_path[,c('Metabolite','log2Fold_change')])

graph_routes <- as_tbl_graph(met_pathtoplot2)

# node_names <- graph_routes %>%
  # activate(nodes) %>%
  # pull(name)


thm <- theme_minimal() +
  theme(
    legend.position = "right",
    axis.title = element_blank(),
    axis.text = element_blank(),
    panel.grid = element_blank(),
    panel.grid.major = element_blank(),
  )

theme_set(thm)


if (colorEdge) {
  igraph::E(graph_routes)$category <- factor(met_pathtoplot2$PATHWAY)
  #edge_color <- geom_edge_diagonal(aes(colour=factor(met_pathtoplot2$PATHWAY)), alpha=.8,show.legend=FALSE)
  #edge_color <- geom_edge_diagonal(alpha=.8, colour='gray')
  edge_color <- geom_edge_diagonal(aes_(colour = ~category), alpha=.8,show.legend=FALSE)
} else {
  edge_color <- geom_edge_diagonal(alpha=.8, colour='gray')
}

fc= met_path_fc$log2Fold_change[which(met_path_fc$Metabolite %in%  V(graph_routes)$name)]
igraph::V(graph_routes)$color = NA
igraph::V(graph_routes)$color[which(V(graph_routes)$name %in% met_path_fc$Metabolite)] =fc
igraph::V(graph_routes)$size = 5 ## set min
tt = as.data.frame(table(met_pathtoplot2$PATHWAY))
#tt$Freq=tt$Freq*3
tt$Freq[tt$Freq == 1] = 1.5 ### inc the size of node to min

for (i in 1:length(igraph::V(graph_routes)$name))
{
  for (j in 1:nrow(tt))
  {
    if (igraph::V(graph_routes)$name[i] == tt$Var1[j])
    {
      igraph::V(graph_routes)$size[i] = tt$Freq[j]
    }
  }
}

igraph::V(graph_routes)$shape=16
pthway= as.data.frame(unique(met_pathtoplot2$PATHWAY))
names(pthway)[1]='pathway'
for (i in 1:length(igraph::V(graph_routes)$name))
{
  for (j in 1:nrow(pthway))
  {
    if (igraph::V(graph_routes)$name[i] == pthway$pathway[j])
    {
      igraph::V(graph_routes)$shape[i] = 15
    }
  }
}


ggraph(graph_routes, layout='kk') +
  edge_color+

  geom_node_point(aes_(color=~as.numeric(as.character(color)), size=~size, shape=~shape)) +

  scale_colour_gradient2(name = "fold change", low = "green", high = "red",na.value = "#E5C494")+

  geom_node_text(aes(label = name), repel=TRUE,size= ifelse(V(graph_routes)$name %in% met_path_fc$Metabolite, 4, 3.5),
                 fontface='bold',
                 color = ifelse(V(graph_routes)$name %in% met_path_fc$Metabolite, "blue", "black")) +
  scale_shape_identity()+theme(legend.text = element_text(color = "black", size = 8))

}
