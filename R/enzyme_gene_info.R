#'Gets gene and enzyme info
#'@param df_metenrichment metenrichment score dataframe having metabolites, metabolite class and their enrichment score obtained
#'from metclassenrichment
#'@param sps species name
#'@param classm class name of the metabolite whether sub_class, main_class or super_class
#'@export
#'@examples
#'met_gene_info = enzyme_gene_info (metenrichment, "hsa","sub_class")
enzyme_gene_info <- function(df_metenrichment,sps, classm)
{
  res= rxninfo(df_metenrichment)
  res_orthology= pathinfo(res, 'ORTHOLOGY')
  met_orthology_reaction = merge(res[,c("Metabolite","KEGG ID","rxn",classm)],res_orthology, by.x = "rxn",by.y = "ENTRY")
  met_orthology_reaction = plyr::rename(met_orthology_reaction, c('NAME'="Rxn_name","id"="orthology_id"))
  #res_orthology=rename(res_orthology, c("NAME"="Rxn_name","ENTRY"="Rxn_id","id"="orthology_id"))
  query = as.vector(met_orthology_reaction$orthology_id)
  query = unique(query)
  ### pass the argument in list of 10s since keggrest takes 10 inputs
  query_split = split(query,  ceiling(seq_along(query)/10))
  info = llply(query_split, function(x)keggGet(x))
  unlist_info <- unlist(info, recursive = F)
  #print(str(unlist_info[1]))
  # BUG FIX: 11/8/2021: keggGET no more returns DEFINITION
  #extract_info <- lapply(unlist_info, '[', c("ENTRY","NAME","DEFINITION","GENES"))
  extract_info <- lapply(unlist_info, '[', c("ENTRY","NAME","GENES"))
  dd=do.call(rbind, extract_info)
  df = as.data.frame(dd)
  print(colnames(df))
  r=df %>% tidyr::unnest('GENES')
  sps_ind = grep(paste0(sps,":"),r$GENES, ignore.case = TRUE)
  r2=r[sps_ind,]
  r2$GENES=gsub(paste0(casefold(sps,upper = TRUE),": "),"", r2$GENES)
  r2$countgenes=NA
  for (i in 1: nrow(r2))
  {
    r2$countgenes[i] = length(strsplit(r2$GENES[i], " ")[[1]])
  }
  r3=r2 %>%
    uncount(weights = countgenes, .id = "n", .remove = F)
  for (j in 1:nrow(r3)){

    r3[j, 'GENES']= strsplit(r3[,"GENES"][[1]][j], " ")[[1]][r3$n[j]]
  }
  r3$GENES=gsub("\\(.+)","", r3$GENES)

  query = as.vector(r3$GENES)
  query = unique(query)
  query = paste0(sps,":",query)
  ### pass the argument in list of 10s since keggrest takes 10 inputs
  query_split = split(query,  ceiling(seq_along(query)/10))


  ### kegg info
  info = llply(query_split, function(x)keggGet(x))
  unlist_info <- unlist(info, recursive = F)
  # BUG FIX: 11/8/2021: keggGET no more returns DEFINITION
  #extract_info <- lapply(unlist_info, '[', c("ENTRY","NAME","DEFINITION","ORTHOLOGY","ORGANISM","PATHWAY","DBLINKS","MOTIF"))
  extract_info <- lapply(unlist_info, '[', c("ENTRY","NAME","ORTHOLOGY","ORGANISM","PATHWAY","DBLINKS","MOTIF"))

  dd=do.call(rbind, extract_info)
  df = as.data.frame(dd)
  df$ORTHOLOGY = as.character(df$ORTHOLOGY)
  df = plyr::rename(df, c("ENTRY"="gene_id","NAME"="gene_name"))
  df$ORTHOLOGY =as.character(df$ORTHOLOGY)

  orthology_info = lapply(unlist_info, '[', c("ORTHOLOGY"))
  orthology_info2=data.frame()
  for (i in 1:length(unlist_info))
    {
    orthology_info2[i,"ORTHOLOGY"] = unlist_info[[i]][["ORTHOLOGY"]]
    orthology_info2[i,"orthology_id"] = names(unlist_info[[i]][["ORTHOLOGY"]])
  }
  #r2$ENTRY=as.character(r2$ENTRY)
  orthology_info2=unique(orthology_info2)
  orthology_gene_id = merge(df,orthology_info2, by="ORTHOLOGY")
  #names(orthology_gene_id)[1]='orthology_id'
  orthology_rxn_gene = merge(orthology_gene_id , met_orthology_reaction, by = "orthology_id")
  # BUG FIX: 11/8/2021: keggGET no more returns DEFINITION
  #namesc= c("orthology_id" , "ORTHOLOGY.x",   "gene_id", "gene_name", "DEFINITION", "ORGANISM", "PATHWAY", "DBLINKS","MOTIF",       "rxn", "Metabolite",  "KEGG ID", "sub_class", "Rxn_name" ,
  namesc= c("orthology_id" , "ORTHOLOGY.x",   "gene_id", "gene_name",  "ORGANISM", "PATHWAY", "DBLINKS","MOTIF",       "rxn", "Metabolite",  "KEGG ID", "sub_class", "Rxn_name" , "RCLASS","ORTHOLOGY.y","EQUATION","EQUATION_more" ,"ENZYME" )

  ### check if any column name is empty
  if (any(is.na(names(orthology_rxn_gene)))){
    names(orthology_rxn_gene)[which(is.na(names(orthology_rxn_gene)))]=setdiff(namesc, names(orthology_rxn_gene))
  }

  orthology_rxn_gene = dplyr::select(orthology_rxn_gene, -ORTHOLOGY.y)
  orthology_rxn_gene=dplyr::rename(orthology_rxn_gene, ORTHOLOGY=ORTHOLOGY.x)
  #orthology_rxn_gene = rename(orthology_rxn_gene, c("gene_id" ="ENTRY","gene_name"="NAME"))
  #orthology_rxn_gene = orthology_rxn_gene[c("Metabolite","KEGG ID",classm,"Rxn_name","rxn","EQUATION","EQUATION_more","RCLASS", "ENZYME","orthology_id","ORTHOLOGY.x","GENES","NAME","DEFINITION")]
  #orthology_rxn_gene = plyr::rename(orthology_rxn_gene,c("ORTHOLOGY.x"="ORTHOLOGY"))
  return(orthology_rxn_gene)
}

