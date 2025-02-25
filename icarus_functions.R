
clean_excel <- function(input_excel,col_numbers){
  input_excel = as.data.frame(input_excel)
  input_excel = input_excel[,col_numbers]
  colnames(input_excel) = input_excel[1,]
  input_excel = input_excel[-1,]
  input_excel = na.omit(input_excel)

  return(input_excel)
}

icarus_ti <- function(reference.network, ti.tool.network,him.threshold=0.02){
  
  # perform checks on the provided reference network(s), the TI tool network, and multilineage scores
  # expects to see "from" nodes in first col,"to" nodes in second col, and an edge score in the third col of .network inputs and rest are edge attributes
  
  # check if cell types in "from" col and "to" col match 
  # check if the third col is numeric -- should be called weights 
  # check if the nodes among the provided networks/scores match
  colnames(reference.network) = c("from","to","length","directed")
  reference.network$length = as.numeric(reference.network$length)
  colnames(ti.tool.network) = c("from","to","length","directed")
  ti.tool.network$length = as.numeric(ti.tool.network$length)
  him_score = calculate_him_score(reference.network,ti.tool.network, ga.threshold = him.threshold)
  
  #him_score = paste("HIM Score:",abs(him_score))
  
  return (him_score)
}

icarus_ti_ham <- function(reference.network, ti.tool.network,him.threshold=0.02){
  
  # perform checks on the provided reference network(s), the TI tool network, and multilineage scores
  # expects to see "from" nodes in first col,"to" nodes in second col, and an edge score in the third col of .network inputs and rest are edge attributes
  
  # check if cell types in "from" col and "to" col match 
  # check if the third col is numeric -- should be called weights 
  # check if the nodes among the provided networks/scores match
  colnames(reference.network) = c("from","to","length","directed")
  reference.network$length = as.numeric(reference.network$length)
  colnames(ti.tool.network) = c("from","to","length","directed")
  ti.tool.network$length = as.numeric(ti.tool.network$length)
  him_score = calculate_him_score_hamming(reference.network,ti.tool.network, ga.threshold = him.threshold)
  
  #him_score = paste("HIM Score:",abs(him_score))
  
  return (him_score)
}

icarus_ti_ipsen <- function(reference.network, ti.tool.network,him.threshold=0.02){
  
  # perform checks on the provided reference network(s), the TI tool network, and multilineage scores
  # expects to see "from" nodes in first col,"to" nodes in second col, and an edge score in the third col of .network inputs and rest are edge attributes
  
  # check if cell types in "from" col and "to" col match 
  # check if the third col is numeric -- should be called weights 
  # check if the nodes among the provided networks/scores match
  colnames(reference.network) = c("from","to","length","directed")
  reference.network$length = as.numeric(reference.network$length)
  colnames(ti.tool.network) = c("from","to","length","directed")
  ti.tool.network$length = as.numeric(ti.tool.network$length)
  him_score = calculate_him_score_ipsen(reference.network,ti.tool.network, ga.threshold = him.threshold)
  
  #him_score = paste("HIM Score:",abs(him_score))
  
  return (him_score)
}

icarus_cell_transition <- function(reference.node.weights, multilineage.scores){
  
  # perform checks on the provided reference network(s), the TI tool network, and multilineage scores
  # expects to see "from" nodes in first col,"to" nodes in second col, and an edge score in the third col of .network inputs and rest are edge attributes
  
  # check if cell types in "from" col and "to" col match 
  # check if the third col is numeric -- should be called weights 
  # check if the nodes among the provided networks/scores match
  
  return (score_list)
}

get_sample_observed_matrix <- function(turn,clones,clone.cell.groups,sampling=0.25,seed=10){
  
  set.seed(seed)
  
  sample_clones = sample(rownames(clones),sampling*nrow(clones))
  sample_clones = clones[sample_clones,]
  
  observed_matrix = as.data.frame(matrix(NA,nrow = length(unique(clone.cell.groups)), ncol =length(unique(clone.cell.groups))))
  rownames(observed_matrix) = unique(clone.cell.groups)
  colnames(observed_matrix) = unique(clone.cell.groups)
  
  for (i in rownames(observed_matrix)){
    for (j in colnames(observed_matrix)){
      shared_clones = length(intersect(sample_clones$clone.id[sample_clones$groups==i], sample_clones$clone.id[sample_clones$groups==j]))
      observed_matrix[i,j] = shared_clones
    }
  }
  
  row_total = rowSums(observed_matrix)
  col_total = colSums(observed_matrix)
  total = sum(col_total)
  
  expected_matrix = as.data.frame(matrix(NA,nrow = length(unique(clone.cell.groups)), ncol =length(unique(clone.cell.groups))))
  rownames(expected_matrix) = unique(clone.cell.groups)
  colnames(expected_matrix) = unique(clone.cell.groups)
  
  for (i in rownames(expected_matrix)){
    for (j in colnames(expected_matrix)){
      e = (row_total[i]*col_total[j])/total
      expected_matrix[i,j] = e
    }
  } 
  
  eij_oij = observed_matrix/expected_matrix
  rm(observed_matrix,expected_matrix,sample_clones,row_total,col_total,total,shared_clones)
  
  return(eij_oij)
  
}

get_clonal_network <- function(cells,clone.id,cell.groups,n=500,sampling=0.25,random.clones=TRUE){
  
  clones = as.data.frame(clone.id, row.names = cells)
  ## add a check here to make sure clone.id is character and not numeric vals 
  clones$groups = cell.groups
  
  ## if random.clones = FALSE else do the following -- add an if statement 
  if (random.clones==FALSE){ 
    eij_oij_list= purrr::map(cli::cli_progress_along(1:5),get_sample_observed_matrix,clones=clones,clone.cell.groups = cell.groups,sampling=1)
    } else{
    eij_oij_list= purrr::map(cli::cli_progress_along(1:n),get_sample_observed_matrix,clones=clones,clone.cell.groups = cell.groups,sampling=sampling)
  }
  
  ## 3-D array of n layers of groups by groups shared clones matrices
  a <- array(unlist(eij_oij_list), c(dim(eij_oij_list[[1]]), length(eij_oij_list))) 
  # na_matrices = which(is.na(a[1,1,]))  ########## CAN YOU HAVE NA VALUES HERE? ##########
  # if (length(na_matrices > 0 )){a = a[,,-na_matrices]} else {a = a}
  median_vals = apply(a, 1:2, median, na.rm=T) ## compute the median values of the groups by groups shared clones matrices across n layers
  
  rownames(median_vals) = unique(cell.groups)
  colnames(median_vals) = unique(cell.groups)

  clonal_network_edges = reshape2::melt(median_vals)
  colnames(clonal_network_edges) = c("from","to","clonal_coupling_score")
  
  return (clonal_network_edges)
}


visualize_clonal_network <- function(clonal.network.edges,edge.weights="clonal_coupling_score",prune=FALSE,prune.threshold=0.6,direction="undirected",weight=TRUE,return.network=FALSE){
  
  median_vals = reshape2::dcast(clonal.network.edges,from ~ to,value.var = edge.weights) ## take melted network edge data frame to make it into a matrix form
  rownames(median_vals) = median_vals[,1]
  median_vals = as.matrix(median_vals[,-1])
  diag(median_vals) = 0
  
  g_median_vals = igraph::graph_from_adjacency_matrix(median_vals, mode=direction,weighted = weight)
  
  set.seed(14)
  
  if (prune == TRUE){
    g_median_vals <- igraph::delete_edges(g_median_vals, which(igraph::E(g_median_vals)$weight < prune.threshold))
    isolated = which(igraph::degree(g_median_vals)==0)
    g_median_vals = igraph::delete.vertices(g_median_vals, isolated)
    
    g.layout<-igraph::layout_with_fr(g_median_vals) ## why not layout_with_fr???
    
    network_plot <- plot(g_median_vals, layout = g.layout, edge.width=igraph::E(g_median_vals)$weight, vertex.size = 10, main="Pruned Network from Adjacency Matrix")
  } else {
    
    g.layout <-igraph::layout_with_fr(g_median_vals) ## fore-directed layout
    
    network_plot <- plot(g_median_vals, layout = g.layout, edge.width=igraph::E(g_median_vals)$weight, vertex.size = 10, main="Constructed from Adjacency Matrix")
  }
  
  if (return.network == TRUE){
    g_median_vals_edges = cbind(igraph::get.edgelist(g_median_vals), igraph::E(g_median_vals)$weight)
    colnames(g_median_vals_edges) = c("from","to","weight")
    return(g_median_vals_edges)} else{
      return(network_plot)
    }
}

visualize_clonal_heatmap <- function(clonal.network.edges){
  
  heatmap_plot <- ggplot2::ggplot(data=clonal.network.edges, ggplot2::aes(x=from, y=to, fill=clonal_coupling_score)) +
    ggplot2::geom_tile() + 
    ggplot2::scale_fill_gradientn(colours=c("white","red","dark red"), limits=c(0,10)) + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1,size=12),aspect.ratio = 1/1,axis.text.y=ggplot2::element_text(size=12))
  
  return(heatmap_plot)
}

calculate_node_weights <- function(network.edges,direction="undirected",numeric_translation){
  
  ## check if the "From" and "to" columns are numeric or character, if numeric then move on, otherwise do the following
  
  colnames(network.edges) = c("from","to","length","directed")
  
  network.edges$from = as.character(network.edges$from)
  network.edges$to = as.character(network.edges$to)
  
  for (i in numeric_translation[,1]){
    if(sum(network.edges$from==i) == 0){
      network.edges[network.edges$to==i,2] = numeric_translation[numeric_translation[,1] == i,2]
    } else if (sum(network.edges$to==i) == 0){
      network.edges[network.edges$from==i,1] = numeric_translation[numeric_translation[,1] == i,2]
    } else {
      network.edges[network.edges$from==i,1] = numeric_translation[numeric_translation[,1] == i,2]
      network.edges[network.edges$to==i,2] = numeric_translation[numeric_translation[,1] == i,2]
    }
  }
  
  
  network.edges = network.edges[,1:3]
  network.edges$from = as.numeric(network.edges$from)
  network.edges$to = as.numeric(network.edges$to)
  network.edges$length = as.numeric(network.edges$length)
  
  network.edges.node.weights = tnet::degree_w(network.edges, alpha=0,type = "out")
  network.edges.node.weights = as.data.frame(network.edges.node.weights)
  num_to_node = numeric_translation
  rownames(num_to_node) = num_to_node[,2]
  network.edges.node.weights$node_name = num_to_node[as.character(network.edges.node.weights$node),1]
  
  if (direction=="undirected"){
    network.edges.node.weights_in = tnet::degree_w(network.edges, alpha=0,type = "in")
    network.edges.node.weights_in = as.data.frame(network.edges.node.weights_in)
    network.edges.node.weights$node_name = num_to_node[as.character(network.edges.node.weights$node),1]
    network.edges.node.weights[,2:3] = network.edges.node.weights[,2:3] + network.edges.node.weights_in[,2:3]
  }
  
  ## assumption is that "from" and "to" columns are numeric values by this point 
  rownames(network.edges.node.weights) = network.edges.node.weights$node_name
  network.edges.node.weights = network.edges.node.weights[,-4]
  return(network.edges.node.weights)
  
}



net_similarity_metrics <- function(net1,net2,method="euclidean",direction=TRUE,weighted=TRUE){
  
  if (method == "euclidean"){
    d_net = sqrt((net1 - net2)^2)
  } else if (method == "manhattan"){
    d_net = sqrt((net1 - net2)^2)
  } else if (method == "jaccard"){
    
  }
  
  
}

calculate_him_score <- function(
    net1,
    net2,
    simplify = TRUE, 
    ga.threshold=0.02
) {
  
  # get the matched adjacencies
  adjacencies <- get_matched_adjacencies(
    net1,
    net2,
    simplify = simplify
  )
  
  # return 0 when the largest length of either graph is 0
  if (max(adjacencies[[2]]) == 0 || max(adjacencies[[1]]) == 0) {
    return(0)
  }
  
  netdist <- nettools::netdist(
    adjacencies[[1]] / sum(adjacencies[[1]]),
    adjacencies[[2]] / sum(adjacencies[[2]]),
    "HIM",
    n.cores = 1,
    ga = ga.threshold
  )["HIM"]
  
  netdist[netdist < 0] <- 0
  
  1 - netdist
}

calculate_him_score_hamming <- function(
    net1,
    net2,
    simplify = TRUE, 
    ga.threshold=0.02
) {
  
  # get the matched adjacencies
  adjacencies <- get_matched_adjacencies(
    net1,
    net2,
    simplify = simplify
  )
  
  # return 0 when the largest length of either graph is 0
  if (max(adjacencies[[2]]) == 0 || max(adjacencies[[1]]) == 0) {
    return(0)
  }
  
  netdist <- nettools::netdist(
    adjacencies[[1]] / sum(adjacencies[[1]]),
    adjacencies[[2]] / sum(adjacencies[[2]]),
    "H",
    n.cores = 1,
    ga = ga.threshold
  )["H"]
  
  netdist[netdist < 0] <- 0
  
  1 - netdist
}


calculate_him_score_ipsen <- function(
    net1,
    net2,
    simplify = TRUE, 
    ga.threshold=0.02
) {
  
  # get the matched adjacencies
  adjacencies <- get_matched_adjacencies(
    net1,
    net2,
    simplify = simplify
  )
  
  # return 0 when the largest length of either graph is 0
  if (max(adjacencies[[2]]) == 0 || max(adjacencies[[1]]) == 0) {
    return(0)
  }
  
  netdist <- nettools::netdist(
    adjacencies[[1]] / sum(adjacencies[[1]]),
    adjacencies[[2]] / sum(adjacencies[[2]]),
    "IM",
    n.cores = 1,
    ga = ga.threshold
  )["IM"]
  
  netdist[netdist < 0] <- 0
  
  1 - netdist
}
