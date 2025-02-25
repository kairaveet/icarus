get_adjacency_lengths <- function(net, nodes = sort(unique(c(net$from, net$to)))) {
  if(nrow(net) == 0) { # special case for circular
    newnet <- matrix(rep(0, length(nodes)))
    dimnames(newnet) <- list(nodes, nodes)
  } else {
    newnet <- net %>%
      mutate(from = factor(from, levels = nodes), to = factor(to, levels = nodes)) %>%
      reshape2::acast(from~to, value.var = "length", fill = 0, drop = FALSE, fun.aggregate = sum)
  }
  
  (newnet + t(newnet))
  #newnet
}

get_adjacency <- function(net) {
  get_adjacency_lengths(net) != 0
}

# add extra rows and columns to matrix
complete_matrix <- function(mat, dim, fill = 0) {
  mat <- rbind(
    mat,
    matrix(
      rep(fill, ncol(mat) * (dim - nrow(mat))),
      ncol = ncol(mat),
      dimnames = list(sample.int(dim-nrow(mat)), colnames(mat))
    )
  )
  mat <- cbind(
    mat,
    matrix(
      rep(fill, nrow(mat) * (dim - ncol(mat))),
      nrow = nrow(mat),
      dimnames = list(rownames(mat), sample.int(dim-nrow(mat)))
    )
  )
  
  mat
}

# get the matched adjacency matrices between two networks
get_matched_adjacencies <- function(net1, net2, simplify = TRUE) {
  if (simplify) {
    directed1 <- any(net1$directed)
    directed2 <- any(net2$directed)
    net1 <- net1 %>%
      rename(weight = length) %>%
      filter(!(from == to))%>% # remove self loop edges with length 0 & weight == 0)
      igraph::graph_from_data_frame(directed = F) %>%
      igraph::as_data_frame() %>%
      rename(length = weight) %>%
      mutate(directed = directed1) 
    net2 <- net2 %>%
      rename(weight = length) %>%
      filter(!(from == to)) %>% # remove self loop edges with length 0 & weight == 0)
      igraph::graph_from_data_frame(directed = F) %>%
      igraph::as_data_frame() %>%
      rename(length = weight) %>%
      mutate(directed = directed2)
  }
  
  adj1 <- get_adjacency_lengths(net1)
  adj2 <- get_adjacency_lengths(net2)
  
  # make the adjacency matrices have the same dimensions
  if (nrow(adj1) > nrow(adj2)) {
    adj2 <- complete_matrix(adj2, nrow(adj1), fill = 0)
  } else {
    adj1 <- complete_matrix(adj1, nrow(adj2), fill = 0)
  }
  
  lst(adj1, adj2)
}


# get the matched adjacency matrices between two networks
get_matched_adjacencies <- function(net1, net2, simplify = TRUE) {
 all_nodes = unique(c(net1$from,net1$to, net2$from, net2$to))
 
 # make a symmetrical matrix full of zeros with rows and cols as all nodes
 adj1 = matrix(0, nrow = length(all_nodes), ncol = length(all_nodes))
 rownames(adj1) = all_nodes
 colnames(adj1) = all_nodes
 
 # fill the matrix with the lengths of the edges
 for (i in rownames(adj1)) {
   for (j in colnames(adj1)) {
       adj1[i, j] = net1[net1$from == i & net1$to == j,"length"]
   }
  }
}


# get the matched adjacency matrices between two networks
get_matched_adjacencies <- function(net1, net2, simplify = TRUE) {
  all_nodes = unique(c(net1$from,net1$to, net2$from, net2$to))
  
  # make a symmetrical matrix full of zeros with rows and cols as all nodes
  adj1 = matrix(0, nrow = length(all_nodes), ncol = length(all_nodes))
  rownames(adj1) = all_nodes
  colnames(adj1) = all_nodes
  
  # fill the matrix with the lengths of the edges
  for (i in rownames(adj1)) {
    for (j in colnames(adj1)) {
      # if net1$from == i & net1$to == j does not exist in net1, then adj1[i, j] = 0
      if (any(net1$from == i & net1$to == j)){
        adj1[i, j] = net1[net1$from == i & net1$to == j,"length"]
      }
    }
  }
  
  
  # make a symmetrical matrix full of zeros with rows and cols as all nodes
  adj2 = matrix(0, nrow = length(all_nodes), ncol = length(all_nodes))
  rownames(adj2) = all_nodes
  colnames(adj2) = all_nodes
  
  # fill the matrix with the lengths of the edges
  for (i in rownames(adj2)) {
    for (j in colnames(adj2)) {
      # if net2$from == i & net2$to == j does not exist in net2, then adj2[i, j] = 0
      if (any(net2$from == i & net2$to == j)){
        adj2[i, j] = net2[net2$from == i & net2$to == j,"length"]
      }
    }
  }
  
  
  return(lst(adj1, adj2))
}
