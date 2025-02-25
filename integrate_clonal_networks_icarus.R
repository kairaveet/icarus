
## integrate the three networks by identifying edges that are seen in at least two networks 

# read invivo network 
invivo_net = read.table("/Users/tha8tf/invivo_cc_scores.txt",sep="\t", header=TRUE, stringsAsFactors = F)

#update cell labels 
cell_type_cols <- c("Var1", "Var2") 

invivo_net[cell_type_cols] <- lapply(invivo_net[cell_type_cols], gsub, pattern = "HSC", replacement = "HSCP-1")
invivo_net[cell_type_cols] <- lapply(invivo_net[cell_type_cols], gsub, pattern = "MPP-MultiLin-2", replacement = "HSCP-2")
invivo_net[cell_type_cols] <- lapply(invivo_net[cell_type_cols], gsub, pattern = "Multi-Lin-2", replacement = "MultiLin-A")
invivo_net[cell_type_cols] <- lapply(invivo_net[cell_type_cols], gsub, pattern = "Multi-Lin-1", replacement = "MultiLin-B")
invivo_net[cell_type_cols] <- lapply(invivo_net[cell_type_cols], gsub, pattern = "Mono_c32", replacement = "Mono-2")
invivo_net[cell_type_cols] <- lapply(invivo_net[cell_type_cols], gsub, pattern = "Mono_c33", replacement = "Mono-3")
invivo_net[cell_type_cols] <- lapply(invivo_net[cell_type_cols], gsub, pattern = "Mono_c3", replacement = "Mono-1")

# read jindal celltag network
celltag_net = read.table("/Users/tha8tf/jindal_lsk_celltag_cc_scores.txt",sep="\t", header=TRUE, stringsAsFactors = F)


# read darlin network
darlin_net = read.table("/Users/tha8tf/darlin_invivo_cc_scores.txt",sep="\t", header=TRUE, stringsAsFactors = F)


# prune each network based on log CC (> 0) and # of shared clones (> 3)
invivo_net = invivo_net[invivo_net$log_cc > 0 | invivo_net$shared_clones > 5,]
celltag_net = celltag_net[celltag_net$log_cc > 0 | celltag_net$shared_clones > 5,]
darlin_net = darlin_net[darlin_net$log_cc > 0 | darlin_net$shared_clones > 5,]

invivo_net$dataset = "invivo"
celltag_net$dataset = "celltag"
darlin_net$dataset = "darlin"

# find the intersections 
invivo_net$edge_id = paste(invivo_net$Var1, invivo_net$Var2, sep = "__")
celltag_net$edge_id = paste(celltag_net$Var1, celltag_net$Var2, sep = "__")
darlin_net$edge_id = paste(darlin_net$Var1, darlin_net$Var2, sep = "__")


# Find intersections between pairs of networks
common_invivo_celltag <- intersect(invivo_net$edge_id, celltag_net$edge_id)
common_invivo_darlin <- intersect(invivo_net$edge_id, darlin_net$edge_id)
common_celltag_darlin <- intersect(celltag_net$edge_id, darlin_net$edge_id)

# Union of edges present in at least two networks
unique_edges <- unique(c(common_invivo_celltag, common_invivo_darlin, common_celltag_darlin))


# compute average log CC score for each edge in the union set

library(dplyr)

# set unique rownames 
rownames(invivo_net) <- invivo_net$edge_id
rownames(celltag_net) <- celltag_net$edge_id
rownames(darlin_net) <- darlin_net$edge_id


# Combine all dataframes into a single one with edge_id as a column
merged_df <- bind_rows(
  invivo_net %>% mutate(edge_id = rownames(.)),
  celltag_net %>% mutate(edge_id = rownames(.)),
  darlin_net %>% mutate(edge_id = rownames(.))
)

# let's remove the self edges first 
# merged_df = merged_df[-which(merged_df$Var1 == merged_df$Var2),]

# Calculate the average score for each edge_id
result <- merged_df %>%
  group_by(edge_id) %>%
  summarise(avg_score = mean(log_cc, na.rm = TRUE), n_shared_clones = sum(shared_clones, na.rm = TRUE))

result$cell_type1 = sapply(strsplit(result$edge_id, "__"), function(x) x[1])
result$cell_type2 = sapply(strsplit(result$edge_id, "__"), function(x) x[2])
result = as.data.frame(result)
rownames(result) <- result$edge_id

result = result[intersect(unique_edges, result$edge_id),]
# result = result[result$avg_score > 0 & result$n_shared_clones > 3,]


library(ggplot2)
pdf("integrated_clonal_log_CC.pdf", width=15,height=8)
p = ggplot(result, aes(cell_type1, cell_type2, fill= avg_score)) + 
  geom_tile() + theme_minimal() +
  scale_fill_gradient2(low="blue", mid = "white", high="red") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
show(p)
dev.off()












