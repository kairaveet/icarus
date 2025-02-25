# analyze darlin dataset from Li et al. 

# read in cell harmony group assignments
darlin_groups = read.table("/Volumes/salomonis2/LabFiles/Kairavee/DARLIN-data-analyses/invivo_dataset_cell_type_projection/cellHarmony/QueryGroups.cellHarmony.txt", sep="\t", stringsAsFactors = F, header = F)

# read in clonal info
clonal_meta <- read.table("/Users/tha8tf/CollaborativeProjects/Grimes/darlin_analysis/DARLIN_clonal_cells_ids_groups.txt", header = T, sep = "\t",stringsAsFactors = F)
colnames(clonal_meta)[3] = "clone.id"

# filter for the cells in the groups data that are clonal cells
rownames(darlin_groups) = darlin_groups$V1
darlin_groups_f = darlin_groups[clonal_meta$cell_id,]
clonal_meta$groups = darlin_groups_f$V2

# visualize the cell types seen in this dataset 

pdf("darlin_clone_cells_cell_type_freqs.pdf", width=15,height=8)
# Assuming your data is in a dataframe called df and the column of interest is 'category_column'
p = ggplot(clonal_meta, aes(x = groups)) +
  geom_bar() +
  theme_minimal() +
  labs(x = "Cell type", y = "Frequency", title = "CLONAL CELLS") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels if needed
show(p)
dev.off()


## perform clonal analysis ## 
source("/Users/tha8tf/CollaborativeProjects/Grimes/shared_clones_functions.R")
source("/Users/tha8tf/CollaborativeProjects/Grimes/clonal_coupling_scores_functions.R")

shared_clones = calculate_modified_shared_clones(clones=clonal_meta,groups=clonal_meta$groups)
cc_scores = calculate_modified_cc(clones=clonal_meta,groups=clonal_meta$groups)

identical(cc_scores$Var1,shared_clones$Var1)
#[1] TRUE
identical(cc_scores$Var2,shared_clones$Var2)
#[1] TRUE
cc_scores$shared_clones = shared_clones$value
colnames(cc_scores)[3] = "cc"

cc_scores[which(cc_scores$cc == "NaN"),"cc"] = NA
cc_scores = na.omit(cc_scores)

cc_scores$log_cc = cc_scores$cc
cc_scores$log_cc[cc_scores$cc >= 30] = 30
cc_scores$log_cc[cc_scores$cc <= 1/30] = 1/30
cc_scores$log_cc = log2(cc_scores$log_cc)

pdf("darlin_log_CC.pdf", width=15,height=8)
p = ggplot(cc_scores, aes(Var1, Var2, fill= log_cc)) + 
  geom_tile() +
  scale_fill_gradient2(low="blue", mid = "white", high="red") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
show(p)
dev.off()

write.table(cc_scores, file = "darlin_invivo_cc_scores.txt", sep = "\t", row.names = F, quote = F,col.names = T)


# save the clonal coupling scores









