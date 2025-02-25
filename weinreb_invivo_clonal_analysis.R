
# reanalyze the invivo dataset for clonal analysis

# read in cell harmony group assignments
groups = read.table("/Volumes/salomonis2/LabFiles/Kairavee/icarus_final_evals/main_files/final_groups_invivo.txt", header = F, sep = "\t", stringsAsFactors = F)

# read in clonal info
clones_meta = read.table("/Users/tha8tf/Desktop/LabWork/MyProjects/inVivo_newLabels_clones/cellbarcode_to_Clone.txt", sep="\t", header=T, stringsAsFactors=F)
colnames(clones_meta) = c("cell.bc", "clone.id")
clones_meta = clones_meta[clones_meta$clone.id != "others",]

# filter for the cells in the groups data that are clonal cells
rownames(groups) = groups$V1
clonal_groups = groups[groups$V1 %in% clones_meta$cell.bc,]
rownames(clones_meta) = clones_meta$cell.bc
clones_meta = clones_meta[clonal_groups$V1,]

identical(rownames(clones_meta), rownames(clonal_groups))
clones_meta$groups = clonal_groups$V2

# visualize the cell types seen in this dataset 

pdf("weinreb_invivo_clone_cells_cell_type_freqs.pdf", width=15,height=8)
# Assuming your data is in a dataframe called df and the column of interest is 'category_column'
p = ggplot(clones_meta, aes(x = groups)) +
  geom_bar() +
  theme_minimal() +
  labs(x = "Cell type", y = "Frequency", title = "CLONAL CELLS") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels if needed
show(p)
dev.off()


## perform clonal analysis ## 
source("/Users/tha8tf/CollaborativeProjects/Grimes/shared_clones_functions.R")
source("/Users/tha8tf/CollaborativeProjects/Grimes/clonal_coupling_scores_functions.R")

shared_clones = calculate_modified_shared_clones(clones=clones_meta,groups=clones_meta$groups)
cc_scores = calculate_modified_cc(clones=clones_meta,groups=clones_meta$groups)

identical(cc_scores$Var1,shared_clones$Var1)
#[1] TRUE
identical(cc_scores$Var2,shared_clones$Var2)
#[1] TRUE
cc_scores$shared_clones = shared_clones$value
colnames(cc_scores)[3] = "cc"

cc_scores$log_cc = cc_scores$cc
cc_scores$log_cc[cc_scores$cc >= 30] = 30
cc_scores$log_cc[cc_scores$cc <= 1/30] = 1/30
cc_scores$log_cc = log2(cc_scores$log_cc)

ggplot(cc_scores, aes(Var1, Var2, fill= log_cc)) + 
  geom_tile() +
  scale_fill_gradient2(low="blue", mid = "white", high="red") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels if needed


write.table(cc_scores, file = "invivo_cc_scores.txt", sep = "\t", row.names = F, quote = F,col.names = T)
