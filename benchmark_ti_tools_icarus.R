#run network similrity metric analysis for all the ti networks with 1) integrated clonal network and 2) expert reference
library(dplyr)

# Load the data
dir = "/Users/tha8tf/MyProjects/ICARUS-files/all_network_files_icarus/"
all_files = list.files(dir)
all_files = all_files[-c(2,3,4)] #removing the groundth truth networks 

expert_net = read.table(paste0(dir, "expert_curated_ref.txt"), sep="\t", header=TRUE, stringsAsFactors = F)
expert_net$directed = FALSE
expert_net = rbind(expert_net, expert_net)
expert_net[1:38, "from"] = expert_net[39:76, "to"]
expert_net[1:38, "to"] = expert_net[39:76, "from"]
expert_net$uid = paste(expert_net$from, expert_net$to, sep="__")

clonal_net = read.table(paste0(dir, "final_integrated_clonal_relaxed_positive.txt"), sep="\t", header=TRUE, stringsAsFactors = F)
clonal_net$uid = paste(clonal_net$from, clonal_net$to, sep="__")

# similarity between expert and clonal network
source("~/CollaborativeProjects/Grimes/Annie_CellTag_ClonalCouplingScores/network_similarity_functions.R")
source("~/MyProjects/dyno_metrics_functions.R")

calculate_jaccard_index(reference.network = expert_net, ti.tool.network = clonal_net,col="uid")
# [1] 0.4429066

icarus_ti(reference.network = expert_net, ti.tool.network = clonal_net, him.threshold = 0.02)

for (i in all_files){
  
  print(i)
  
  ti_net = read.table(paste0(dir, i), sep="\t", header=TRUE, stringsAsFactors = F)
  ti_net$directed = FALSE
  ti_net = rbind(ti_net, ti_net)
  half_mark = nrow(ti_net)/2
  ti_net[1:half_mark, "from"] = ti_net[(half_mark+1):nrow(ti_net), "to"]
  ti_net[1:half_mark, "to"] = ti_net[(half_mark+1):nrow(ti_net), "from"]
  
  ti_net$uid = paste(ti_net$from, ti_net$to, sep="__")
  
  print(paste("jaccard with expert for:", i))
  print(calculate_jaccard_index(reference.network = expert_net, ti.tool.network = ti_net,col="uid"))
  print(paste("jaccard with clonal for:", i))
  print(calculate_jaccard_index(reference.network = clonal_net, ti.tool.network = ti_net,col="uid"))
  
}

colnames(expert_net)[3] = "length"
combined_expert_clonal = rbind(expert_net, clonal_net)


for (i in all_files){
  
  print(i)
  
  ti_net = read.table(paste0(dir, i), sep="\t", header=TRUE, stringsAsFactors = F)
  ti_net$directed = FALSE
  ti_net = rbind(ti_net, ti_net)
  half_mark = nrow(ti_net)/2
  ti_net[1:half_mark, "from"] = ti_net[(half_mark+1):nrow(ti_net), "to"]
  ti_net[1:half_mark, "to"] = ti_net[(half_mark+1):nrow(ti_net), "from"]
  
  ti_net$uid = paste(ti_net$from, ti_net$to, sep="__")
  
  print(paste("jaccard with expert+clonal for:", i))
  print(calculate_jaccard_index(reference.network = combined_expert_clonal, ti.tool.network = ti_net,col="uid"))
  
}

expert_hamming_vector = c()
expert_ipsen_vector = c()

clonal_hamming_vector = c()
clonal_ipsen_vector = c()


for (i in all_files){
  
  print(i)
  
  ti_net = read.table(paste0(dir, i), sep="\t", header=TRUE, stringsAsFactors = F)
  if (any(ti_net$directed)){
    print("directed!")
  } else {
    ti_net$directed = FALSE
    ti_net = rbind(ti_net, ti_net)
    half_mark = nrow(ti_net)/2
    ti_net[1:half_mark, "from"] = ti_net[(half_mark+1):nrow(ti_net), "to"]
    ti_net[1:half_mark, "to"] = ti_net[(half_mark+1):nrow(ti_net), "from"]

    ti_net$uid = paste(ti_net$from, ti_net$to, sep="__")
    
    expert_net$directed = FALSE
    expert_net = rbind(expert_net, expert_net)
    expert_net[1:38, "from"] = expert_net[39:76, "to"]
    expert_net[1:38, "to"] = expert_net[39:76, "from"]
    }

  
  expert_ham_i = icarus_ti(reference.network = expert_net[,1:4], ti.tool.network = ti_net[,1:4], him.threshold = NULL)
  expert_ipsen_i = icarus_ti_ipsen(reference.network = expert_net[,1:4], ti.tool.network = ti_net[,1:4], him.threshold = NULL)
  
  expert_hamming_vector = c(expert_hamming_vector, expert_ham_i)
  expert_ipsen_vector = c(expert_ipsen_vector, expert_ipsen_i)
  
  expert_net = expert_net_og
}




