# Load necessary libraries
library(data.table)
library(httr)
library(igraph)
library(RColorBrewer)

# Step 1: Combine Protein IDs Across Multiple Files
file_paths <- c(
  "/Users/ning/Downloads/1vs0_selected_proteins_with_info_bootstrap (6).csv",
  #"/Users/ning/Downloads/2vs0_selected_proteins_with_info_bootstrap (1).csv",
  "/Users/ning/Downloads/2vs1_selected_proteins_with_info_bootstrap (6).csv",
  #"/Users/ning/Downloads/3vs0_selected_proteins_with_info_bootstrap (1).csv",
  #"/Users/ning/Downloads/3vs1_selected_proteins_with_info_bootstrap (1).csv",
  "/Users/ning/Downloads/3vs2_selected_proteins_with_info_bootstrap (6).csv"
  #"/Users/ning/Downloads/4vs0_selected_proteins_with_info_bootstrap (1).csv",
  #"/Users/ning/Downloads/4vs1_selected_proteins_with_info_bootstrap (1).csv",
  #"/Users/ning/Downloads/4vs2_selected_proteins_with_info_bootstrap (1).csv",
  #"/Users/ning/Downloads/4vs3_selected_proteins_with_info_bootstrap (3).csv"
)

# Combine protein IDs and remove ".x" suffix
all_proteins <- unique(unlist(lapply(file_paths, function(file) {
  data <- fread(file)
  gsub("\\.x$", "", data$Protein)  # Standardize protein IDs
})))

# Step 2: Map Protein IDs to STRING IDs
protein_list <- paste(all_proteins, collapse = "%0d")
string_api_url <- "https://string-db.org/api/tsv/get_string_ids"
response <- POST(
  url = string_api_url,
  body = list(
    identifiers = protein_list,
    species = 9606,  # Homo sapiens
    limit = 1
  ),
  encode = "form"
)
string_ids_content <- content(response, as = "text", encoding = "UTF-8")
string_ids <- fread(text = string_ids_content, sep = "\t", header = TRUE)

# Step 3: Retrieve Network Interactions
string_network_url <- "https://string-db.org/api/tsv/network"
response_network <- POST(
  url = string_network_url,
  body = list(
    identifiers = protein_list,
    species = 9606,
    required_score = 700  # High confidence threshold
  ),
  encode = "form"
)
network_content <- content(response_network, as = "text", encoding = "UTF-8")
network_data <- fread(text = network_content, sep = "\t", header = TRUE)
write.csv(network_data, "/Users/ning/Downloads/networkno4vs3/full_network_interactions_new.csv", row.names = FALSE)


# Step 4: Full Interaction Plot and Enhanced Cluster Plot with Larger Distances and Labels on the Right
if (exists("network_data") && nrow(network_data) > 0) {
  g <- graph_from_data_frame(network_data[, c("preferredName_A", "preferredName_B")], directed = FALSE)
  
  # Save the full interaction network plot
  full_network_plot <- "/Users/ning/Downloads/networkno4vs3/full_network_plot_labels_on_right_new.png"
  png(filename = full_network_plot, width = 3000, height = 3000, res = 400)
  
  # Generate layout for the full network with larger distances
  layout_full <- layout_with_fr(
    g,
    niter = 5000,                    # Increase the number of iterations
    area = vcount(g)^3,              # Expand the area for node placement
    repulserad = vcount(g)^3.5       # Increase the repulsion radius
  )
  
  # Scale the layout to further increase distances
  layout_full <- layout_full * 3    # Scale the layout to triple the size
  
  # Plot the full network
  plot(
    g,
    layout = layout_full,
    main = "Full Interaction Network with Labels on Right",
    vertex.size = 3,                 # Adjust node size
    vertex.label.cex = 0.2,          # Smaller text size for labels
    vertex.label.family = "Helvetica",
    vertex.label.color = "black",
    vertex.color = "lightblue",
    vertex.label.dist = 0.5,         # Distance of labels from nodes
    vertex.label.degree = 0,         # Place labels to the right
    edge.width = 0.2,
    edge.color = "gray"
  )
  dev.off()
  print(paste("Full interaction network plot saved to:", full_network_plot))
  
  # Set seed for reproducibility
  set.seed(12345)
  
  # Louvain community detection
  clusters <- cluster_louvain(as.undirected(g))
  membership <- membership(clusters)
  unique_clusters <- unique(membership)
  
  print(paste("Number of clusters detected:", length(unique_clusters)))
  
  # Save cluster info and plots
  output_folder_csv <- "/Users/ning/Downloads/networkno4vs3/Cluster_CSVs_new/"
  output_folder_plots <- "/Users/ning/Downloads/networkno4vs3/Cluster_Plots_Labels_Right_new/"
  output_folder_centrality <- "/Users/ning/Downloads/networkno4vs3/Cluster_Centrality_new/"
  dir.create(output_folder_csv, showWarnings = FALSE)
  dir.create(output_folder_plots, showWarnings = FALSE)
  dir.create(output_folder_centrality, showWarnings = FALSE)
  
  # Loop over clusters
  for (cluster_id in unique_clusters) {
    nodes_in_cluster <- names(membership[membership == cluster_id])
    subgraph <- induced_subgraph(g, vids = nodes_in_cluster)
    
    # Save cluster edges to CSV
    edges_list <- igraph::as_edgelist(subgraph)
    cluster_edges <- data.frame(From = edges_list[,1], To = edges_list[,2])
    
    write.csv(cluster_edges, file = paste0(output_folder_csv, "cluster_", cluster_id, "_edges.csv"), row.names = FALSE)
    
    # Calculate centrality measures
    degree_centrality <- degree(subgraph)
    betweenness_centrality <- betweenness(subgraph)
    closeness_centrality <- closeness(subgraph)
    
    centrality_data <- data.frame(
      Protein = names(V(subgraph)),
      Degree = degree_centrality,
      Betweenness = betweenness_centrality,
      Closeness = closeness_centrality
    )
    
    # Sort by degree centrality (or any preferred metric) and save the data
    centrality_data <- centrality_data[order(-centrality_data$Degree), ]
    write.csv(centrality_data, file = paste0(output_folder_centrality, "cluster_", cluster_id, "_centrality.csv"), row.names = FALSE)
    
    # Print the top proteins
    print(paste("Top proteins for Cluster", cluster_id, ":"))
    print(head(centrality_data, n = 5))  # Display top 5 proteins by degree centrality
    
    # Generate the layout with larger spacing for clusters
    layout_large <- layout_with_fr(
      subgraph,
      niter = 3000,                   # Increase iterations for better layout
      area = vcount(subgraph)^3,      # Increase area for node placement
      repulserad = vcount(subgraph)^3.5 # Increase repulsion radius
    )
    
    # Scale the layout to make it even larger
    layout_large <- layout_large * 3  # Triple the distances
    
    # Save cluster plot as PNG
    cluster_plot <- paste0(output_folder_plots, "cluster_", cluster_id, "_labels_on_right.png")
    png(filename = cluster_plot, width = 2000, height = 2000, res = 400)
    
    plot(
      subgraph,
      layout = layout_large,
      main = paste("Cluster", cluster_id),
      vertex.size = 4,                # Adjust node size
      vertex.label.cex = 0.3,         # Smaller text size for labels
      vertex.label.family = "Helvetica",
      vertex.label.color = "black",
      vertex.color = "lightblue",
      vertex.label.dist = 0.5,        # Distance of labels from nodes
      vertex.label.degree = 0,        # Place labels to the right
      edge.width = 0.5,
      edge.color = "gray"
    )
    dev.off()
    print(paste("Cluster plot with labels on the right saved for Cluster", cluster_id, "to:", cluster_plot))
  }
} else {
  stop("Error: Network data is empty.")
}



# Step 5: Weighted Sum Calculation for Clusters
all_nodes <- names(V(g))
weighted_sums_list <- list()

for (cluster_id in unique_clusters) {
  nodes_in_cluster <- names(membership[membership == cluster_id])
  subgraph <- induced_subgraph(g, vids = nodes_in_cluster)
  weights <- degree(subgraph) / sum(degree(subgraph))  # Normalize weights
  
  # Align weights with all nodes
  weights_full <- setNames(rep(0, length(all_nodes)), all_nodes)
  weights_full[nodes_in_cluster] <- weights
  weighted_sums_list[[paste0("comm", cluster_id, "_weights")]] <- weights_full
}

weighted_sums <- do.call(cbind, weighted_sums_list)
write.csv(weighted_sums, "/Users/ning/Downloads/networkno4vs3/weighted_sums_new.csv", row.names = TRUE)

# Step 6: Enrichment Analysis for Clusters
string_enrichment_url <- "https://string-db.org/api/tsv/enrichment"
output_folder_enrichment <- "/Users/ning/Downloads/networkno4vs3/Cluster_Enrichment_Results_new/"
dir.create(output_folder_enrichment, showWarnings = FALSE)
for (cluster_id in unique_clusters) {
  nodes_in_cluster <- names(membership[membership == cluster_id])
  protein_list_cluster <- paste(nodes_in_cluster, collapse = "%0d")
  
  response_enrichment <- POST(
    url = string_enrichment_url,
    body = list(
      identifiers = protein_list_cluster,
      species = 9606
    ),
    encode = "form"
  )
  
  enrichment_content <- content(response_enrichment, as = "text", encoding = "UTF-8")
  enrichment_data <- fread(text = enrichment_content, sep = "\t", header = TRUE)
  
  if (exists("enrichment_data") && nrow(enrichment_data) > 0) {
    write.csv(enrichment_data, file = paste0(output_folder_enrichment, "cluster_", cluster_id, "_enrichment.csv"), row.names = FALSE)
  } else {
    warning(paste("Error retrieving enrichment data for Cluster", cluster_id))
  }
}

