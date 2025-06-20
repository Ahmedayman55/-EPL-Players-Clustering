# Load data
data = read.csv("EPL dataset.csv")

# Select numeric columns only
x = data[, c("Age", "MP", "Starts", "Min", "Gls", "Ast", "CrdY", "CrdR")]

# Scale the data
x = scale(x, center = TRUE, scale = TRUE)
rownames(x) <- 1:nrow(x)

# R² Function
R2 = function(x, clusters, k){
  n = nrow(x)
  tss = var(x)
  tss = (n-1) * sum(diag(tss))
  wss = 0
  for(j in 1:k){
    cj = x[clusters == j, , drop=FALSE]
    nj = nrow(cj)
    if (nj > 1) {
      vj = var(cj)
      wssj = (nj-1) * sum(diag(vj))
      wss = wss + wssj
    }
  }
  r2 = 1 - wss / tss
  return(r2)
}

# Distance and linkage methods
distance_methods = c("euclidean", "manhattan", "binary", "canberra", "minkowski")
linkage_methods = c("single", "complete", "average", "centroid", "ward.D2")

# Empty results table
results <- list()

# Main analysis loop
for(dist_method in distance_methods) {
  for(link_method in linkage_methods) {
    # Skip incompatible combinations (Ward methods work best with Euclidean)
    if((link_method %in% c("ward.D", "ward.D2")) && dist_method != "euclidean") next
    
    cat("\n\n=== Analyzing with:", dist_method, "distance and", link_method, "linkage ===\n")
    
    # Calculate distance matrix
    d <- dist(x, method = dist_method)
    
    # Hierarchical clustering
    hc <- tryCatch({
      hclust(d, method = link_method)
    }, error = function(e) {
      cat("Error in hclust for", dist_method, "+", link_method, ":", conditionMessage(e), "\n")
      next
    })
    
    # Cut tree into 3 clusters
    clusters <- cutree(hc, k = 4)
    
    # Plot dendrogram in a new window
    dev.new()
    plot(hc, main = paste("Dendrogram:", dist_method, "distance +", link_method, "linkage"))
    rect.hclust(hc, k = 4, border = "red")
    
    # Calculate cluster centers
    centers <- aggregate(x, list(clusters), mean)
    centers <- centers[, -1]
    k <- nrow(centers)
    
    # Compute R²
    r2 <- R2(x, clusters, 4)
    
    # Store results
    results[[paste(dist_method, link_method, sep = "_")]] <- list(
      distance = dist_method,
      linkage = link_method,
      clusters = clusters,
      centers = centers,
      r_squared = r2,
      cluster_sizes = table(clusters)
    )
    
    # Print summary
    cat("Cluster sizes:\n")
    print(table(clusters))
    cat("R-squared:", r2, "\n")
    cat("Cluster centers:\n")
    print(centers)
  }
}

# After all combinations are run, analyze the results
r2_values <- sapply(results, function(x) x$r_squared)
best_method <- names(which.max(r2_values))

cat("\n\n=== BEST PERFORMING COMBINATION ===")
cat("\nMethod:", best_method)
cat("\nR-squared:", results[[best_method]]$r_squared)
cat("\nCluster sizes:\n")
print(results[[best_method]]$cluster_sizes)

# K-means clustering elbow method
r2_values <- numeric(10)
wss <- numeric(10)

for(k in 1:10){
  kmc <- kmeans(x, centers = k)
  r2_values[k] <- R2(x, kmc$cluster, k)
  wss[k] <- sum(kmc$withinss)
}

# Elbow plot
dev.new()
plot(1:10, wss, type = "b", pch = 19,
     xlab = "k", ylab = "WSS", main = "WSS vs Number of Clusters")

# K-means clustering with k = 4
set.seed(2003)
kmc <- kmeans(x, 4)

# Comparison
hc_best <- results[[best_method]]
hc_r2 <- hc_best$r_squared
km_r2 <- R2(x, kmc$cluster, 3)

cat("\n=== CLUSTERING PERFORMANCE COMPARISON (k=3) ===\n")
cat("Best hierarchical R²:", hc_r2, "(method:", best_method, ")\n")
cat("K-means R²:", km_r2, "\n")

cat("\n=== CLUSTER SIZE COMPARISON ===\n")
cat("Hierarchical cluster sizes:\n")
print(hc_best$cluster_sizes)
cat("\nK-means cluster sizes:\n")
print(table(kmc$cluster))
