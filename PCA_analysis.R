
######### PCA(Must be validated) ###########

# Apply PCA with scaling (standardizes variables)
pca_result <- prcomp(df_colinearity, center = TRUE, scale. = TRUE)

# Summary of explained variance
summary(pca_result)

# Scree plot: proportion of variance explained
plot(pca_result, type = "l", main = "Scree Plot (PCA)")

# Biplot: samples and variables together
biplot(pca_result,
       xlim = c(-0.2, 0.3),   
       ylim = c(-0.2, 0.6),
       cex = 0.5,
       expand = 1.2,
       cex.lab = 0.7,         
       cex.axis = 0.7,        
       col = c("grey", "blue"))  # color: primer per a mostres, segon per a variables

# Eigenvalues and proportion of variance
explained_var <- pca_result$sdev^2
prop_var <- explained_var / sum(explained_var)
cumulative_var <- cumsum(prop_var)

# Display proportion of variance explained
print(data.frame(
  PC = paste0("PC", 1:length(prop_var)),
  Variance = round(prop_var, 3),
  Cumulative = round(cumulative_var, 3)
))

# Contributions to each PC
loadings <- pca_result$rotation # Principal loadings
sorted_loadings <- lapply(1:ncol(loadings), function(i) {
  sort(abs(loadings[, i]), decreasing = TRUE)}) # Put in order names and contributions
names(sorted_loadings) <- colnames(loadings) # Asign names to each PC
sorted_loadings # Show contributions in order for each PC 
