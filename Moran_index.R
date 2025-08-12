# Authors: Militão, T., González-Solís, J., Bellmunt A.,
# Tittle: Predictive nesting habitat modelling of Pterodroma feae in the Cabo Verde Archipelago
# R code for developing habitat modelling based on Pterodroma feae nesting locations in Cabo Verde archipélago
# Last updated version: 09/06/2025
# For doubts and questions about the code: adria.bellmunt.ribas@gmail.com

# Code explanation: This code shown below performs a pure Moran index based on
# first and robust version deleting and not introducing all the variables that
# not accomplish all Moran's assessments


######### **Moran Index ########

# Radius in meters to compare on Moran’s I
radii <- seq(0, 10000, by = 3000)

# Coordinates directly from projected Longitude and Latitude
coords <- df[, c("Longitude", "Latitude")] |> as.matrix()

# Results container
results <- data.frame()

# Loop through variables and radii
for (var in selected_vars) {
  for (r in radii) {
    
    nb <- dnearneigh(coords, 0, r)
    
    # Identify non-isolated points
    non_isolated <- card(nb) > 0
    
    # Skip if not enough valid points
    if (sum(non_isolated) < 2) {
      cat("Too few non-isolated points for", var, "at radius", r, "\n")
      next
    }
    
    nb_filtered <- subset.nb(nb, subset = non_isolated)
    lw <- nb2listw(nb_filtered, style = "W")
    
    x <- as.numeric(df[[var]])[non_isolated]
    if (any(is.na(x))) next
    
    moran <- moran.test(x, lw)
    
    results <- rbind(results, data.frame(
      Variable = var,
      Radius_m = r,
      Moran_I = moran$estimate[1],
      p_value = moran$p.value
    ))
  }
}
# Plot Moran's I by radius for each variable with p-value
ggplot(results, aes(x = Radius_m, y = Moran_I, group = Variable, color = Variable)) +
  geom_line(linewidth = 1) +
  geom_point() +
  labs(
    title = "Moran's Index across distance (m)",
    x = "Radius (m)",
    y = "Moran's Index"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Explicació de perquè algunes variables no s'acaben representant en l'Índex de Moran:
#         -Variables amb molt poca variació
#         -Variables molt extrems 
#         -Variables amb molta asimetria
#         -Variables amb una distribució anormal

