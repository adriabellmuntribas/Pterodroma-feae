# Authors: Militão, T., González-Solís, J., Bellmunt A.,
# Tittle: Predictive nesting habitat modelling of Pterodroma feae in the Cabo Verde Archipelago
# R code for developing habitat modelling based on Pterodroma feae nesting locations in Cabo Verde archipélago
# Last updated version: 17/08/2025
# For doubts and questions about the code: adria.bellmunt.ribas@gmail.com

# Code explanation: This code shown below performs a colinearity analysis, and 
# based on this analysis, allows you to remove which correlations higher than 
# thresholds of 0.7-1 are interesting to extract from the model. In addition, this 
# code also extracts graphs in relation to the previous analyses mentioned.strip width of the survey´+0oo


######### Libraries ########### 

if (!requireNamespace("librarian", quietly = TRUE)) {
  install.packages("librarian")
}
suppressPackageStartupMessages(library(librarian))

librarian::shelf(
  dplyr,
  tidyr,
  readr,
  readxl,
  stringr,
  data.table,
  car,
  tibble,
  corrplot,
  spdep,
  tidyverse,
  sf,
  terra,
  Hmisc,
  vapour,
  ggplot2,
  gam,
  mgcv,
  car,
  writexl,
  RColorBrewer,
  pROC)

######### Data loading ###########
df <- ("/Users/adriabellmuntribas/Desktop/TFM/QGIS/Data analysis/Colinearity analysis/finaldatabasefiltered_v7_20250810_2319.csv")
df <- read_csv(df, show_col_types = FALSE)
View(df)
head(df)

######### Preliminary analysis (Pearson) ###########

# Variables to include in the analysis
df <- read.csv("/Users/adriabellmuntribas/Desktop/TFM/QGIS/Data analysis/Colinearity analysis/df_final_random_points.csv")
selected_vars <- c("dem_mean", "northeastness_mean", 
                     "slope_mean", "flowacc_mean_scaled", 
                     "distcoast_mean_scaled", "distcourtship_mean_scaled", 
                     "distdrainage_mean_scaled", "VNL_mean", 
                     "GHI_mean", "WIND_mean", 
                     "NDVI_mean", "NDMI_mean")
df_numeric <- df %>% # Select only specified numeric columns
  select(all_of(selected_vars))

# Compute Pearson correlation matrix
cor_matrix <- cor(df_numeric, method = "pearson")
print(round(cor_matrix, 2))
#write.csv(round(cor_matrix, 2), "cor_matrix.csv")

# Visualize correlation matrix
corrplot(cor_matrix, method = "color", type = "lower",
         tl.col = "black", tl.srt = 45, addCoef.col = "black",
         number.cex = 0.7, diag = FALSE)

# Correlation test with p-values and visualize significant correlations 
cor_test <- rcorr(as.matrix(df_numeric))
p.mat <- cor_test$P
corrplot(cor_matrix,
         p.mat = p.mat,
         method = "circle",
         type = "lower",
         insig = "blank",
         tl.col = "black",
         tl.srt = 45,
         addCoef.col = "black",
         number.cex = 0.75,
         order = "AOE",
         diag = FALSE)
######### Filter high correlation variables (|r| ≥ 0.7) (interactive) ###########

# Identify highly correlated pairs (excluding diagonal and duplicates)
high_corr_pairs <- which(abs(cor_matrix) >= 0.7 & abs(cor_matrix) < 1, arr.ind = TRUE)
high_corr_pairs <- high_corr_pairs[high_corr_pairs[, 1] < high_corr_pairs[, 2], , drop = FALSE]

# If no such pairs, report and skip
if (nrow(high_corr_pairs) == 0) {
  cat("No variable pairs with correlation ≥ 0.7 found.\n")
} else {
  cat("Highly correlated variable pairs (|r| ≥ 0.7):\n")
  excluded_vars <- c()
  
  for (i in 1:nrow(high_corr_pairs)) {
    var1 <- rownames(cor_matrix)[high_corr_pairs[i, 1]]
    var2 <- colnames(cor_matrix)[high_corr_pairs[i, 2]]
    corr_val <- round(cor_matrix[var1, var2], 2)
    
    cat(paste0("\nPair: ", var1, " & ", var2, " | Correlation = ", corr_val, "\n"))
    response <- readline(prompt = paste0("Which variable do you want to EXCLUDE (", var1, "/", var2, ")? Press Enter to skip: "))
    
    if (response %in% c(var1, var2)) {
      excluded_vars <- c(excluded_vars, response)
      cat("Variable", response, "marked for exclusion.\n")
    } else if (response != "") {
      cat("'", response, "' is not a valid option. Skipping this pair.\n")}}
  
  # Exclude selected variables
  if (length(excluded_vars) > 0) {
    df_colinearity <- df_numeric %>% select(-all_of(unique(excluded_vars)))
    cat("\nVariables excluded:", paste(unique(excluded_vars), collapse = ", "), "\n")
  } else {
    cat("\nNo variables were excluded.\n")
    df_colinearity <- df_numeric}}

# Data base with colinearity filter
View(df_colinearity)

# Database without colinearity filter
View(df_numeric)

# Database original 
View(df)

######### **Moran Index ########
source("/Users/adriabellmuntribas/Desktop/TFM/QGIS/Data analysis/Colinearity analysis/moran_warnings.R")

# Radius in meters to compare on Moran’s I
radii <- seq(0, 10000, by = 50)

# Coordinates directly from projected Longitude and Latitude
coords <- df[, c("Longitude", "Latitude")] |> as.matrix()

# Contenidors buits
results <- tibble::tibble(
  Variable = character(), Radius_m = numeric(),
  Moran_I  = numeric(),   p_value  = numeric(),
  method   = character())
issues  <- tibble::tibble(
  Variable = character(), Radius_m = numeric(),
  reason   = character(), detail   = character())

# Bucle
for (var in selected_vars) {
  for (r in radii) {
    out <- safe_moran_test(var, r, df, coords, issues, results)
    issues  <- out$issues
    results <- out$results}}

View(results)
View(issues)

# Plot Moran's I by radius for each variable with p-value
ggplot(results, aes(x = Radius_m, y = Moran_I, group = Variable, color = Variable)) +
  geom_line(linewidth = 0.8) +
  labs(
    title = "Moran's Index across distance (m)",
    x = "Radius (m)",
    y = "Moran's Index"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    text = element_text(family = "Times New Roman", size = 12),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_text(face = "bold"))



# Explicació de perquè algunes variables no s'acaben representant en l'Índex de Moran:
#         -Variables amb molt poca variació
#         -Variables molt extrems 
#         -Variables amb molta asimetria
#         -Variables amb una distribució anormal
#
#         ---ACTUALITZACIÓ 11/08/2025---       
#
#         -Si hi han NAs en algun rang x es salta el rang i continua avaluant
#         -Ara permet que alguns casos caiguin per veïns zero i permet obtenir Moran’s I 
#         -Mínim de punts aillats abans 2 i ara 3: Això, sol, en trauria; però 
#              l’efecte dominant és que entren molts casos que abans petaven per errors/NA.

######### Retrieve categorical and complementary variables ########

# Variables to restore
vars_to_restore <- c("id_temp", "crop_majority", "soil_majority", "srd_majority", "Latitude", "Longitude", "presence","active")

# Variables missing in df_colinearity but present in df
vars_to_add <- vars_to_restore[!(vars_to_restore %in% colnames(df_colinearity)) & 
                                 vars_to_restore %in% colnames(df)]
cat("Available variables that have been included in the current df but were not present:
    \n", vars_to_add, "\n")

# Variables that cannot be added because they do not exist in either df_colinearity or df
missing_vars <- vars_to_restore[!(vars_to_restore %in% colnames(df_colinearity)) & 
                                  !(vars_to_restore %in% colnames(df))]

# Check that the number of rows matches
if (nrow(df) != nrow(df_colinearity)) {
  stop("Row mismatch between df and df_colinearity. Cannot safely merge.")}

# Extract columns to add
cols_to_add <- df[, vars_to_add, drop = FALSE]

# Remove if any of these columns accidentally exist (for safety)
df_final <- df_colinearity[, !(colnames(df_colinearity) %in% vars_to_add), drop = FALSE]

# Order: id_temp first (if needed), the rest at the end
if ("id_temp" %in% vars_to_add) {
  id_temp_col <- cols_to_add[, "id_temp", drop = FALSE]
  other_cols <- cols_to_add[, setdiff(vars_to_add, "id_temp"), drop = FALSE]
  df_final <- cbind(id_temp_col, df_final, other_cols)
} else {df_final <- cbind(df_final, cols_to_add)}

# Identify columns containing any NA
cols_with_na <- colnames(df_final)[colSums(is.na(df_final)) > 0]

# Remove those columns
df_final <- df_final[, !(colnames(df_final) %in% cols_with_na), drop = FALSE]


View(df_final)
write.csv(df_final, file = "/Users/adriabellmuntribas/Desktop/TFM/QGIS/Data analysis/Colinearity analysis/df_final.csv", row.names = FALSE)

######### Analyzing ALL GAM model's possibility ########

df_final <- read.csv("C:/Users/abellmun/Downloads/df_final.csv")
df_final <- df_final[, !(names(df_final) %in% "id_temp")]
#df_final <- subset(df_final, active == 1) #OPTIONAL: TO AVOID INCLUD CONFLICTIVE NESTS
View(df_final)

continuous_vars <- c("dem_mean", "northeastness_mean", 
                     "slope_mean", "flowacc_mean_scaled", 
                     "distcoast_mean_scaled", "distcourtship_mean_scaled", 
                     "distdrainage_mean_scaled", "VNL_mean", 
                     "GHI_mean", "WIND_mean", 
                     "NDVI_mean", "NDMI_mean")

factor_vars <- c("crop_majority", "soil_majority")

# Create all the continuous variable combinations 
all_combinations <- unlist(
  lapply(0:length(continuous_vars), function(k) combn(continuous_vars, k, simplify = FALSE)),
  recursive = FALSE)

model_idx <- 1
results <- list() #create a list for the results

for (cont_combo in all_combinations) {
  
  # Iterate for all categoric variables
  for (j in 0:length(factor_vars)) {
    factor_subsets <- combn(factor_vars, j, simplify = FALSE)
    
    for (fact_combo in factor_subsets) {
      
      if (length(cont_combo) == 0 && length(fact_combo) == 0) next
      
      cont_terms <- if (length(cont_combo) > 0) paste0("s(", cont_combo, ", k=3)", collapse = " + ") else "" #degrees of freedom k=4
      factor_terms <- if (length(fact_combo) > 0) paste(fact_combo, collapse = " + ") else ""
      
      # Construct the formula
      formula_str <- paste("presence ~", paste(c(cont_terms, factor_terms)[c(cont_terms, factor_terms) != ""], collapse = " + "))
      formula <- as.formula(formula_str)
      
      gam_model <- tryCatch(
        gam(formula, family = binomial, data = df_final, select = TRUE),
        error = function(e) NULL)
      
      if (!is.null(gam_model)) {
        results[[model_idx]] <- list(
          vars = c(cont_combo, fact_combo),
          formula = formula,
          aic = AIC(gam_model),
          deviance_explained = summary(gam_model)$dev.expl)
        model_idx <- model_idx + 1}}}}

# Convert the results in a data.frame

results_df <- do.call(rbind, lapply(results, function(x) {
  formula_clean <- gsub("\\s+", " ", paste(deparse(x$formula), collapse = " ")) #delete spaces in formula
  formula_clean <- trimws(formula_clean)
  
  data.frame(
    vars = paste(x$vars, collapse = ", "), #add vars separated in comas
    formula = formula_clean, #add a string showing respectively formula 
    aic = x$aic,
    deviance_explained = x$deviance_explained,
    stringsAsFactors = FALSE)}))


# Order by AIC
results_df <- results_df[order(results_df$aic), ]

# Add ΔAIC
results_df$delta_AIC <- results_df$aic - min(results_df$aic, na.rm = TRUE)

# Add variables number
results_df$num_vars <- sapply(strsplit(results_df$vars, ","), function(x) length(trimws(x)))

# Add AIC weight
w <- exp(-0.5 * results_df$delta_AIC)
results_df$AIC_weight <- w / sum(w, na.rm = TRUE)

# Order just by ΔAIC
results_df <- results_df[order(results_df$delta_AIC), ]

# Order by variable numbers and ΔAIC
results_df <- results_df[order(results_df$num_vars, results_df$delta_AIC), ]


View(results_df)

write.csv(results_df,
          file = "/Users/adriabellmuntribas/Desktop/TFM/QGIS/Data analysis/Colinearity analysis/results_df_sorted4.csv",
          row.names = FALSE)

write_xlsx(results_df, path = "/Users/adriabellmuntribas/Desktop/TFM/QGIS/Data analysis/Colinearity analysis/results_df_sorted4.xlsx")


#---------------------------RELOAD OPTION----------------------------#
df_final <- read.csv("/Users/adriabellmuntribas/Desktop/TFM/QGIS/Data analysis/Colinearity analysis/results_df_sorted5.csv")
View(df_final)
#---------------------------RELOAD OPTION----------------------------#


######### Analyzing ALL GAM model's possibility with 80-20% split data ########

df_final <- read.csv("/Users/adriabellmuntribas/Desktop/TFM/QGIS/Data analysis/Colinearity analysis/df_final_random_points.csv")
df_final <- df_final[, !(names(df_final) %in% "id_temp")]
#df_final <- subset(df_final, active == 1)
summary(df_final)



#YOU MUST DELETE VARIABLES THAT CONTAINS COLINEARITY
continuous_vars <- c("northeastness_mean", 
                     "slope_mean", "flowacc_mean_scaled", 
                     "distcoast_mean_scaled", "distcourtship_mean_scaled", 
                     "distdrainage_mean_scaled", "VNL_mean", 
                     "GHI_mean", "WIND_mean", "NDMI_mean")

factor_vars <- c("crop_majority", "soil_majority")

# Optional section to split data into training (80%) and validation (20%) datasets,
# stratified by presence and pseudo-absence classes.
# This is used to train models on training_set and evaluate on validation_set.

set.seed(123)  # For reproducibility

# Unique id for each row
df_final$row_id <- 1:nrow(df_final)

# Split presence data
pres <- df_final %>% filter(presence == 1)
pres_train <- pres %>% sample_frac(0.8)
pres_valid <- anti_join(pres, pres_train, by = "row_id")

# Split pseudo-absence data
pseudo <- df_final %>% filter(presence == 0)
pseudo_train <- pseudo %>% sample_frac(0.8)
pseudo_valid <- anti_join(pseudo, pseudo_train, by = "row_id")

# Combine into final training and validation sets
training_set <- bind_rows(pres_train, pseudo_train)
validation_set <- bind_rows(pres_valid, pseudo_valid)

# Delete id before adjusting models
training_set <- training_set %>% select(-row_id)
validation_set <- validation_set %>% select(-row_id)


# BEGIN MODEL FITTING WITH GAM ON training_set

# Create all continuous variable combinations 
all_combinations <- unlist(
  lapply(0:length(continuous_vars), function(k) combn(continuous_vars, k, simplify = FALSE)),
  recursive = FALSE)

model_idx <- 1
results <- list() #create a list for the results

for (cont_combo in all_combinations) {
  
  # Iterate for all categorical variable combinations
  for (j in 0:length(factor_vars)) {
    factor_subsets <- combn(factor_vars, j, simplify = FALSE)
    
    for (fact_combo in factor_subsets) {
      
      if (length(cont_combo) == 0 && length(fact_combo) == 0) next
      
      cont_terms  <- if (length(cont_combo) > 0) paste0("s(", cont_combo, ", k=4)", collapse = " + ") else "" # degrees of freedom k=4
      factor_terms <- if (length(fact_combo) > 0) paste0("factor(", fact_combo, ")", collapse = " + ") else ""
      
      # Formula construction
      formula_str <- paste("presence ~",paste(c(cont_terms, factor_terms)[c(cont_terms, factor_terms) != ""], collapse = " + "))
      formula <- as.formula(formula_str)
      
      gam_model <- tryCatch(
        mgcv::gam(
          formula,
          family = binomial(link = "logit"),
          data = training_set,
          method = "REML",
          select = TRUE), # INTRODUCE FALSE IF YOU DON'T WANT PENALITZATION
        error = function(e) NULL)
      
      if (!is.null(gam_model)) {
        # Predictions on validation set to assess model fit
        pred_probs <- tryCatch(
          predict(gam_model, newdata = validation_set, type = "response"),
          error = function(e) rep(NA, nrow(validation_set)))
        
        observed <- validation_set$presence
        epsilon <- 1e-15
        pred_probs <- pmax(pmin(pred_probs, 1 - epsilon), epsilon)  # avoid log(0)
        
        # Compute external deviance as goodness-of-fit measure on validation set
        deviance_external <- -2 * sum(observed * log(pred_probs) + (1 - observed) * log(1 - pred_probs))
        
        results[[model_idx]] <- list(
          vars = c(cont_combo, fact_combo),
          formula = formula,
          aic = AIC(gam_model),
          deviance_explained = summary(gam_model)$dev.expl,
          deviance_external = deviance_external)
        model_idx <- model_idx + 1}}}}


# Convert the results into a data.frame
results_clean <- Filter(Negate(is.null), results)

results_df <- do.call(rbind, lapply(results_clean, function(x) {
  formula_clean <- gsub("\\s+", " ", paste(deparse(x$formula), collapse = " ")) #delete spaces in formula name
  formula_clean <- trimws(formula_clean)
  
  data.frame(
    vars = paste(x$vars, collapse = ", "),
    formula = formula_clean,
    aic = x$aic,
    deviance_explained = x$deviance_explained,
    deviance_external = x$deviance_external,
    stringsAsFactors = FALSE)}))


# Order results by AIC, delta AIC, AIC weights and number of variables
results_df <- results_df[order(results_df$aic), ]
results_df$delta_AIC <- results_df$aic - min(results_df$aic, na.rm = TRUE)
w <- exp(-0.5 * results_df$delta_AIC)
results_df$AIC_weight <- w / sum(w, na.rm = TRUE)
results_df$num_vars <- sapply(strsplit(results_df$vars, ","), function(x) length(trimws(x)))
results_df <- results_df[order(results_df$num_vars, results_df$delta_AIC), ]

# show data.frame
View(results_df)

# Save results and datasets
write.csv(results_df,
          file = "/Users/adriabellmuntribas/Desktop/TFM/QGIS/Data analysis/Colinearity analysis/FINAL RESULTS/results_df1.csv",
          row.names = FALSE)

write.csv(training_set,
          file = "/Users/adriabellmuntribas/Desktop/TFM/QGIS/Data analysis/Colinearity analysis/FINAL RESULTS/training1_set.csv",
          row.names = FALSE)

write.csv(validation_set,
          file = "/Users/adriabellmuntribas/Desktop/TFM/QGIS/Data analysis/Colinearity analysis/FINAL RESULTS/validation1_set.csv",
          row.names = FALSE)

View(results_df)
View(training_set)
View(validation_set)

######### GAM testing and filtering ##########

resultats_df1 <- read.csv("/Users/adriabellmuntribas/Desktop/TFM/QGIS/Data analysis/Colinearity analysis/results_df1.csv")
resultats_df2 <- read.csv("/Users/adriabellmuntribas/Desktop/TFM/QGIS/Data analysis/Colinearity analysis/results_df2.csv")
resultats_df3 <- read.csv("/Users/adriabellmuntribas/Desktop/TFM/QGIS/Data analysis/Colinearity analysis/results_df3.csv")

validation1_set <- read.csv("/Users/adriabellmuntribas/Desktop/TFM/QGIS/Data analysis/Colinearity analysis/validation1_set.csv")
validation2_set <- read.csv("/Users/adriabellmuntribas/Desktop/TFM/QGIS/Data analysis/Colinearity analysis/validation2_set.csv")
validation3_set <- read.csv("/Users/adriabellmuntribas/Desktop/TFM/QGIS/Data analysis/Colinearity analysis/validation3_set.csv")

training1_set <- read.csv("/Users/adriabellmuntribas/Desktop/TFM/QGIS/Data analysis/Colinearity analysis/training1_set.csv")
training2_set <- read.csv("/Users/adriabellmuntribas/Desktop/TFM/QGIS/Data analysis/Colinearity analysis/training2_set.csv")
training3_set <- read.csv("/Users/adriabellmuntribas/Desktop/TFM/QGIS/Data analysis/Colinearity analysis/training3_set.csv")

# Filter best models from each results dataframe with delta_AIC <= 2
best_models_df1 <- subset(resultats_df1, delta_AIC <= 2)
best_models_df2 <- subset(resultats_df2, delta_AIC <= 2)
best_models_df3 <- subset(resultats_df3, delta_AIC <= 2)

# Check the number of best models in each set
cat("Number of best models in df1:", nrow(best_models_df1), "\n")
cat("Number of best models in df2:", nrow(best_models_df2), "\n")
cat("Number of best models in df3:", nrow(best_models_df3), "\n")

# Combine the filtered dataframes into one consolidated dataframe
combined_best_models <- rbind(best_models_df1, best_models_df2, best_models_df3)

# Remcombined_best_models# Remove duplicate formulas (if any) and sort by number of variables and delta_AIC
combined_best_models_unique <- combined_best_models[!duplicated(combined_best_models$formula), ]
combined_best_models_unique <- combined_best_models_unique[order(combined_best_models_unique$num_vars, combined_best_models_unique$delta_AIC), ]
View(combined_best_models_unique)

# Save combined best models for record
write.csv(combined_best_models_unique, "combined_best_models.csv", row.names = FALSE)


######### GAM AUC and TSS analysis ##########

# Function for TSS analysis
calc_tss <- function(obs, pred_probs, threshold = 0.5) {
  pred_bin <- ifelse(pred_probs >= threshold, 1, 0)
  conf_mat <- table(obs, pred_bin)
  if (nrow(conf_mat) < 2 || ncol(conf_mat) < 2) return(NA)
  
  sensitivity <- conf_mat[2, 2] / sum(conf_mat[2, ])
  specificity <- conf_mat[1, 1] / sum(conf_mat[1, ])
  tss <- sensitivity + specificity - 1
  return(tss)}

# Function for processing each result_df with their train/validation set
process_results <- function(results_df, training_set, validation_set) {
  
  # Filtering models with delta_AIC <= 2
  best_models <- subset(results_df, delta_AIC <= 2)
  
  # Create a list to save the metrics
  metrics_list <- list()
  
  for (i in seq_len(nrow(best_models))) {
    formula_str <- best_models$formula[i]
    form <- as.formula(formula_str)
    
    # Adjust model
    gam_model <- tryCatch(
      gam(form, family = binomial, data = training_set, select = TRUE),
      error = function(e) NULL)
    
    if (is.null(gam_model)) next
    
    # Predictions
    pred_probs <- tryCatch(
      predict(gam_model, newdata = validation_set, type = "response"),
      error = function(e) rep(NA, nrow(validation_set)))
    
    # Observed one's
    obs <- validation_set$presence
    
    # Analysis of AUC
    auc_val <- tryCatch(
      auc(obs, pred_probs),
      error = function(e) NA)
    
    # Analysis of TSS
    tss_val <- calc_tss(obs, pred_probs)
    
    # Save results
    metrics_list[[i]] <- data.frame(
      formula = formula_str,
      aic = best_models$aic[i],
      delta_AIC = best_models$delta_AIC[i],
      AUC = auc_val,
      TSS = tss_val,
      num_vars = best_models$num_vars[i])}
  
  # Combine results
  metrics_df <- do.call(rbind, metrics_list)
  return(metrics_df)}

# Process every result_df with their training and validation set
metrics_df1 <- process_results(results_df, training_set, validation_set)
metrics_df2 <- process_results(resultats_df2, training2_set, validation2_set)
metrics_df3 <- process_results(resultats_df3, training3_set, validation3_set)

# Coerce AUC to plain numeric and bind safely (if you don't do it you can't bind rows)
dfs <- Filter(Negate(is.null), list(metrics_df1, metrics_df2, metrics_df3))
dfs <- Filter(Negate(is.null), list(metrics_df1))


# Coerce AUC (and TSS, just in case) to numeric in each df
dfs <- lapply(dfs, function(x) {
  if ("AUC" %in% names(x))  x$AUC  <- as.numeric(x$AUC)
  if ("TSS" %in% names(x))  x$TSS  <- as.numeric(x$TSS)
  x})

final_metrics <- dfs %>%
  dplyr::bind_rows() %>%
  dplyr::distinct(formula, .keep_all = TRUE) %>%
  dplyr::left_join(
    results_df %>% dplyr::select(formula, vars),
    by = "formula") %>%
  dplyr::arrange(num_vars, delta_AIC)


View(final_metrics)

write_xlsx(final_metrics,
           path = "/Users/adriabellmuntribas/Downloads/final_metrics.xlsx")

######### Analyzing individually GAM model's (interactive and VIF analysis) ##########

df_final <- read.csv("/Users/adriabellmuntribas/Desktop/TFM/QGIS/Data analysis/Colinearity analysis/df_final_random_points.csv") #CAUTION you must include the whole dataset (NOT TRAINING one's)
df_final <- df_final[, !(names(df_final) %in% "id_temp")]

continuous_vars <- c("dem_mean", "northeastness_mean", 
                     "slope_mean", "flowacc_mean_scaled", 
                     "distcoast_mean_scaled", "distcourtship_mean_scaled", 
                     "distdrainage_mean_scaled", "VNL_mean", 
                     "GHI_mean", "WIND_mean", 
                     "NDVI_mean", "NDMI_mean")

factor_vars <- c("crop_majority", "soil_majority") 

#REMINDER: before this step, check firstly the PEARSON CORRELATION:
#----------0.7 or > 0.8, depending on the threshold you want!!

repeat {
  vars <- names(df_final)
  cat("Available variables:\n", paste(vars, collapse = ", "), "\n")
  input_vars <- readline("Enter the variables you want to include (comma-separated): ")
  if (!nzchar(input_vars)) {
    cat("No variables entered. Please try again.\n\n")
    next}
  selected_vars <- trimws(unlist(strsplit(input_vars, ",")))
  
  # Variable validation
  missing_vars <- setdiff(selected_vars, vars)
  if (length(missing_vars) > 0) {
    cat("Error: The following variables don't exist in the dataset:\n")
    cat(paste0(" - ", missing_vars), sep = "\n")
    cat("Please try again.\n\n")} else {
      cat("All the variables are valid.\n")
      break}}

# Build formula dynamically with smoothing terms
formula_lm <- as.formula(paste("presence ~", paste(selected_vars, collapse = " + ")))
df_final$presence_num <- as.numeric(df_final$presence) #VIF don't needs a binominal response, just numerical
lm_model <- lm(update(formula_lm, presence_num ~ .), data = df_final)
vif_values <- vif(lm_model)
print(vif_values) #Variables with VIF > 5 or > 10 are usually problematic.


# Convert to factor BEFORE model CROP or SOIL
df_final$crop_majority <- as.factor(df_final$crop_majority)
# Adjust GAM Model
formula_gam <- presence ~ s(slope_mean, k = 4) + s(flowacc_mean_scaled, k = 4) + 
  s(distcoast_mean_scaled, k = 9) + s(distcourtship_mean_scaled, k = 4) + 
  s(distdrainage_mean_scaled, k = 10) + s(GHI_mean, k = 12) + s(WIND_mean, k = 12) + 
  s(NDMI_mean, k = 4) + crop_majority

gam_full <- mgcv::gam(
  formula_gam,
  family = binomial,
  data = df_final,
  method = "REML") # Select = TRUE forbids analyze just not significative results

summary(gam_full) #The significance of each predictor and the degree 
# of smoothing (edf) can be observed: predictors with high p-value (>0.05) 
# and very low edf (<0.5) may not contribute anything and be candidates 
# for elimination.


# Plot smooth terms
cat("\n===== MODEL SUMMARY =====\n")
plot(gam_full, pages = 1, residuals = TRUE, shade = TRUE, scale = 0, family = "Times New Roman")
title("GAM MODEL", outer = TRUE, line = -1,family = "Times New Roman")
print(summary(gam_full))

# Suavitzation curves
par(mfrow = c(ceiling(length(selected_vars)/2), 2))  # automatic dispose 
plot(gam_full, pages = 1, residuals = TRUE, shade = TRUE, seWithMean = TRUE)
mtext(paste("GAM model with:", paste(selected_vars, collapse = ", ")), outer = TRUE, line = -2, cex = 1.2)




######### Habitat mapping Tiff creator ##########

if (!requireNamespace("librarian", quietly = TRUE)) {
  install.packages("librarian")
}
suppressPackageStartupMessages(library(librarian))

librarian::shelf(
  terra, #terra shows some problems with other packages (shares some commands)
  mgcv)


output <- setwd("/Users/adriabellmuntribas/Desktop/TFM/QGIS/Gis project/26Nreprojected/fogo/final_layers")
list.files(output)

pred1 <- rast("slope_fogo_26N_bilinear.tif")                         # slope_mean
pred2 <- rast("flowaccumulation_26N_fogo_aligned.tif")              # flowacc_mean_scaled
pred3 <- rast("distance_fogo_26N_aligned.tif")                      # distcoast_mean_scaled
pred4 <- rast("distancecourtship_26N_fogo_aligned.tif")             # distcourtship_mean_scaled
pred5 <- rast("distance_drainage_fogo_26N_aligned.tif")             # distdrainage_mean_scaled
pred6 <- rast("GHI_fogo_26N_bilinear_aligned.tif")                  # GHI_mean
pred7 <- rast("wind10_fogo_26N_bilinear_aligned.tif")               # WIND_mean
pred8 <- rast("NDMI_26N_fogo_aligned.tif")                          # NDMI_mean
pred9 <- rast("carta_agricola_26N_fogo_aligned.tif")                # crop_majority

preds <- c(pred1, pred2, pred3, pred4, pred5, pred6, pred7, pred8, pred9) # Apile them

# Rename to be same names with GAM models
names(preds) <- c("slope_mean", 
                  "flowacc_mean_scaled", 
                  "distcoast_mean_scaled", 
                  "distcourtship_mean_scaled", 
                  "distdrainage_mean_scaled", 
                  "GHI_mean", 
                  "WIND_mean", 
                  "NDMI_mean", 
                  "crop_majority")


# Ensure crop_majority is treated as factor with correct levels
preds[["crop_majority"]] <- as.factor(preds[["crop_majority"]])
levels(preds[["crop_majority"]]) <- data.frame(
  value = 1:6,
  label = as.character(1:6))

preds <- scale_rasters_like_model(preds, df_final)

# Prediction of adjusted GAM
map_pred <- terra::predict(preds, gam_full, type="link")

map_pred <- 1 / (1 + exp(-map_pred))
# Save as a TIFF
writeRaster(map_pred, "gam_prediction1.tif", overwrite=TRUE)
rmdir(output)
