install.packages("inflection")
# Paquets
library(dplyr)
library(tidyverse)
library(ggplot2)
library(readr)
library(inflection)     # Per detectar punt d'inflexió
library(car)            # Per test d'homogeneïtat (Levene)
library(stats)          # ANOVA i Kruskal-Wallis


# Ruta i fitxers
ruta <- "/Users/adriabellmuntribas/Desktop/TFM/QGIS/Gis project/26Nreprojected/santoantao/produced_layers/stream_NDVI/taules_NDVI"
arxius <- c(
  "4_5_santoantao_26N_cor_est.csv",
  "5_6_santoantao_26N_cor_est.csv",
  "6_7_santoantao_26N_cor_est.csv",
  "7_8_santoantao_26N_cor_est.csv",
  "8_9_santoantao_26N_cor_est.csv",
  "9_10_santoantao_26N_cor_est.csv"
)

# Carrega i unifica dades, filtrant area = 0
dades <- lapply(arxius, function(fitxer) {
  read_csv(file.path(ruta, fitxer)) %>%
    filter(area > 0) %>%  # 💥 FILTRE APLICAT AQUÍ
    mutate(threshold_range = gsub("_santoantao_26N_cor_est.csv", "", fitxer))
}) %>% bind_rows()


view(dades)

# Assigna valor numèric de threshold mitjà (per regressió i inflexió)
dades <- dades %>%
  mutate(threshold_numeric = as.numeric(sub("_.*", "", threshold_range)))

# 1. Regressió lineal NDVI mitjà ~ threshold
model <- lm(`_NDVImean` ~ threshold_numeric, data = dades)
summary(model)

# Ajust del model lineal
mod <- lm(`_NDVImean` ~ threshold_numeric, data = dades)
pred_df <- data.frame(
  threshold_numeric = unique(dades$threshold_numeric)
)
pred_df$fit <- predict(mod, newdata = pred_df)

# Línia de regressió amb valors etiquetats en posició superior dreta
ggplot(dades, aes(x = threshold_numeric, y = `_NDVImean`)) +
  geom_point(color = "grey20", size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", color = "grey40", se = TRUE, linewidth = 1.2) +
  geom_text(data = pred_df,
            aes(x = threshold_numeric, y = fit, label = round(fit, 3)),
            hjust = -0.2, vjust = -0.8,
            size = 3.5, family = "Times New Roman", color = "grey20") +
  scale_x_continuous(limits = c(4, 10.2), breaks = 4:9.5) +  # Aquí afegim marge dret
  labs(
    x = "Threshold value",
    y = "Mean NDVI"
  ) +
  theme_minimal(base_family = "Times New Roman") +
  theme(
    text = element_text(family = "Times New Roman"),
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14, face = "italic"),
    axis.title.x = element_text(size = 14, margin = margin(t = 12)),
    axis.title.y = element_text(size = 14, margin = margin(r = 12)),
    axis.text = element_text(size = 12),
    plot.caption = element_text(size = 10, hjust = 1),
    panel.grid.major = element_line(color = "grey80")
  )

# 2. Detecció de punt d’inflexió
# Cal mitjanitzar per threshold
ndvi_mitjanes <- dades %>%
  group_by(threshold_numeric) %>%
  summarise(mean_ndvi = mean(`_NDVImean`, na.rm = TRUE))

# Mètode de Kneedle (inflection package)
inf <- findiplist(ndvi_mitjanes$threshold_numeric, ndvi_mitjanes$mean_ndvi, 0)
inflexio <- inf[2]  # Segon punt d'inflexió (últim); recordem que detecta la posició
# Converteix a valor real de threshold
inflexio_real <- ndvi_mitjanes$threshold_numeric[inflexio]

cat("Punt d'inflexió detectat a threshold aproximat:", inflexio_real, "\n")

ggplot(ndvi_mitjanes, aes(x = threshold_numeric, y = mean_ndvi)) +
  geom_line(color = "grey50", size = 1) +                     # Línia grisa
  geom_point(size = 2) +
  scale_x_continuous(breaks = 4:10, limits = c(4, 10)) +     # Thresholds de 5 a 10
  labs(
    x = "Threshold",
    y = "NDVI mitjà"
  ) +
  theme_minimal(base_family = "Times New Roman") +           # Font Times
  theme(
    axis.title.x = element_text(hjust = 0.5, size = 12, color = "black", family = "Times New Roman"),
    axis.title.y = element_text(hjust = 0.5, size = 12, color = "black", family = "Times New Roman"),
    axis.text = element_text(color = "black", family = "Times New Roman"),
    plot.title = element_text(hjust = 0.5, family = "Times New Roman", size = 14)
  )

# 3. Test d’hipòtesi
# Prova d’homogeneïtat de variància
leveneTest(`_NDVImean` ~ threshold_range, data = dades)

# Si NO hi ha homogeneïtat -> Kruskal-Wallis
kruskal.test(`_NDVImean` ~ threshold_range, data = dades)

# Si SÍ hi ha homogeneïtat -> ANOVA
anova_model <- aov(`_NDVImean` ~ threshold_range, data = dades)
summary(anova_model)

