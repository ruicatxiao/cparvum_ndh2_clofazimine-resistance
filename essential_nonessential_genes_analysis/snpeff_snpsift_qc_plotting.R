# R plotting for sra2vcf and snpeff output

library(chromoMap)
library(tidyverse)
library(forcats)
ibrary(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(ggpubr)
library(ggpmisc) 


chromoMap("cp_bgf_chromoMap_genome.txt","4chromomap_11_gene.txt",
          data_based_color_map = T,
          data_type = "categorical",
          canvas_height = 5000,
          canvas_width = 1000,
          data_colors = list(c("yellow","green"))
          )


chromoMap("cp_bgf_chromoMap_genome.txt","4chromomap_sample_cov_per_high_impact.txt",
          data_based_color_map = T,
          data_type = "numeric",
          plots = "scatter",
          heat_map = F,
          canvas_height = 3000,
          canvas_width = 3000,
          plot_height = 100,
          anno_col = c("black"),
          ref_line = T,
          refl_pos = 75,
          plot_ticks = 4,
          plot_y_domain = list(c(0,100)),
          plot_filter = list(c("lt",20,"grey"))
          )



d1 <- read.table("snpeff_variant_effects_impact.txt", header=TRUE)

d1long <- d1 %>%
  pivot_longer(
    cols = c(Moderate, Low, Modifier, High),
    names_to = "Category",
    values_to = "Count"
  )

sample_order <- d1long %>%
  group_by(Sample, Species) %>%
  summarise(Total = sum(Count, na.rm = TRUE)) %>%
  arrange(desc(Total)) %>%
  ungroup()

d1long <- d1long %>%
  left_join(sample_order, by = c("Sample", "Species")) %>%
  mutate(Sample = fct_reorder(Sample, Total, .desc = TRUE))

ggplot(d1long, aes(x = Sample, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  facet_wrap(~ Species, scales = "free_y") +
  labs(title = "Hominis_vs_Parvum_Variant_Impact_Counts",
       x = "Sample",
       y = "Count") +
  theme(legend.position = "bottom")


ggplot(d1long, aes(x = Sample, y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = "fill") +
  coord_flip() +
  facet_wrap(~ Species, scales = "free_y") +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Hominis_vs_Parvum_Variant_Impact_Profile",
       x = "Sample",
       y = "Percentage") +
  theme(legend.position = "bottom")

d2 <- read.table("snpeff_variant_effects_region.txt", header=TRUE)

d2long <- d2 %>%
  pivot_longer(
    cols = c(Exon,Downstream,Upstream,Intergenic,UTR5,UTR3,Intron,SpliceSite),
    names_to = "Category",
    values_to = "Count"
  )

sample_order2 <- d2long %>%
  group_by(Sample, Species) %>%
  summarise(Total = sum(Count, na.rm = TRUE)) %>%
  arrange(desc(Total)) %>%
  ungroup()

d2long <- d2long %>%
  left_join(sample_order2, by = c("Sample", "Species")) %>%
  mutate(Sample = fct_reorder(Sample, Total, .desc = TRUE))


ggplot(d2long, aes(x = Sample, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  facet_wrap(~ Species, scales = "free_y") +
  labs(title = "Hominis_vs_Parvum_Variant_Region_Counts",
       x = "Sample",
       y = "Count") +
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position = "bottom")


ggplot(d2long, aes(x = Sample, y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = "fill") +
  coord_flip() +
  facet_wrap(~ Species, scales = "free_y") +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Hominis_vs_Parvum_Variant_Region_Profile",
       x = "Sample",
       y = "Percentage") +
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position = "bottom")




  data <- read_tsv("high_impact_indel.txt", na = c("nan", "NA"))

# -------------------------------
# Data Preprocessing
# -------------------------------

# Define columns of interest for Set 1 and Set 2
af_columns <- c("cpbgf_7001900_AF", "cpbgf_3002080_AF", "cpbgf_4004390_AF")
dp_columns <- c("cpbgf_7001900_DP", "cpbgf_3002080_DP", "cpbgf_4004390_DP")

# For Set 1 (Violin + Boxplot), retain rows without NA in AF columns
violin_data <- data %>%
  select(Species, all_of(af_columns)) %>%
  drop_na()

# For Set 2 (Scatter Plots), retain rows without NA in corresponding DP and AF columns
scatter_data <- data %>%
  select(Species, all_of(dp_columns), all_of(af_columns)) %>%
  drop_na()

# -------------------------------
# Set 1: Violin + Boxplot with Statistical Testing
# -------------------------------

# Reshape data to long format for AF measurements
violin_long <- violin_data %>%
  pivot_longer(
    cols = all_of(af_columns),
    names_to = "Condition",
    values_to = "AF"
  ) %>%
  mutate(
    Condition = factor(Condition, 
                       levels = af_columns,
                       labels = c("7001900_AF", "3002080_AF", "4004390_AF"))
  )

# Define comparisons: list of pairs to compare between Species
comparisons <- list(c("hominis", "parvum"))

# Determine a fixed y-axis position for p-value annotations
fixed_label_y <- max(violin_long$AF, na.rm = TRUE) * 1.05

# Create the violin + boxplot with statistical annotations using Student's t-test
ggplot(violin_long, aes(x = Species, y = AF, fill = Species)) +
  geom_violin(trim = TRUE, position = position_dodge(width = 0.9)) +
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), outlier.shape = NA) +
  facet_wrap(~ Condition, nrow = 1) +
  labs(title = "Violin + Boxplot for AF Measurements by Species",
       x = "Species",
       y = "AF Value") +
  theme(legend.position = "none",
        strip.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_compare_means(
    comparisons = comparisons, 
    method = "t.test", 
    label = "p.format", 
    label.y = fixed_label_y) +
  ylim(0,1)

# -------------------------------
# Set 2: Scatter Plots with Regression
# -------------------------------

# Split scatter_data by Species
scatter_split <- split(scatter_data, scatter_data$Species)

# Function to create individual scatter plots with regression and annotations
create_scatter_plot <- function(df, dp_col, af_col, title) {
  # Fit linear model
  model <- lm(as.formula(paste(af_col, "~", dp_col)), data = df)
  
  # Extract coefficients and R-squared
  intercept <- round(coef(model)[1], 4)
  slope <- round(coef(model)[2], 4)
  r_squared <- round(summary(model)$r.squared, 4)
  
  # Compute Pearson correlation coefficient
  correlation <- round(cor(df[[dp_col]], df[[af_col]], method = "pearson"), 4)
  
  # Create annotation text
  eq_text <- paste0("y = ", slope, "x + ", intercept)
  r2_text <- paste0("R² = ", r_squared)
  r_text <- paste0("r = ", correlation)
  
  # Generate scatter plot
  p <- ggplot(df, aes_string(x = dp_col, y = af_col)) +
    geom_point(color = "steelblue") +
    geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "grey") +  # Dashed grey regression line
    labs(title = title, x = dp_col, y = af_col) +
    theme_minimal() +
    # Add regression equation, R², and correlation coefficient annotations
    annotate("text", 
             x = Inf, y = -Inf, 
             label = paste(eq_text, r2_text, r_text, sep = "\n"),
             hjust = 1.1, vjust = -0.1,
             size = 3) +
    theme(
      plot.title = element_text(size = 10, face = "bold"),
      axis.title = element_text(size = 9)
    )
  
  return(p)
}

# Initialize a list to store all scatter plots
scatter_plots <- list()

# Loop through each species and create scatter plots for each DP-AF pair
for (species in names(scatter_split)) {
  species_data <- scatter_split[[species]]
  
  # Create scatter plots for each DP-AF pair
  plot_a <- create_scatter_plot(species_data, "cpbgf_7001900_DP", "cpbgf_7001900_AF", "7001900_DP vs 7001900_AF")
  plot_b <- create_scatter_plot(species_data, "cpbgf_3002080_DP", "cpbgf_3002080_AF", "3002080_DP vs 3002080_AF")
  plot_c <- create_scatter_plot(species_data, "cpbgf_4004390_DP", "cpbgf_4004390_AF", "4004390_DP vs 4004390_AF")
  
  # Add plots to the list with names indicating the species
  scatter_plots[[paste0(species, "_1")]] <- plot_a
  scatter_plots[[paste0(species, "_2")]] <- plot_b
  scatter_plots[[paste0(species, "_3")]] <- plot_c
}

# Arrange scatter plots: two rows (one per species), three columns
plot2 <- grid.arrange(
  scatter_plots[["hominis_1"]], scatter_plots[["hominis_2"]], scatter_plots[["hominis_3"]],
  scatter_plots[["parvum_1"]], scatter_plots[["parvum_2"]], scatter_plots[["parvum_3"]],
  nrow = 2,
  top = "Scatter Plots with Regression Lines by Species"
)
