# File path specifications needed:
# 1. Path to the directory containing the input file "CFZ.SNP.BSA.filter.table"
# 2. Output directory for all CSV and PDF files
# 3. The script should be run from the directory containing the input file, or specify the full path to the input file

library(magrittr)
library(dplyr)
library("QTLseqr")
library(ggplot2)
library("ggpubr")
library(reshape2)
library(ggridges)
library(doBy)
library(gridExtra)

format_genomic <- function(...) {
  function(x) {
    limits <- c(1e0,   1e3, 1e6)
    #prefix <- c("","Kb","Mb")
    # Vector with array indices according to position in intervals
    i <- findInterval(abs(x), limits)
    # Set prefix to " " for very small values < 1e-24
    i <- ifelse(i==0, which(limits == 1e0), i)
    paste(format(round(x/limits[i], 1),
                 trim=TRUE, scientific=FALSE, ...)
          #  ,prefix[i]
    )
  }
}

# INPUT: Specify path to your input file
BSA <- read.table("CFZ.SNP.BSA.filter.table", sep = '\t', header = TRUE)

# Check if BSA data was loaded successfully
if (!exists("BSA") || nrow(BSA) == 0) {
  stop("Failed to load BSA data. Please check the file path and format.")
}

# Verify expected dimensions
cat("BSA dimensions:", dim(BSA), "\n")

# Extract AD and DP matrices
AD <- BSA[, seq(5, 36, 4)]
DP <- BSA[, seq(6, 36, 4)]

# Function to extract reference depth from AD field
ref.DP <- function(X) {
  if (is.na(X) || X == "") return(NA)
  tryCatch({
    as.numeric(strsplit(as.character(X), ",")[[1]])[1]
  }, error = function(e) NA)
}

# Validate AD dimensions
if (ncol(AD) != 8) {
  stop("Expected 8 columns in AD matrix, got ", ncol(AD))
}

# Calculate reference allele frequencies
refFre.AD <- matrix(ncol = 8, nrow = nrow(AD))
for(j in 1:8) {
  refFre.AD[,j] <- sapply(AD[,j], ref.DP, simplify = "array") / as.numeric(DP[,j])
}

# Set NA for low coverage positions
refFre.AD[as.numeric(DP) < 10] <- NA

# Add back the genomic coordinates and rename columns
colnames(refFre.AD) <- colnames(AD)	
refFre.AD <- cbind(BSA[, 1:4], refFre.AD)
colnames(refFre.AD) <- gsub(".AD", "", colnames(refFre.AD))

# Verify chromosome names exist
unique_chroms <- unique(refFre.AD$CHROM)
cat("Unique chromosomes:", paste(unique_chroms, collapse = ", "), "\n")

# Rename chromosomes to standard format
refFre.AD <- refFre.AD %>%
  mutate(CHROM = case_when(
    CHROM == "CP141118.1" ~ "Chr1",
    CHROM == "CP141119.1" ~ "Chr2",
    CHROM == "CP141120.1" ~ "Chr3",
    CHROM == "CP141121.1" ~ "Chr4",
    CHROM == "CP141122.1" ~ "Chr5",
    CHROM == "CP141123.1" ~ "Chr6",
    CHROM == "CP141124.1" ~ "Chr7",
    CHROM == "CP141125.1" ~ "Chr8",
    TRUE ~ as.character(CHROM)  # Keep other values unchanged
  ))

# OUTPUT: Specify path for output file
write.csv(cbind(refFre.AD, DP[1:nrow(refFre.AD), ]), 
          file = "CFZ.BSA.refFre.AD.csv", 
          row.names = FALSE)

# Filter for extreme allele frequencies (likely to remove outliers)
refFre.AD <- refFre.AD[
  refFre.AD$X04_Seb_Cp_Tom_S4. > 0.9 & 
  refFre.AD$X05_Seb_Cp_Isra_S5. < 0.1, ]

# OUTPUT: Specify path for cleaned output file
write.csv(refFre.AD, 
          file = "CFZ.BSA.clean.refFre.AD.csv", 
          row.names = FALSE)

cat("Dimensions after filtering:", dim(refFre.AD), "\n")

# Function to remove outliers using MAD (Median Absolute Deviation)
outliersMAD <- function(data, MADCutOff = 2.5, replace = NA, values = FALSE, bConstant = 1.4826, digits = 2) {
  # Handle NA values
  clean_data <- data[!is.na(data)]
  if (length(clean_data) < 2) return(data)
  
  absMADAway <- abs((data - median(data, na.rm = TRUE)) / mad(data, constant = bConstant, na.rm = TRUE))
  
  data[absMADAway > MADCutOff] <- replace
  
  if (values == TRUE) { 
    return(round(absMADAway, digits)) #if values == TRUE, return number of mads for each value
  } else {
    return(round(data, digits)) #otherwise, return values with outliers replaced
  }
}

# Function to apply outlier removal in sliding windows
outlierByMAD <- function(x, k) {
  n <- length(x)
  y <- x
  
  for (i in (k + 1):(n - k)) {
    if (!is.na(x[i])) {  # Only process non-NA values
      data <- x[(i - k):(i + k)]
      if (sum(!is.na(data)) > 1) {  # Need at least 2 non-NA values
        y[i] <- outliersMAD(data)[k+1]
      }
    }
  }
  return(y)
}

# Apply outlier filtering to each sample
tryCatch({
  CFZ_1.filter <- outlierByMAD(refFre.AD$CFZ_D_6_10_S7., 50)
  CFZ_2.filter <- outlierByMAD(refFre.AD$CFZ_D_13_16_S8., 50)
  CFZ_3.filter <- outlierByMAD(refFre.AD$CFZ_D_26_28_S9., 50)
  Vehicle_1.filter <- outlierByMAD(refFre.AD$Vehicle_D_6_10_S10., 50)
  Vehicle_2.filter <- outlierByMAD(refFre.AD$Vehicle_D_13_16_S11., 50)
  Vehicle_3.filter <- outlierByMAD(refFre.AD$Vehicle_D_26_28_S12., 50)
}, error = function(e) {
  stop("Error applying outlier filtering: ", e$message)
})

# Combine filtered data
refFre.AD.filter <- cbind(refFre.AD[, 1:4], 
                          CFZ_1.filter,
                          CFZ_2.filter,
                          CFZ_3.filter,
                          Vehicle_1.filter,
                          Vehicle_2.filter,
                          Vehicle_3.filter)

# Tricube smoothing function
tricubeStat <- function(POS, Stat, windowSize = 1e5, ...) {
  if (windowSize <= 0) {
    stop("A positive smoothing window is required")
  }
  
  # Remove NA values for smoothing
  valid_idx <- !is.na(Stat) & !is.na(POS)
  if (sum(valid_idx) < 2) {
    return(rep(NA, length(Stat)))
  }
  
  tryCatch({
    smoothed <- stats::predict(locfit::locfit(Stat[valid_idx] ~ locfit::lp(POS[valid_idx], h = windowSize, deg = 0), ...), POS)
    return(smoothed)
  }, error = function(e) {
    warning("Smoothing failed: ", e$message)
    return(rep(NA, length(Stat)))
  })
}

# Function to apply tricube smoothing to all samples
SAnalysis.1 <- function(SNPset, windowSize = 1e5, ...) {
  SNPset <- SNPset %>%
    dplyr::group_by(CHROM) %>%
    dplyr::mutate(CFZ_1.tricube = tricubeStat(POS = POS, Stat = CFZ_1.filter, windowSize, ...)) %>%
    dplyr::mutate(CFZ_2.tricube = tricubeStat(POS = POS, Stat = CFZ_2.filter, windowSize, ...)) %>%
    dplyr::mutate(CFZ_3.tricube = tricubeStat(POS = POS, Stat = CFZ_3.filter, windowSize, ...)) %>%
    dplyr::mutate(Vehicle_1.tricube = tricubeStat(POS = POS, Stat = Vehicle_1.filter, windowSize, ...)) %>%
    dplyr::mutate(Vehicle_2.tricube = tricubeStat(POS = POS, Stat = Vehicle_2.filter, windowSize, ...)) %>%
    dplyr::mutate(Vehicle_3.tricube = tricubeStat(POS = POS, Stat = Vehicle_3.filter, windowSize, ...)) %>%
    return(as.data.frame(SNPset))	
}

# Apply tricube smoothing
refFre.AD.tri <- SAnalysis.1(refFre.AD.filter)

# Remove first row if needed (based on original code)
if (nrow(refFre.AD.tri) > 1) {
  refFre.AD.tri <- refFre.AD.tri[-1, ]
}

cat("Final dimensions after tricube smoothing:", dim(refFre.AD.tri), "\n")

# OUTPUT: Specify path for final output file
write.csv(refFre.AD.tri, 
          file = "CFZ.clean.BSA.refFre.AD.tri.csv", 
          row.names = FALSE)

# Prepare data for whole genome analysis
refFre.LC <- setNames(melt(refFre.AD.tri[, 5:16]), c('BSAs', 'AlleleFrequency'))
refFre.LC.filter <- refFre.LC[rowSums(is.na(refFre.LC)) == 0, ]

# Calculate summary statistics
sum <- summaryBy(AlleleFrequency ~ BSAs, data = refFre.LC.filter, FUN = list(mean, median))

# OUTPUT: Specify path for summary file
write.csv(sum, 
          file = "Whole.Genome.AF.sum_add.drug.csv", 
          row.names = FALSE)

# Create density plot
p <- ggplot(refFre.LC.filter, aes(x = AlleleFrequency, y = BSAs, fill = "")) +
  geom_density_ridges(scale = 1, alpha = 0.5) + 
  theme_ridges() + 
  theme(axis.text.x = element_blank(), 
        legend.position = c(0.8, 0.8)) +
  scale_y_discrete(expand = c(0.01, 0)) +  
  scale_x_continuous(expand = c(0, 0)) + 
  coord_flip()

# OUTPUT: Specify path for PDF output
pdf('CpxCp_whole_genome_allele_frequency.pdf', width = 8, height = 4)
print(p)
dev.off()

# Create genomic plots
p0 <- ggplot(data = refFre.AD.tri) + 
  scale_x_continuous(
    breaks = seq(from = 0, to = max(refFre.AD.tri$POS, na.rm = TRUE), by = 400000),  # Line every 400,000 base pairs
    labels = function(x) x / 1e6  # Convert to Mb without adding "Mb" to labels
  ) + 
  labs(x = "Genomic Position (Mb)") +  # Set x-axis title
  facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
  theme_classic() + 
  ylim(0, 1) +
  theme(
    strip.background = element_blank(),  # Remove the box around facet labels
    strip.text = element_text(size = 10)  # Customize facet label text
  )

# Plot for CFZ samples
pCpCp_CFZ <- p0 + ylab("BGF allele frequency") +
  geom_line(aes_string(x = "POS", y = "CFZ_1.tricube"), color = "#eba184", size = 0.5) + 
  geom_line(aes_string(x = "POS", y = "CFZ_2.tricube"), color = "#c96958", size = 0.5) +
  geom_line(aes_string(x = "POS", y = "CFZ_3.tricube"), color = "#a72e3f", size = 0.5)

# Plot for Vehicle samples
pCpCp_Vehicle <- p0 + ylab("BGF allele frequency") +
  geom_line(aes_string(x = "POS", y = "Vehicle_1.tricube"), color = "#77b9e8", size = 0.5) + 
  geom_line(aes_string(x = "POS", y = "Vehicle_2.tricube"), color = "#118dd4", size = 0.5) +
  geom_line(aes_string(x = "POS", y = "Vehicle_3.tricube"), color = "#0d6bba", size = 0.5)

# Define highlight region for Chr7
highlight_region <- data.frame(
  CHROM = "Chr7",  # Use the exact label "Chr7"
  xmin = 444049, 
  xmax = 1163985,
  ymin = -Inf,
  ymax = Inf
)

# Combined plot with highlight
pCpCp_Vehicle_CFZ <- p0 + ylab("BGF allele frequency") +
  # Highlight region on chromosome Chr7
  geom_rect(data = highlight_region, 
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
            fill = "#D6C0B3", alpha = 0.5) +  # Highlight with transparency
  geom_point(aes_string(x = "POS", y = "Vehicle_3.filter"), color = "#77b9e8", size = 0.2) + 
  geom_line(aes_string(x = "POS", y = "Vehicle_3.tricube"), color = "#0d6bba", size = 0.5) + 
  geom_point(aes_string(x = "POS", y = "CFZ_3.filter"), color = "#eba184", size = 0.2) + 
  geom_line(aes_string(x = "POS", y = "CFZ_3.tricube"), color = "#a72e3f", size = 0.5) + 
  facet_grid(~ CHROM, scales = "free_x", space = "free_x")

# Create combined plot
allele.frequency_summary <- ggarrange(
  pCpCp_Vehicle, pCpCp_CFZ, pCpCp_Vehicle_CFZ,
  ncol = 1, nrow = 3
)

# OUTPUT: Specify paths for final plots
pdf('pCpCp.allele.frequency_summary.pdf', width = 7, height = 5)
print(ggarrange(
  pCpCp_Vehicle, pCpCp_CFZ, pCpCp_Vehicle_CFZ,
  ncol = 1, nrow = 3
))
dev.off()

# Save as SVG
ggsave(
  filename = "pCpCp.allele.frequency_summary.svg",
  plot = allele.frequency_summary,
  width = 7, height = 5
)

cat("Analysis completed successfully!\n")