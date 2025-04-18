install.packages("markdown")
install.packages("Cairo")


file_path_outcome_1 <- 
  "rh_volume_Cerebellum-White-Matter.tsv"
file_path_exposure_1 <- 
  "data_exposure_smoking.txt"

# Read the exposure data into R
exposure_data_1 <- read.table(file_path_exposure_1, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# Read the outcome data into R
outcome_data_1 <- read.table(file_path_outcome_1, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Rename the 'rsids' column to 'SNP'
colnames(outcome_data_1)[colnames(outcome_data_1) == "rsids"] <- "SNP"

# Find common SNPs between exposure and outcome datasets
common_snps_1 <- intersect(exposure_data_1$SNP, outcome_data_1$SNP)

# Filter the exposure dataset to include only common SNPs
exposure_data_filtered_common_1 <- exposure_data_1 %>%
  filter(SNP %in% common_snps_1)

# Filter the outcome dataset to include only common SNPs
filtered_outcome_data_1 <- outcome_data_1 %>%
  filter(SNP %in% common_snps_1)


filtered_outcome_data_1$sample_size <- 33224
outcome_data_merged_filtered <- filtered_outcome_data_1


# Assuming your filtered outcome dataset is named 'filtered_outcome_data'
write.table(outcome_data_merged_filtered, file = "outcome_data.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)

outcome_data_formatted_1 <- read_outcome_data(
  snps = NULL,
  filename="outcome_data.txt",
  sep = "\t",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af",
  pval_col = "pval",
  chr_col = "chrom",
  samplesize_col = "sample_size",
  pos_col = "pos"
)

#Harmonize data
data_harmonized_1 <- harmonise_data(
  exposure_dat = exposure_data_filtered_common_1, 
  outcome_dat = outcome_data_formatted_1
)

write.table(data_harmonized_1, file = "harmonized_data.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)
#Prune data
# data_pruned_1 <- power_prune(data_harmonized_1, method = 2, dist.outcome = "continuous")

#Do straight report
report_1 <- mr_report(
  data_harmonized_1,
  output_path = "Quick report",
  output_type = "html",
  path = system.file("reports", package = "TwoSampleMR")
)

head (data_pruned_1)
#Perform MR
#res <- mr(data_pruned_1, method_list = c("mr_egger_regression", "mr_ivw", "mr_weighted_median"))
#p1 <- mr_scatter_plot(res, data_pruned_1)
#p1[[1]]
print("hello")
