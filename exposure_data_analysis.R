#Establish the local variables
smoking_exposure_filepath <- "GSCAN_SmkInit_2022_GWAS_SUMMARY_STATS_EUR.txt"
jwt_token = "token"
bfile_directory = "/EUR"
plink = "plinkbinr/bin/plink_Windows.exe" 
.libPaths()
.libPaths("/renv/library")
install.packages("LDlinkR")

Sys.getenv("HOME")
install.packages("remotes")
install.packages("devtools") 
install.packages("TwoSampleMR", repos = c("https://mrcieu.r-universe.dev", "https://cloud.r-project.org"))

library(remotes)
install_github("MRCIEU/TwoSampleMR")

devtools::install_github("explodecomputer/genetics.binaRies")
genetics.binaRies::get_plink_binary()

if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("MRCIEU/genetics.binaRies")

#Now to read the exposure data set
smoking_exposure_data_file <- read.table(smoking_exposure_filepath, 
                   header = TRUE, 
                   sep = "\t",  
                   stringsAsFactors = FALSE)

#Transform dataframe into correct format  
smoking_exposure_data <- format_data(
  smoking_exposure_data_file,
  type = "exposure",
  snp_col  = "RSID",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "AF_1000G",
  effect_allele_col = "EFFECT_ALLELE",
  other_allele_col = "OTHER_ALLELE",
  pval_col = "P",
  samplesize_col = "N",
  chr_col = "CHR",
  pos_col = "POS",
) 
head(smoking_exposure_data)

filtered_data <- smoking_exposure_data[smoking_exposure_data$pval.exposure <= 5E-08, ]

filtered_data_2 <- filtered_data %>%
  filter(eaf.exposure >= 0.01)
# Some snps removed due to duplications
head(filtered_data)

# Filter SNPs based on MAF >= 1%
smoking_exposure_data <- smoking_exposure_data %>%
  filter(eaf.exposure >= 0.01)


# Perform LD clumping
smoking_exposure_data_clumped  <-  ld_clump(
 dplyr::tibble(rsid=smoking_exposure_data$SNP, 
               pval=smoking_exposure_data$pval.exposure, 
               id=smoking_exposure_data$id.exposure),
 clump_kb = 10000,
 clump_r2 = 0.001,
 clump_p = 5e-8,
 pop = "EUR",
 opengwas_jwt = jwt_token,
 bfile = bfile_directory,
 plink_bin = plink
)

# Now merge back all the information about rsid to the clunped data set
# Rename columns in the first dataset

smoking_exposure_data_clumped <- smoking_exposure_data_clumped %>%
  rename(
    SNP = rsid,                 # Rename 'rsid' to 'SNP'
    pval.exposure = pval,       # Rename 'pval' to 'pval.exposure'
    id.exposure = id            # Rename 'id' to 'id.exposure'
  )

# Perform the inner join with the second dataset
smoking_exposure_data_combined <- smoking_exposure_data_clumped %>%
  inner_join(smoking_exposure_data, by = "SNP")

# Rename columns in the dataset
smoking_exposure_data_combined <- smoking_exposure_data_combined %>%
  rename(
    pval.exposure = pval.exposure.x,  # Rename pval.exposure.x to pval.exposure
    id.exposure = id.exposure.x        # Rename id.exposure.x to id.exposure
  )


# Remove columns named 'id.exposure.y' and 'pval.exposure.y'
smoking_exposure_data_combined <- smoking_exposure_data_combined[, !(names(smoking_exposure_data_combined) %in% c("id.exposure.y", "pval.exposure.y"))]

# Calculate the F-statistic for each SNP
smoking_exposure_data_combined_F <- smoking_exposure_data_combined %>%
  mutate(F_statistic = (beta.exposure^2) / (se.exposure^2))

smoking_exposure_data_combined_F <- smoking_exposure_data_combined_F %>%
  dplyr::filter(F_statistic >= 10)

# Function to check if a SNP is palindromic
is_palindromic <- function(allele1, allele2) {
  (allele1 == "A" & allele2 == "T") |
    (allele1 == "T" & allele2 == "A") |
    (allele1 == "C" & allele2 == "G") |
    (allele1 == "G" & allele2 == "C")
}

# Apply the function to filter out palindromic SNPs
filtered_data_exposure_final <- smoking_exposure_data_combined_F %>%
  dplyr::filter(!is_palindromic(effect_allele.exposure, other_allele.exposure))

# Remove column named 'F_statistic'
filtered_data_exposure_final <- filtered_data_exposure_final[, !(names(filtered_data_exposure_final) %in% "F_statistic")]

# Save the filtered exposure data as a space-separated text file
write.table(filtered_data_exposure_final, file = "data_exposure_smoking.txt", sep = "\t", row.names = FALSE, quote = FALSE)


