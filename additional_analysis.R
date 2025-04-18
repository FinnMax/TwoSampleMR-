install_github("phenoscanner/phenoscanner")

res_phenoscanner <- phenoscanner(snpquery="rs...")
head(res$results)
res$snps

head(data_harmonized_1)


if (!require("devtools")) { install.packages("devtools") } else {}
  devtools::install_github("rondolab/MR-PRESSO")
#Do MRPRESSO analysis
results_PRESSO_1 <- mr_presso(BetaOutcome = "beta.outcome",
                     BetaExposure = "beta.exposure",
                     SdOutcome = "se.outcome",
                     SdExposure = "se.exposure",
                     data = data_harmonized_1,
                     OUTLIERtest = TRUE,  # Set to TRUE if you want to test for outliers
                     DISTORTIONtest = TRUE,  # Set to TRUE if you want to test for distortion
                     SignifThreshold = 0.05,  # Significance threshold
                     NbDistribution = 10000,  # Number of distributions to simulate
                     seed = 123)  # Set seed for reproducibility (optional)
print(results_PRESSO_1)

# Define the RSIDs of the outliers
outliers <- c("rs....")

data_harmonized_cleaned <- data_harmonized_1 %>%
  filter(!SNP %in% outliers)

# Save the dataset
write.table(data_harmonized_cleaned, 
            file = "harmonized_datatxt", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

report_2 <- mr_report(
  data_harmonized_cleaned,
  output_path = "After PRESSO",
  output_type = "html",
  path = system.file("reports", package = "TwoSampleMR")
)

steiger_data <- steiger_filtering(data_harmonized_1)
