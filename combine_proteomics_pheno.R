install.packages("R.methodsS3", lib = "/n/home_fasse/txia/R/libs", repos = "https://cloud.r-project.org")


# Install R.oo if not already installed
if (!requireNamespace("R.oo", quietly = TRUE)) {
  install.packages("R.oo", lib = "/n/home_fasse/txia/R/libs", repos = "https://cloud.r-project.org")
}

# Install R.utils if not already installed
if (!requireNamespace("R.utils", quietly = TRUE)) {
  install.packages("R.utils", lib = "/n/home_fasse/txia/R/libs", repos = "https://cloud.r-project.org")
}

# Load the necessary libraries from the specified library location
library(R.methodsS3, lib.loc = "/n/home_fasse/txia/R/libs")
library(data.table, lib.loc = "/n/home_fasse/txia/R/libs")
library(R.oo, lib.loc = "/n/home_fasse/txia/R/libs")
library(R.utils, lib.loc = "/n/home_fasse/txia/R/libs")
# Continue with your existing code


# Load the necessary library from the specified library location
library(data.table, lib.loc = "/n/home_fasse/txia/R/libs")

# Define the file paths for the proteomics data
path <- "/n/holylfs05/LABS/liang_lab_l3/Projects/biobank_genomics_l3/UKBiobank/phen45052/Olink_proteomics/"

# Read in the main proteomics dataset
Olink_NPX_Baseline <- fread(paste0(path, "Olink_NPX_Basline_for_analysis_11.10.2023.txt.gz"))
# Rename 'eid' to 'f.eid'
setnames(Olink_NPX_Baseline, old = "eid", new = "f.eid")

# Optionally read in batch and LODs data if needed for preprocessing
Olink_Batch_Baseline <- fread(paste0(path, "Olink_Batch_Basline_for_analysis_11.10.2023.txt.gz"))
Olink_LODs_Baseline <- fread(paste0(path, "Olink_LODs_Basline_for_analysis_11.10.2023.txt.gz"))
setnames(Olink_Batch_Baseline, old = "eid", new = "f.eid")
setnames(Olink_LODs_Baseline, old = "eid", new = "f.eid")

# Assuming 'f.eid' is the common identifier
# Merge batch and LODs data with NPX data if necessary
proteomics_combined <- merge(Olink_NPX_Baseline, Olink_Batch_Baseline, by = "f.eid", all.x = TRUE)
proteomics_combined <- merge(proteomics_combined, Olink_LODs_Baseline, by = "f.eid", all.x = TRUE)

# Load the phenotype data (assumed to be already loaded as 'prot_pheno_diet')
prot_pheno_diet <- fread("/n/home_fasse/txia/UKB_CKM/Data/CKMdefine_export.csv")

# Merge the combined proteomics data with phenotype data
combined_dataset <- merge(proteomics_combined, prot_pheno_diet, by = "f.eid", all.x = TRUE)

# Save the combined dataset
fwrite(combined_dataset, "/n/home_fasse/txia/UKB_CKM/Data/combined_proteomics_pheno.csv")

# Display structure and summary for verification
str(combined_dataset)
summary(combined_dataset)



