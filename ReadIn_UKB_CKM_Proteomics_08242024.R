##################################################
# Load the necessary package from the specified library location
library(data.table, lib.loc = "/n/home_fasse/txia/R/libs")

# Function to select UK Biobank fields
f_select_UKB_field <- function(x, file_output, version = "July_2023") {
  cat("\nCurrent version: ", version, "\n")
  f.number = paste0("f.", x, sep = "\\.", collapse = "|")
  
  if (version == "July_2023") {
    df_head <- paste0("/n/home_fasse/txia/UKB_CKM/Data/", version, "_col.txt")
    if (!file.exists(df_head)) {
      print("header file not exist")
      df <- fread(paste0("/n/holylfs05/LABS/liang_lab_l3/Projects/biobank_genomics_l3/UKBiobank/phen45052/July_2023_full_data/ukb674044.tab"), nrows = 1, header = FALSE)
      df <- data.table(t(df), keep.rownames = TRUE)
      df[, rn := substr(rn, 2, 100)]
      fwrite(df, df_head, sep = "\t", col.names = FALSE)
    }
    cmd <- paste0("egrep   \'", f.number, "\' ", df_head, " > ~/temp.txt")
  }
  print(cmd)
  system(cmd)
  
  temp <- fread("~/temp.txt", header = FALSE, sep = "\t")
  
  if (version == "July_2023") {
    cmd <- paste0("cut -f ", paste(c(1, temp$V1), collapse = ","), " /n/holylfs05/LABS/liang_lab_l3/Projects/biobank_genomics_l3/UKBiobank/phen45052/July_2023_full_data/ukb674044.tab >", file_output)
  }
  print(cmd)
  system(cmd)
  system("rm ~/temp.txt")
}

# Specify input field IDs and output file
x <- c(
  # Beverage exposures
  100160, 100180, 100170, 100190, 26002, 100210, 100220, 100230, 100240,
  100370, 100380, 100390, 100490, 100500, 100530, 100540,
  
  # Covariates
  31, 54, 3166, 100002, 100022, 105010, 20116, 21022, 21000,
  6138, 21822, #education, education measurement time (Touchscreen sign-off timestamp)
  
  21001, 21834, #BMI, BMI measurement time (Biometrics sign-off timestamp)
  48, 21834, #WC, WC measurement time (Biometrics sign-off timestamp)
  4079, 21831, #Diastolic blood pressure, automated reading, BP measurement time (Verbal nterview sign-off timestamp)
  94, 21831,#Diastolic blood pressure, manual reading, BP measurement time (Verbal nterview sign-off timestamp)
  4080, 21831,#Systolic blood pressure, automated reading, BP measurement time (Verbal nterview sign-off timestamp)
  93, 21831,#Systolic blood pressure, manual reading, BP measurement time (Verbal nterview sign-off timestamp)
  30740, 30741, #glucose, Glucose assay date
  30750, 30751, #HbA1c, Glycated haemoglobin (HbA1c) assay date, 
  30760, 30761, #HDL, HDL assay date
  30870, 30871, #TG, TG assay date
  30690, 30691, #cholesteroal, cholesteroal assay date
  30700, 30701, # Serum creatinine, Serum creatinine assay date
  30500, 30502, #urinary albumin, urinary albumin ssay date
  30510, 30512, #URINE creatinine
  20002, 21831, 20008, 20009, #self-reported diagnosis of hypertension, Verbal nterview sign-off timestamp, Interpolated Year when non-cancer illness first diagnosed, Interpolated Age of participant when non-cancer illness first diagnosed
  41270, 41280, #ICD-10, Date of first in-patient diagnosis - ICD10
  41271, 41281, #ICD-9, Date of first in-patient diagnosis - ICD9
  6153, 21822, #Diabetes/blood pressure-lowering medications, Touchscreen sign-off timestamp
  132032, #Date N18 first reported (chronic renal failure)
 "f.eid", 26002, 100200, 100470, 100520, 100550, 100560, 100770, 100800, 100810, 100820, 100830, 100840, 100850, 100860, 102810, 102820, 102830, 102840,
  102850, 102860, 102870, 102880, 102890, 102900, 102910, 102940, 102950, 102960,
  102970, 102980, 100950, 101020,  101090, 101160,
  101230, 101240, 101250, 101260, 101270, 101310, 101350, 101390,
  101430, 101470, 101510, 101550, 
  102710, 102720, 102730, 102740, 102750, 102760, 102770, 102780, 101970, 101980, 101990, 102000,
  102010, 102020, 102030, 102040, 102050, 102060, 102070, 102080, 102090, 102120,
  102130, 102140, 102150, 102170, 102180, 102190, 102200, 102210, 102220, 102230,
  102260, 102270, 102280, 102290, 102300, 102310, 102320, 102330, 102340, 102350, 102360,
  102370, 102380, 102410, 102420, 102430, 102440, 102450, 102460, 102470, 102480, 102490,
  102500, 102530, 102540, 102620, 103010, 103020, 103030,
  103040, 103050, 103060, 103070, 103080, 103090, 103100, 103150,
  103160, 103170, 103180, 103190, 103200, 103210, 103220, 103230, 103250, 103260, 103270,
  103280, 103290, 104000, 104010, 104020,
  104030, 104050, 104060, 104070, 104080, 104090, 104100, 104110, 104120, 104130,
  104140, 104150, 104160, 104170, 104180, 104190, 104200, 104210, 104220, 104230, 104240,
  104250, 104260, 104270, 104280, 104290, 104300, 104310, 104320, 104330, 104340, 104350,
  104360, 104370, 104380, 104410, 104420, 104430, 104440, 104450, 104460, 104470,
  104480, 104490, 104500, 104510, 104520, 104530, 104540, 104550, 104560, 104570, 104580, 
  104590,
  21842, #Date of blood draw: YYYY-MM-DD
  40005, 40006, #Cancer: date of diagnosis; situ carcinoma and benign cancers
  131298, 131300, 131302, 131304, 131306, 131360, 131362, 131366, 131368, 131354, 131364, #date CVD diagnosis
131296, 131058, 131056,   #131298, 131300, 131302, 131304 are Heart attack, 131354 are Heart failure, 131296, 131306 are Chronic ischaemic heart disease, 131360, 131362, 131364, 131366, 131368, 131058 are Stroke, 131056 are Transient ischaemic attack
  131386, # Date I73 first reported (other peripheral vascular diseases)
  131350, # Date I48 first reported (atrial fibrillation and flutter)
  132032, # Date N18 first reported (chronic renal failure)
  74, # #Fasting hours
  22189, #Townsend deprivation index
  130708, #date t2d diagnosis
  131286, 131294, #date hypertension diagnosis
  130814, #date hypercholesterolemia diagnosis
  2724, 2814, #menopausal status & hormone therapy
  20116, 20161, #smoking & pack-years
  22040, #Physical activity
  6155, #multivitamin use
  1568, 1578, 1588, 1598, 1608, #alcohol intake
  100020, #diet on a typical day
  100240, #coffee
  100190, 100200, 100210, 100880, #fruit intake (whole fruit + fruit juices) #fruit added to cereal variable
  26014, 26032,  #MUFA, SFA
  26015, 26016, #n3, n6 pufa
  100890, #dairy
  20088, #nuts/seeds/legumes #hummus variable
  20091, 20092, 20093, 20094, #refined grains
  20003, #Treatment/medication code
42000, 42006, 42008, 42010  ,
6177 ,
2714, 3581, 2734,3829, 2814, 3591, 2824, 2834, 3882 )

# Output file paths
file_output <- "/n/home_fasse/txia/UKB_CKM/Data/phenotype_CKM"
csv_output <- "/n/home_fasse/txia/UKB_CKM/Data/UKB_dataset_08232024.csv"

# Apply function to select fields and save to phenotype_SSB_ASB
f_select_UKB_field(x, file_output, version = "July_2023")

# Read the selected data into R
df <- fread(file_output)

# Write the dataframe to a CSV file
write.table(df, file = csv_output, row.names = TRUE, quote = FALSE, col.names = TRUE, sep = ",")



