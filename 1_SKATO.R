library(MetaSKAT)
library(SKAT)
library(magrittr)
library(data.table)
library(dplyr)


LoadData <- function(inputPath){
  ## Load information of the cohort to test
  
  print(sprintf('Loading data...'))
  
  # Plink and Variant annotation files
  File.Fam <- sprintf("%s/cohort.fam", inputPath)
  File.Cov <- sprintf("%s/cohort.cov", inputPath)
  File.Var <- sprintf("%s/cohort.raw", inputPath)
  File.SetID <- sprintf("%s/cohort.SetID", inputPath)
  
  # Genotype
  df_geno_raw <- fread(File.Var, header=TRUE, stringsAsFactors=FALSE)
  df_geno <- df_geno_raw[,7:(ncol(df_geno_raw)), drop = FALSE] # select only variant columns (from the 7th col to the end)
  colnames(df_geno) <- sapply(colnames(df_geno), function(colname) {
    parts <- strsplit(colname, "_")[[1]]
    return(parts[1])
  })
  
  # Phenotype and Convariates
  df_famcov <- Read_Plink_FAM_Cov(File.Fam, File.Cov, Is.binary=TRUE, flag1=0, cov_header=TRUE)
  
  # Data Check
  if (identical(df_geno_raw$IID, df_famcov$IID)) {
    print("Succ")
  } else {
    print("Error")
  }
  
  # ElementSet
  df_set <- fread(File.SetID, header=FALSE, stringsAsFactors=FALSE)
  
  return(list(df_geno, df_famcov, df_set))
}


RunSKATO <- function(df_famcov, df_geno_element) {
  ## Perform SKAT-O Analysis
  
  if ((ncol(df_geno_element) > 0)) {
    
	# Phenotype
    Cov <- as.matrix(df_famcov[,7:(ncol(df_famcov))]) # Covariances
    y <- df_famcov[,6] # Phenotype
    
    # SKATO
    obj <- SKAT_Null_Model(y ~ Cov, out_type = 'D') # Build obj
    res_skat <- SKAT(as.matrix(df_geno_element), obj, method="optimal.adj", max_maf = 0.05)
  } 

  return(res_skat)
}


RunLogit <- function(df_famcov, df_geno_element) {
  ## Perform Logit regression
  
  if (ncol(df_geno_element) > 0) {
    
    # Data Preprocess
    var_count_ <- rowSums(df_geno_element, na.rm = TRUE)
    df_reg <- df_famcov
    df_reg$VarCount <- var_count_
    df_reg <- df_reg[!is.na(df_reg$Phenotype), ] # remove sample not in AD and NC group
    
    # Logist regression
    if (length(table(df_reg$VarCount)) > 1){
      if ("Edu_yrs" %in% colnames(df_reg)){
        res_logit <- glm(Phenotype ~ VarCount + Gender + Age + Age2 + Edu_yrs + PC1 + PC2 + PC3 + PC4 + PC5 + Num_e4, data = df_reg, family = binomial(link = "logit"))
      } else {
        res_logit <- glm(Phenotype ~ VarCount + Gender + Age + Age2 + PC1 + PC2 + PC3 + PC4 + PC5 + Num_e4, data = df_reg, family = binomial(link = "logit"))
      }
  }
  
  return(res_logit)
}



########## Main ##########
args = unlist(strsplit(commandArgs(trailingOnly = TRUE)," "))

inputPath <- args[1]
outputFile <- args[2]
varclass_ <- args[3]

## For test
# inputPath <- ""
# varclass_ <- 'All'


## Get the cohort information
data.loder <- LoadData(inputPath)
df_geno <- data.loder[[1]]
df_famcov <- data.loder[[2]]
df_set <- data.loder[[3]]


## Get the no-duplication element list
df_set_uniq <- df_set %>% group_by(V1, V2) %>% filter(row_number() == 1) %>% ungroup()
element.list <- unique(df_set_uniq$V1)


## Dataframe to store results
df_res <- data.frame()


## Loop by Element
for (ne_ in 1:length(element.list)){
  element_ <- element.list[ne_]
  print(element_)
  
  # Variant set in the element
  var.list <- unique(df_set_uniq[df_set_uniq['V1'] == element_, ]$V2)

  # Genotype matrix for the variant set
  var.list.exist <- intersect(var.list, colnames(df_geno))
  df_geno_element <- df_geno[, ..var.list.exist]
  
  # SKATO and regression
  df_res_skato <- RunSKATO(df_famcov, df_geno_element) # SKATO
  df_res_logit <- RunLogit(df_famcov, df_geno_element) # Logit
  
}
