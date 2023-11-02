library('MetaSKAT')
library(SKAT)


# Parameters
args = unlist(strsplit(commandArgs(trailingOnly = TRUE)," "))
inputPath <- args[1]
outputPath <- args[2]
list_cohort <- args[3:length(args)]


# Work path
setwd(inputPath)


# Data path
n_cohort <- length(list_cohort)
File.Mat.vec <- rep("", n_cohort)
File.Info.vec <- rep("", n_cohort)
Group_Idx <- c()


# Load cohort data
for (nc_ in 1:length(list_cohort)){
  
  cohort_ <- list_cohort[nc_]
  print(cohort_)
  
  # Data path
  File.Mat <- sprintf("./%s.cohort_%s.MSSD", cohort_, cohort_)
  File.SetInfo <- sprintf("./%s.cohort_%s.MInfo", cohort_, cohort_)	

  # List
  File.Mat.vec[nc_] <- File.Mat
  File.Info.vec[nc_] <- File.SetInfo
  
  # Group_Idx
  if (cohort_ == 'ZIB_AD'){
    Group_Idx <- c(Group_Idx, 1)
  } else{
    Group_Idx <- c(Group_Idx, 2)
  }
  
}


# Execute the SKAT-O analysis
Cohort.Info <- Open_MSSD_File_2Read(File.Mat.vec, File.Info.vec)
print('Start analysis')
table_res <- MetaSKAT_MSSD_ALL(Cohort.Info, combined.weight=FALSE, is.separate = TRUE, MAF.cutoff=0.05, method='optimal', Group_Idx=Group_Idx)

# Write
write.table(table_res, outputPath, quote=F, row.names=F, sep=",")


