#INPUT ARGUMENTS
args <- commandArgs(trailingOnly = T)

set1_directory <- args[1]
set2_directory <- args[2]
set1_mutatedSampleSize <- args[3]
set2_mutatedSampleSize <- args[4]

print("--------------------------------------------------")
print("Step 5: Cross Train Preprocessing")
print("--------------------------------------------------")

#READ CSVS
print("Reading feature and labels.")
feature_and_labels1 <- read.csv(paste0(set1_directory, "/feature_and_labels.csv"))
feature_and_labels2 <- read.csv(paste0(set2_directory, "/feature_and_labels.csv"))

print("Reading overlap tests.")
overlap1 <- read.csv(paste0(set1_directory, "/overlap_test_fdr_05_RNASeq.csv"))
overlap2 <- read.csv(paste0(set2_directory, "/overlap_test_fdr_05_RNASeq.csv"))

print("Reading normalized tables.")
norm.table1 <- read.csv(paste0(set1_directory, "/normalized_table.csv"), row.names = 1)
norm.table2 <- read.csv(paste0(set2_directory, "/normalized_table.csv"), row.names = 1)

print("Reading reference counts per million.")
ref1 <- read.csv(paste0(set1_directory, "/ref_cpm.csv"), row.names = 1)
ref2 <- read.csv(paste0(set2_directory, "/ref_cpm.csv"), row.names = 1)

print("Finding intersection of genes if the first directory is used as a classifier.")
#FIND INTERSECT OF FEATURELABELS1 WITH OVERLAP2
test <- overlap2[overlap2$row.names.WNT_cpm. %in% colnames(feature_and_labels1), ]
test <- norm.table2[row.names(norm.table2) %in% test$row.names.WNT_cpm., ]

print("Creating corresponding feature and labels for matching test data.")
for (i in 1:nrow(test)){
  for (j in 1:ncol(test)){
    if (as.numeric(test[i,j]) <= as.numeric(ref2[1,j])){
      test[i,j] = 0
    }
    else if (as.numeric(test[i,j]) > as.numeric(ref2[1,j]) && as.numeric(test[i,j]) <= as.numeric(ref2[2,j])){
      test[i,j] = 1
    }
    else if (as.numeric(test[i,j]) > as.numeric(ref2[2,j]) && as.numeric(test[i,j]) <= as.numeric(ref2[3,j])){
      test[i,j] = 2
    }
    else if (as.numeric(test[i,j]) > as.numeric(ref2[3,j]) && as.numeric(test[i,j]) <= as.numeric(ref2[4,j])){
      test[i,j] = 3
    }
    else {
      test[i,j] = 4
    }
  }
}

print("Binding tables.")
groupings <- rep("group2",ncol(test))
groupings[1:set2_mutatedSampleSize] <- "group1"
test <- cbind(t(test),label=groupings)

train <- feature_and_labels1[,colnames(feature_and_labels1) %in% colnames(test)]

train_size <- nrow(train)
complete <- rbind(train, test)

print("Writing file if the first directory is the training set.")
groupings <- rep("set2",nrow(complete))
groupings[1:train_size] <- "set1"

complete <- cbind(complete,set=groupings)

write.csv(complete, paste0(set1_directory, "-asClassifier.csv"), row.names = F)

#THIS PROCESS IS REPEATING BUT INVERTED TO TRAIN IN THE OTHER DIRECTION. 
#THE BETTER OF THE TWO RESULTS WILL EVENTUALLY BE ACCEPTEd

#FIND INTERSECT OF FEATURELABELS2 WITH OVERLAP1
print("Finding intersection of genes if the second directory is used as a classifier.")
test <- overlap1[overlap1$row.names.WNT_cpm %in% colnames(feature_and_labels2), ]
test <- norm.table1[row.names(norm.table1) %in% test$row.names.WNT_cpm, ]

print("Creating corresponding feature and labels for matching test data.")
for (i in 1:nrow(test)){
  for (j in 1:ncol(test)){
    if (as.numeric(test[i,j]) <= as.numeric(ref1[1,j])){
      test[i,j] = 0
    }
    else if (as.numeric(test[i,j]) > as.numeric(ref1[1,j]) && as.numeric(test[i,j]) <= as.numeric(ref1[2,j])){
      test[i,j] = 1
    }
    else if (as.numeric(test[i,j]) > as.numeric(ref1[2,j]) && as.numeric(test[i,j]) <= as.numeric(ref1[3,j])){
      test[i,j] = 2
    }
    else if (as.numeric(test[i,j]) > as.numeric(ref1[3,j]) && as.numeric(test[i,j]) <= as.numeric(ref1[4,j])){
      test[i,j] = 3
    }
    else {
      test[i,j] = 4
    }
  }
}

print("Binding tables.")
groupings <- rep("group2",ncol(test))
groupings[1:set1_mutatedSampleSize] <- "group1"
test <- cbind(t(test),label=groupings)

train <- feature_and_labels2[,colnames(feature_and_labels2) %in% colnames(test)]

train_size <- nrow(train)
complete <- rbind(train, test)

print("Writing file if the second directory is the training set.")
groupings <- rep("set2",nrow(complete))
groupings[1:train_size] <- "set1"

complete <- cbind(complete,set=groupings)

write.csv(complete, paste0(set2_directory, "-asClassifier.csv"), row.names = F)
