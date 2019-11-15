#####EdgeR###############
library(edgeR)

###### INPUT ARGUEMENTS#########
args <- commandArgs(trailingOnly = T)

inputFilePath <- args[1]
data <- args[2] #DATAYPE: RNASeq or microarray
isNormalized <- as.logical(args[3]) #IS THE INPUT DATA ALREADY NORMALIZED? True or False 
usesSymbols <- as.logical(args[4]) #DOES THE INPUT DATA USE SYMBOLS FOR GENE IDENTIFICATION. FALSE INDICATES IT USES ENSG IDENTIFICATION
PYTHON.PATH <- args[5]
group1.sampleSize <- as.numeric(args[6]) #NUMBER OF ELEMENTS IN MUTATED GROUP (FIRST GROUP). REMAINING ELEMENTS GO TO GROUP 2
n <- as.numeric(args[7]) #NUMBER OF GENES TO BE ANALYZED AT END OF SEQUENCE. INCLUDES TOP 'n' AND BOTTOM 'n'

#PRIMARY INPUT FILE. ASSUMES THAT FIRST N COLUMNS ARE OF 1 GROUP, AND REMAINING ARE STORED AS OTHER
print("--------------------------------------------------")
print("Step 1 - Calculate Stats")
print("--------------------------------------------------")
complete <- read.csv(inputFilePath, row.names=1)
print("Input File Read Successfuly - Dimensions of file:")
print(dim(complete))

#CONVERT ENSG TO SYMBOL IF APPLICABLE
if(usesSymbols == F){

  print("Converting ENSG to Gene IDs.")
  geneConversion <- read.csv("data/ensg2symbol.csv", row.names = 1, header = F)
  genes <- row.names(complete)
  newGenes <- c()

  for (i in genes){
    counter <- 0
    for (j in row.names(geneConversion)){
      if (i == j){
        newGenes <- c(newGenes, geneConversion[counter, 'V2'])
      }
      counter <- counter + 1
    }
  }

  row.names(complete) <- newGenes
}

#PROCESSING
if(data != "RNASeq"){
  print("Data already normalized. Normalized table complete.")
  norm.table <- complete
}else{
  print("Normalizing input table.")
  grp <- as.factor(rep(1, ncol(complete)))
  y <- DGEList(complete, group = grp)
  y <- calcNormFactors(y, method="TMM")
  norm.table <- cpm(y)
}

#GENERATES CORE STATISTICS FOR THE NORMALIZED TABLE -- USED IN DETERMINING FORMATION OF 5 GROUPS
print("Determining core metrics")
cov <- apply(norm.table, 1, sd)/apply(norm.table, 1, mean)
mean <- apply(norm.table, 1, mean)
std <- apply(norm.table, 1, sd)
MFC <- apply(norm.table, 1, max)/apply(norm.table, 1, min)

#WRITES THE NORMALIZED TABLE
write.csv(norm.table,"normalized_table.csv")
out <- data.frame(cov, mean, std, MFC)
out <- out[!is.infinite(out$MFC),]
out <- out[!is.na(out$MFC),] #REMOVES NA VALUES WHICH PRODUCE ERRORS IN READING

out <- cbind("genes"= rownames(out), out)

#BINDS WITH GENE NAMES
print("Writing core metrics file - final_out.csv")
write.csv(out, "final_out.csv", row.names=F)

print("--------------------------------------------------")
print("Step 2 - Find Reference Genes")
print("--------------------------------------------------")

print("Locating reference genes.")

refRange <- quantile(as.numeric(data.matrix(out)), probs = c(0.2, 0.4, 0.6, 0.8))

sink("data/referenceRanges.csv")
print(refRange)
sink()

#DATA FORMATTING
if(usesSymbols){
  reference_gene = system2(PYTHON.PATH, args="src/ref_sel_test_symbol.py",stdout=T)
} else{
  reference_gene = system2(PYTHON.PATH, args="src/ref_sel_test.py",stdout=T)
}

print(paste0("Reference genes: ",paste(reference_gene, collapse = ",")))
ref_cpm <- norm.table[reference_gene,]

#CREATES THE GUIDELINE COUNTS PER MILLION TO ASSOCIATE WITH 1-5 CATEGORIZATION
write.csv(ref_cpm, "ref_cpm.csv")
print("Saving a bin for the reference genes.")
save(reference_gene, file = "ref.gene.bin")

print("--------------------------------------------------")
print("Step 3 - Fisher Overlap Tests")
print("--------------------------------------------------")
if(data == "microarray"){
  TN_cpm <- norm.table[, 1:group1.sampleSize]
  OT_cpm <- norm.table[, (group1.sampleSize+1):ncol(norm.table)]
  
  #################Caterories################################
  TN_ref <- TN_cpm[reference_gene, ] 
  OT_ref <- OT_cpm[reference_gene, ] 
  
  ##############Counting distribution in each category#################
  ##################Computing Categorized No. #####################
  print("Creating groupings: 0% completed")
  TN_C0 <- c()
  for (i in 1:nrow(TN_cpm)){
    test <- c()
    for (j in 1:ncol(TN_cpm)){
      
      if(as.numeric(TN_cpm[i,j]) < as.numeric(TN_ref[1,j])){ 
        test[j] <- 1
      }
      else{
        test[j] <- 0
      }
    }
    TN_C0[i] <- sum(test)
  }
  
  print("Creating groupings: 8% completed")
  TN_C1 <- c()
  for (i in 1:nrow(TN_cpm)){
    test <- c()
    for (j in 1:ncol(TN_cpm)){
      
      if(as.numeric(TN_cpm[i,j]) >= as.numeric(TN_ref[1,j]) && as.numeric(TN_cpm[i,j])< as.numeric(TN_ref[2,j])){
        test[j] <- 1
      }
      else{
        test[j] <- 0
      }
    }
    TN_C1[i] <- sum(test)
  }
  
  print("Creating groupings: 16% completed")
  TN_C2 <- c()
  #out <- c()
  for (i in 1:nrow(TN_cpm)){
    test <- c()
    for (j in 1:ncol(TN_cpm)){
      
      if(as.numeric(TN_cpm[i,j]) >= as.numeric(TN_ref[2,j]) && as.numeric(TN_cpm[i,j])< as.numeric(TN_ref[3,j])){
        test[j] <- 1
      }
      else{
        test[j] <- 0
      }
    }
    TN_C2[i] <- sum(test)
  }
  
  print("Creating groupings: 25% completed")
  TN_C3 <- c()
  for (i in 1:nrow(TN_cpm)){
    test <- c()
    for (j in 1:ncol(TN_cpm)){
      
      if(as.numeric(TN_cpm[i,j]) >= as.numeric(TN_ref[3,j]) && as.numeric(TN_cpm[i,j])< as.numeric(TN_ref[4,j])){
        test[j] <- 1
      }
      else{
        test[j] <- 0
      }
    }
    TN_C3[i] <- sum(test)
  }
  
  print("Creating groupings: 33% completed")
  print(dim(TN_ref))
  TN_C4 <- c()
  for (i in 1:nrow(TN_cpm)){
    test <- c()
    for (j in 1:ncol(TN_cpm)){
      
      if(as.numeric(TN_cpm[i,j]) >= as.numeric(TN_ref[4,j]) && as.numeric(TN_cpm[i,j])< as.numeric(TN_ref[5,j])){
        test[j] <- 1
      }
      else{
        test[j] <- 0
      }
    }

    TN_C4[i] <- sum(test)
  }
  
  print("Creating groupings: 41% completed")
  TN_C5 <- c()
  for (i in 1:nrow(TN_cpm)){
    test <- c()
    for (j in 1:ncol(TN_cpm)){
      
      if(as.numeric(TN_cpm[i,j]) >= as.numeric(TN_ref[5,j])){
        test[j] <- 1
      }
      else{
        test[j] <- 0
      }
    }
    TN_C5[i] <- sum(test)
  }
  
  print("Creating groupings: 50% completed")
  OT_C0 <- c()
  for (i in 1:nrow(OT_cpm)){
    test <- c()
    for (j in 1:ncol(OT_cpm)){
      
      if(as.numeric(OT_cpm[i,j]) < as.numeric(OT_ref[1,j])){ 
        test[j] <- 1
      }
      else{
        test[j] <- 0
      }
    }
    OT_C0[i] <- sum(test)
  }
  
  print("Creating groupings: 58% completed")
  OT_C1 <- c()
  for (i in 1:nrow(OT_cpm)){
    test <- c()
    for (j in 1:ncol(OT_cpm)){
      
      if(as.numeric(OT_cpm[i,j]) >= as.numeric(OT_ref[1,j]) && as.numeric(OT_cpm[i,j])< as.numeric(OT_ref[2,j])){
        test[j] <- 1
      }
      else{
        test[j] <- 0
      }
    }
    OT_C1[i] <- sum(test)
  }
  
  print("Creating groupings: 66% completed")
  OT_C2 <- c()
  for (i in 1:nrow(OT_cpm)){
    test <- c()
    for (j in 1:ncol(OT_cpm)){
      
      if(as.numeric(OT_cpm[i,j]) >= as.numeric(OT_ref[2,j]) && as.numeric(OT_cpm[i,j])< as.numeric(OT_ref[3,j])){
        test[j] <- 1
      }
      else{
        test[j] <- 0
      }
    }
    OT_C2[i] <- sum(test)
  }
  
  print("Creating groupings: 75% completed")
  OT_C3 <- c()
  for (i in 1:nrow(OT_cpm)){
    test <- c()
    for (j in 1:ncol(OT_cpm)){
      
      if(as.numeric(OT_cpm[i,j]) >= as.numeric(OT_ref[3,j]) && as.numeric(OT_cpm[i,j])< as.numeric(OT_ref[4,j])){
        test[j] <- 1
      }
      else{
        test[j] <- 0
      }
    }
    OT_C3[i] <- sum(test)
  }
  
  print("Creating groupings: 83% completed")
  OT_C4 <- c()
  for (i in 1:nrow(OT_cpm)){
    test <- c()
    for (j in 1:ncol(OT_cpm)){
      
      if(as.numeric(OT_cpm[i,j]) >= as.numeric(OT_ref[4,j]) && as.numeric(OT_cpm[i,j])< as.numeric(OT_ref[5,j])){
        test[j] <- 1
      }
      else{
        test[j] <- 0
      }
    }
    OT_C4[i] <- sum(test)
  }
  
  print("Creating groupings: 91% completed")
  OT_C5 <- c()
  for (i in 1:nrow(OT_cpm)){
    test <- c()
    for (j in 1:ncol(OT_cpm)){
      
      if(as.numeric(OT_cpm[i,j]) >= as.numeric(OT_ref[5,j])){
        test[j] <- 1
      }
      else{
        test[j] <- 0
      }
    }
    OT_C5[i] <- sum(test)
  }
  
  print("Creating groupings: 100% completed")
  cat_count_total <- data.frame(row.names(TN_cpm), TN_C0, TN_C1, TN_C2, TN_C3, TN_C4,TN_C5, OT_C0, OT_C1, OT_C2, OT_C3, OT_C4,OT_C5)
  
  overlap_p <- c()
  for (i in 1: nrow(cat_count_total)){
    A <- as.numeric(cat_count_total[i,2:7])
    B <- as.numeric(cat_count_total[i,8:13])
    tab=as.table(rbind(A,B))
    row.names(tab)=c('TN','OT')
    c <- fisher.test(tab, workspace=2e+09,hybrid=TRUE)
    overlap_p[i] <- c$p.value
  }
  
  ED <- c()
  for (i in 1: nrow(cat_count_total)){
    A <- as.numeric(cat_count_total[i,2:7])
    B <- as.numeric(cat_count_total[i,8:13])
    D <- sum(A * seq_along(A))/ncol(TN_cpm)
    E <- sum(B * seq_along(B))/ncol(OT_cpm)
    c <- (D-E)
    ED[i] <- c
  }
  
  print("Performing overlap tests.")
  overlap_test <- data.frame(row.names(TN_cpm),overlap_p, ED )
  overlap_test_padjust <- p.adjust(overlap_p, method = "fdr", n = length(overlap_p))
  overlap_test_fdr <- data.frame(row.names(TN_cpm),overlap_test_padjust, ED, overlap_p)
  overlap_pvalue_fdr <- overlap_test_fdr[order(overlap_test_padjust), ]
  keep <- (overlap_pvalue_fdr$overlap_test_padjust <= 0.1) 
  write.csv(overlap_pvalue_fdr[keep, ], "overlap_test_fdr_1_RNASeq.csv")
  
  keep <- (overlap_pvalue_fdr$overlap_test_padjust <= 0.05) 
  write.csv(overlap_pvalue_fdr[keep, ], "overlap_test_fdr_05_RNASeq.csv")
} else if (data == "RNASeq") {
  WNT_cpm <- norm.table[, 1:group1.sampleSize]
  OT_cpm <- norm.table[, (group1.sampleSize+1):ncol(norm.table)]
  
  #################Caterories################################
  WNT_ref <- WNT_cpm[reference_gene, ] 
  OT_ref <- OT_cpm[reference_gene, ] 
  
  print("Creating groupings: 0% completed")
  WNT_C0 <- c()
  #out <- c()
  for (i in 1:nrow(WNT_cpm)){
    test <- c()
    for (j in 1:ncol(WNT_cpm)){
      
      if(as.numeric(WNT_cpm[i,j]) < as.numeric(WNT_ref[1,j])){ #&& as.numeric(WNT_cpm[i,j])<=as.numeric(WNT_ref[4,j])){
        test[j] <- 1
      }
      else{
        test[j] <- 0
      }
    }
    WNT_C0[i] <- sum(test)
  }
  
  print("Creating groupings: 10% completed")
  WNT_C1 <- c()
  #out <- c()
  for (i in 1:nrow(WNT_cpm)){
    test <- c()
    for (j in 1:ncol(WNT_cpm)){
      
      if(as.numeric(WNT_cpm[i,j]) >= as.numeric(WNT_ref[1,j]) && as.numeric(WNT_cpm[i,j])< as.numeric(WNT_ref[2,j])){
        test[j] <- 1
      }
      else{
        test[j] <- 0
      }
    }
    WNT_C1[i] <- sum(test)
  }
  
  print("Creating groupings: 20% completed")
  WNT_C2 <- c()
  #out <- c()
  for (i in 1:nrow(WNT_cpm)){
    test <- c()
    for (j in 1:ncol(WNT_cpm)){
      
      if(as.numeric(WNT_cpm[i,j]) >= as.numeric(WNT_ref[2,j]) && as.numeric(WNT_cpm[i,j])< as.numeric(WNT_ref[3,j])){
        test[j] <- 1
      }
      else{
        test[j] <- 0
      }
    }
    WNT_C2[i] <- sum(test)
  }
  
  print("Creating groupings: 30% completed")
  WNT_C3 <- c()
  #out <- c()
  for (i in 1:nrow(WNT_cpm)){
    test <- c()
    for (j in 1:ncol(WNT_cpm)){
      
      if(as.numeric(WNT_cpm[i,j]) >= as.numeric(WNT_ref[3,j]) && as.numeric(WNT_cpm[i,j])< as.numeric(WNT_ref[4,j])){
        test[j] <- 1
      }
      else{
        test[j] <- 0
      }
    }
    WNT_C3[i] <- sum(test)
  }
  
  print("Creating groupings: 40% completed")
  WNT_C4 <- c()
  #out <- c()
  for (i in 1:nrow(WNT_cpm)){
    test <- c()
    for (j in 1:ncol(WNT_cpm)){
      
      if(as.numeric(WNT_cpm[i,j]) >= as.numeric(WNT_ref[4,j])){#&& as.numeric(WNT_cpm[i,j])< as.numeric(WNT_ref[5,j])){
        test[j] <- 1
      }
      else{
        test[j] <- 0
      }
    }
    WNT_C4[i] <- sum(test)
  }
  
  print("Creating groupings: 50% completed")
  OT_C0 <- c()
  #out <- c()
  for (i in 1:nrow(OT_cpm)){
    test <- c()
    for (j in 1:ncol(OT_cpm)){
      
      if(as.numeric(OT_cpm[i,j]) < as.numeric(OT_ref[1,j])){ #&& as.numeric(OT_cpm[i,j])<=as.numeric(OT_ref[4,j])){
        test[j] <- 1
      }
      else{
        test[j] <- 0
      }
    }
    OT_C0[i] <- sum(test)
  }
  
  print("Creating groupings: 60% completed")
  OT_C1 <- c()
  #out <- c()
  for (i in 1:nrow(OT_cpm)){
    test <- c()
    for (j in 1:ncol(OT_cpm)){
      
      if(as.numeric(OT_cpm[i,j]) >= as.numeric(OT_ref[1,j]) && as.numeric(OT_cpm[i,j])< as.numeric(OT_ref[2,j])){
        test[j] <- 1
      }
      else{
        test[j] <- 0
      }
    }
    OT_C1[i] <- sum(test)
  }
  
  print("Creating groupings: 70% completed")
  OT_C2 <- c()
  #out <- c()
  for (i in 1:nrow(OT_cpm)){
    test <- c()
    for (j in 1:ncol(OT_cpm)){
      
      if(as.numeric(OT_cpm[i,j]) >= as.numeric(OT_ref[2,j]) && as.numeric(OT_cpm[i,j])< as.numeric(OT_ref[3,j])){
        test[j] <- 1
      }
      else{
        test[j] <- 0
      }
    }
    OT_C2[i] <- sum(test)
  }
  
  print("Creating groupings: 80% completed")
  OT_C3 <- c()
  #out <- c()
  for (i in 1:nrow(OT_cpm)){
    test <- c()
    for (j in 1:ncol(OT_cpm)){
      
      if(as.numeric(OT_cpm[i,j]) >= as.numeric(OT_ref[3,j]) && as.numeric(OT_cpm[i,j])< as.numeric(OT_ref[4,j])){
        test[j] <- 1
      }
      else{
        test[j] <- 0
      }
    }
    OT_C3[i] <- sum(test)
  }
  
  print("Creating groupings: 90% completed")
  OT_C4 <- c()
  #out <- c()
  for (i in 1:nrow(OT_cpm)){
    test <- c()
    for (j in 1:ncol(OT_cpm)){
      
      if(as.numeric(OT_cpm[i,j]) >= as.numeric(OT_ref[4,j])){#&& as.numeric(OT_cpm[i,j])< as.numeric(OT_ref[5,j])){
        test[j] <- 1
      }
      else{
        test[j] <- 0
      }
    }
    OT_C4[i] <- sum(test)
  }
  
  print("Creating groupings: 100% completed")
  cat_count_total <- data.frame(row.names(WNT_cpm), WNT_C0, WNT_C1, WNT_C2, WNT_C3, WNT_C4, OT_C0, OT_C1, OT_C2, OT_C3, OT_C4)

  print("Performing Fisher tests.") 
  overlap_p <- c()
  for (i in 1: nrow(cat_count_total)){
    A <- as.numeric(cat_count_total[i,2:6])
    B <- as.numeric(cat_count_total[i,7:11])
    #D <- as.numeric(cat_count_total[i,12:16])
    #E <- as.numeric(cat_count_total[i,17:21])
    tab=as.table(rbind(A,B))
    row.names(tab)=c('G4','OTHERS')
    c <- fisher.test(tab, workspace=2e+07,hybrid=TRUE)
    overlap_p[i] <- c$p.value
  }
  
  #####FOLD CHANGE########
  print("Calculating fold change.")
  ED <- c()
  for (i in 1: nrow(cat_count_total)){
    A <- as.numeric(cat_count_total[i,2:6])
    B <- as.numeric(cat_count_total[i,7:11])
    D <- sum(A * seq_along(A))/ncol(WNT_cpm)
    E <- sum(B * seq_along(B))/ncol(OT_cpm)
    #D <- as.numeric(cat_count_total[i,12:16])
    #E <- as.numeric(cat_count_total[i,17:21])
    c <- (D-E)
    ED[i] <- c
  }
  
  print("Performing overlap tests")
  overlap_test <- data.frame(row.names(WNT_cpm),overlap_p, ED )
  overlap_test_padjust <- p.adjust(overlap_p, method = "fdr", n = length(overlap_p))
  overlap_test_fdr <- data.frame(row.names(WNT_cpm),overlap_test_padjust, ED)
  overlap_pvalue_fdr <- overlap_test_fdr[order(ED, decreasing = TRUE), ]
  keep <- (overlap_pvalue_fdr$overlap_test_padjust <= 0.05)
  write.csv(overlap_pvalue_fdr[keep, ], "overlap_test_fdr_05_RNASeq.csv")
}

#FIND SIGNIFICANT GENES
print("Fingding list of significant genes")
geneList <- read.csv("overlap_test_fdr_05_RNASeq.csv", row.names = 1)
geneList <- geneList[order(geneList$ED), ]

print("--------------------------------------------------")
print("Step 4: Label Assignment")
print("--------------------------------------------------")

if (usesSymbols == F){
  sigLow <- as.character(head(geneList$row.names.TN_cpm., n = n))
  #sigLow <- geneList$row.names.TN_cpm.[0:n]

  sigHigh <- as.character(tail(geneList$row.names.TN_cpm., n = n))
  #sigHigh <- geneList$row.names.TN_cpm.[(length(geneList$ED) - (n-1)) : (length(geneList$ED))]
}

if (usesSymbols == T){
  sigLow <- as.character(head(geneList$row.names.WNT_cpm., n = n))
  #sigLow <- geneList$row.names.TN_cpm.[0:n]

  sigHigh <- as.character(tail(geneList$row.names.WNT_cpm., n = n))
  #sigHigh <- geneList$row.names.TN_cpm.[(length(geneList$ED) - (n-1)) : (length(geneList$ED))]
}

sigGenes <- c(as.character(sigLow), as.character(sigHigh))

#SUBSETS THE NORMALIZED TABLE INTO FEATURES WE CARE ABOUT
nameSearch <- row.names(norm.table)

feature <- norm.table[row.names(norm.table) %in% sigGenes, ]
#rownames(feature)
#sigGenes
#feature_old <- feature

#SEARCH THROUGH CATEGORIZATIONS TO FIND PROPER ASSIGNMENT
print("Assigning features to category")
for (i in 1:nrow(feature)){
  for (j in 1:ncol(feature)){
    if (as.numeric(feature[i,j]) <= as.numeric(ref_cpm[1,j])){
      feature[i,j] = 0
    }
    else if (as.numeric(feature[i,j]) > as.numeric(ref_cpm[1,j]) && as.numeric(feature[i,j]) <= as.numeric(ref_cpm[2,j])){
      feature[i,j] = 1
    }
    else if (as.numeric(feature[i,j]) > as.numeric(ref_cpm[2,j]) && as.numeric(feature[i,j]) <= as.numeric(ref_cpm[3,j])){
      feature[i,j] = 2
    }
    else if (as.numeric(feature[i,j]) > as.numeric(ref_cpm[3,j]) && as.numeric(feature[i,j]) <= as.numeric(ref_cpm[4,j])){
      feature[i,j] = 3
    }
    else {
      feature[i,j] = 4
    }
  }
}

#ADD LABELS OF GROUPINGS
print("Writing features and labels.")
groupings <- rep("group2",ncol(feature))
groupings[1:group1.sampleSize] <- "group1"

feature <- cbind(t(feature),label=groupings)
write.csv(feature, "feature_and_labels.csv",row.names = F)
