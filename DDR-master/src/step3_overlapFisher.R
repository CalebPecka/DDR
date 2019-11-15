###### INPUT ARGUEMENTS#########
args <- commandArgs(trailingOnly = T)

group1.sampleSize <- as.numeric(args[1])

print("--------------------------------------------------")
print("Step 3 - Fisher Overlap Tests")
print("--------------------------------------------------")

print("Reading gene bin.")
load("ref.gene.bin")
print("Reading normalized table.")
norm.table <- read.csv("results/normalized_table.csv",row.names = 1)

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
  
print("Writing overlap tests")
overlap_test <- data.frame(row.names(WNT_cpm),overlap_p, ED )
overlap_test_padjust <- p.adjust(overlap_p, method = "fdr", n = length(overlap_p))
overlap_test_fdr <- data.frame(row.names(WNT_cpm),overlap_test_padjust, ED)
overlap_pvalue_fdr <- overlap_test_fdr[order(ED, decreasing = TRUE), ]
keep <- (overlap_pvalue_fdr$overlap_test_padjust <= 0.05)
write.csv(overlap_pvalue_fdr[keep, ], "results/overlap_test_fdr_05_RNASeq.csv")