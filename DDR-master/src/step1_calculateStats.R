#####EdgeR###############
library(edgeR)

###### INPUT ARGUEMENTS#########
args <- commandArgs(trailingOnly = T)

inputFilePath <- args[1]
usesSymbols <- as.logical(args[2]) #DOES THE INPUT DATA USE SYMBOLS FOR GENE IDENTIFICATION. FALSE INDICATES IT USES ENSG IDENTIFICATION

#PRIMARY INPUT FILE. ASSUMES THAT FIRST N COLUMNS ARE OF 1 GROUP, AND REMAINING ARE STORED AS OTHER
print("--------------------------------------------------")
print("Step 1 - Calculate Stats")
print("--------------------------------------------------")
print("Reading input file.")
complete <- read.csv(inputFilePath, row.names = 1)
print("Input file read successfuly - Dimensions of file:")
print(dim(complete))

#CONVERT ENSG TO SYMBOL IF APPLICABLE
if(usesSymbols == F){
  print("Converting ENSG to Gene IDs.")
  geneConversion <- read.csv("data/ensg2symbol.csv", row.names = 1, header = F)
  genes <- row.names(complete)
  symbols <- geneConversion[genes,1]
  symbID <- symbols[-which(duplicated(symbols))] #REMOVES DUPLICATES
  complete <- complete[-which(duplicated(symbols)),]
  row.names(complete) <- symbID #ASSIGNS NEW GENE ROW NAMES
}

#PROCESSING
print("Normalizing input table.")
grp <- as.factor(rep(1, ncol(complete)))
y <- DGEList(complete, group = grp)
y <- calcNormFactors(y, method="TMM")
norm.table <- cpm(y)

#GENERATES CORE STATISTICS FOR THE NORMALIZED TABLE
print("Determining core metrics")
cov <- apply(norm.table, 1, sd)/apply(norm.table, 1, mean)
mean <- apply(norm.table, 1, mean)
std <- apply(norm.table, 1, sd)
MFC <- apply(norm.table, 1, max)/apply(norm.table, 1, min)

#WRITES THE NORMALIZED TABLE
write.csv(norm.table,"results/normalized_table.csv")
out <- data.frame(cov, mean, std, MFC)
out <- out[!is.infinite(out$MFC),]
out <- out[!is.na(out$MFC),] #REMOVES NA VALUES WHICH PRODUCE ERRORS IN READING

out <- cbind("genes"= rownames(out), out)

#BINDS WITH GENE NAMES
print("Writing finalized normal table - final_out.csv")
write.csv(out, "final_out.csv", row.names=F)
