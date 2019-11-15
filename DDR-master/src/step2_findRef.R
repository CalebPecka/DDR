###### INPUT ARGUEMENTS#########
args <- commandArgs(trailingOnly = T)

PYTHON.PATH <- args[1]

print("--------------------------------------------------")
print("Step 2 - Find Reference Genes")
print("--------------------------------------------------")

print("Reading input files.")
norm.table <- read.csv("results/normalized_table.csv", row.names = 1)
out <- read.csv("final_out.csv", row.names = 1)

intermediate <- norm.table[row.names(norm.table) %in% row.names(out), ]

print("Determining reference range.")
refRange <- quantile(as.numeric(data.matrix(intermediate)), probs = c(0.2, 0.4, 0.6, 0.8))

sink("data/referenceRanges.csv")
print(refRange)
sink()

#DATA FORMATTING
reference_gene = system2(PYTHON.PATH, args="src/ref_sel_test_symbol.py",stdout=T)

print(paste0("Reference genes: ",paste(reference_gene, collapse = ",")))
print("*WARNING*: Fewer than 4 results will result in a future error. Modify range constant of ref_sel_test_symbols.py in this event.")
ref_cpm <- norm.table[reference_gene,]

#CREATES THE GUIDELINE COUNTS PER MILLION TO ASSOCIATE WITH 1-5 CATEGORIZATION
write.csv(ref_cpm, "results/ref_cpm.csv")
print("Saving a bin for the reference genes.")
save(reference_gene, file = "ref.gene.bin")
