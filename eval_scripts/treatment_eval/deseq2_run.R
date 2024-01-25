library(DESeq2)

##########
### Run DESeq2
##########

# Read and format data
count_mat <- read.table("data/mouse_treatment/raw_counts.txt", header=TRUE)

response_data <- read.csv("data/mouse_treatment/responder_data.csv")
rownames(count_mat) <- count_mat$gene_id
colnames(count_mat) <- sub("^X", "", colnames(count_mat))
colnames(count_mat) <- sub("__", "_", colnames(count_mat))

response_data <- subset(response_data, response_data$Response != "4IgG") # Remove controls
response_data$Response[response_data$Response == "2S"] <- "1R" # Include stable disease in responders

rownames(response_data) <- response_data$SYMBOL
count_mat <- count_mat[, rownames(response_data)]

row_names <- rownames(count_mat)
col_names <- colnames(count_mat)

count_mat <- matrix(as.integer(unlist(count_mat)), nrow=length(row_names), ncol=length(col_names))
colnames(count_mat) <- col_names
rownames(count_mat) <- row_names

response_data$Response <- factor(response_data$Response, levels = c("3NR","1R")) # This means the logFold will be responder / non-responder

print(head(count_mat))

dds <- DESeqDataSetFromMatrix(countData = count_mat, colData = response_data, design = ~ Response)
print(dds)
res <- DESeq(dds)
print(res)
res <- results(res)
print(res)
write.csv(res, "results/mouse_treatment/differential_expression/deseq2_results.csv")

