load("data/lihc_raw_counts.RData")
raw_counts = lihc_raw_counts[,1:100]
dim(raw_counts)
DESeq2_math = log(raw_counts)
for(i in 1:length(raw_counts[,1])){
  for(j in 1:length(raw_counts[1,])){
    DESeq2_math[i,j] = log(raw_counts[i,j])
  }
}
# DESeq2_math

# Step 2 : Averages of rows

DESeq2_math_avg_gene = rowMeans(DESeq2_math)

# Step 3 : -Inf filtering

DESeq2_math_avg_gene = as.matrix(DESeq2_math_avg_gene)[!is.infinite(as.matrix(DESeq2_math_avg_gene)),]

DESeq2_math = DESeq2_math[!is.infinite(rowSums(DESeq2_math)),]

# Step 4 : average subtraction per gene

for(i in 1:length(DESeq2_math[,1])){
  for(j in 1:length(DESeq2_math[1,])){
    DESeq2_math[i,j] = DESeq2_math[i,j] - DESeq2_math_avg_gene[i]
  }
}

# Step 5 : Average per sample

DESeq2_math_avg_sample = colMeans(DESeq2_math)

# Step 6 : Convert column averages to decimal numbers to obtain scaling factors

for(i in 1:length(DESeq2_math_avg_sample)){
  DESeq2_math_avg_sample[i] = exp(DESeq2_math_avg_sample[i])
}
DESeq2_manual_scaling_factors = DESeq2_math_avg_sample
DESeq2_manual_scaling_factors

# Step 7 : Divide original counts by scaling factors

norm_DESeq2_counts = raw_counts
for(i in 1:length(raw_counts[,1])){
  for(j in 1:length(raw_counts[1,])){
    norm_DESeq2_counts[i,j] = norm_DESeq2_counts[i,j] / DESeq2_manual_scaling_factors[i]
  }
}

norm_DESeq2_counts
