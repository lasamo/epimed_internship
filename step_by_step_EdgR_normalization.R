load("data/lihc_raw_counts.RData")
raw_counts = lihc_raw_counts[,1:100]
dim(raw_counts)
scaling_factors = c()

## Step 1: Remove all 0 count genes

raw_counts = raw_counts[!apply(raw_counts == 0 , 1, all),]

# !!!
# apply(raw_counts, 1, function(l){
#   print(l)
#   sum(l)
#   }

dim(raw_counts)

## Step 2: Select reference sample

# Substep 1: Scale each sample to total read counts 

scaled_counts = raw_counts/colSums(raw_counts)
dim(scaled_counts)

# Substep 2: Calculate the 75th percentile

quantiles = c()
for(i in 1:ncol(scaled_counts)){
  quantiles[i] = quantile(scaled_counts[,i], .75)
}
names(quantiles) = colnames(scaled_counts)

# Substep 3: Select sample with quantile closest to quantile average

reference_sample = which.min(abs(quantiles-mean(quantiles)))
reference_sample_name = colnames(scaled_counts)[reference_sample]
print(paste0("Reference Sample: ", reference_sample_name, " | Index: ", reference_sample))

# Substep 1: calculate log2 of (reference sample/sample n)

for(i in 1:length(raw_counts[1,])){
  log_ratio = log2(scaled_counts[,reference_sample_name]/scaled_counts[,i])

    # Substep 2: filter out -Inf values
  
  log_ratio[log_ratio==-Inf] = NA
  log_ratio = na.omit(log_ratio)
  length(log_ratio)
  
  # Substep 3: calculate geometric mean
  
  gmeans = (log2(scaled_counts[,reference_sample_name])+log2(scaled_counts[,i]))/2
  
  # Substep 4: remove -Inf values
  
  gmeans[gmeans==-Inf] = NA
  gmeans = na.omit(gmeans)
  length(gmeans)
  
  # Substep 5: order tables from low to high
  # Substep 6 : remove top and bottom 30% of ratio table (log_ratio) and top and bottom 5% of means table (means)
  qs =  quantile(log_ratio, c(0.3,0.7))
  log_ratio = log_ratio[log_ratio > qs[1] & log_ratio < qs[2]]
  qs =  quantile(gmeans, c(0.05, 0.95))
  gmeans = gmeans[gmeans > qs[1] & gmeans < qs[2]]
  length(log_ratio)
  length(gmeans)
  
  # Substep 7 : Common genes between the tables are genes used for normalization
  
  genes = intersect(names(log_ratio), names(gmeans))
  length(genes)
  
  ## Step 4 : Calculate weighted average of the remaining log2 ratios with selected genes
  
  scaling_factors[i] = 2^(weighted.mean(log_ratio[genes], raw_counts[genes,i]))
  scaling_factors[i] = 2^(mean(log_ratio[genes]))
}
names(scaling_factors) = colnames(raw_counts)
scaling_factors[reference_sample_name] = 1
scaling_factors
