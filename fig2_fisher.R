##########################
# Fisher's test for Fig 2

mcma <- c(7, 1, 0) 
mcda <- c(1, 5, 1)
dcda <- c(1, 4, 2)


## Build 3x3 table
tbl <- rbind(DCDA = dcda, MCMA = mcma, MCDA = mcda)
colnames(tbl) <- c('full', 'sub', 'para')

##############################################
# Run Fisher–Freeman–Halton test for the total 3x3 table

result_3x3 = fisher.test(tbl)
# Fisher's Exact Test for Count Data
# 
# data:  tbl
# p-value = 0.01134
# alternative hypothesis: two.sided

##############################################
# 2 x 3 pairwise comparisons
## Pairs (2x3 subtables)
pairs <- list(
  c("MCMA", "MCDA"),
  c("MCMA", "DCDA"),
  c("MCDA", "DCDA")
)


ffh_exact <- function(M) {
  res <- fisher.test(M, simulate.p.value = FALSE, workspace = 2e8)
  return(res)
}

## get raw p-values
raw_p <- sapply(pairs, function(grp) {
  sub_tbl <- tbl[grp, , drop = FALSE]
  ffh_exact(sub_tbl)$p.value
})

# Bonferroni correction
adj_p <- p.adjust(raw_p, method = "bonferroni")

# summarize results
cmp_labels <- sapply(pairs, function(x) paste(x, collapse = " vs "))
result <- data.frame(
  Comparison = cmp_labels,
  P_value = signif(raw_p, 4),
  P_adj_Bonferroni = signif(adj_p, 4),
  row.names = NULL
)

print(result)
# Comparison P_value P_adj_Bonferroni
# 1 MCMA vs MCDA 0.01445          0.04336
# 2 MCMA vs DCDA 0.01445          0.04336
# 3 MCDA vs DCDA 1.00000          1.00000
# 


##############################################
# Amnionicity
## MCMA vs MCDA +DCDA
## Combine MCDA + DCDA
mcda_dcda <- mcda + dcda

## Build 2x3 table
tbl_combined <- rbind(
  MCMA = mcma,
  MCDA_DCDA = mcda_dcda
)
colnames(tbl_combined) <- c('full', 'sub', 'para')

print(tbl_combined)
#            full sub para
# MCMA         7   1    0
# MCDA_DCDA    2   9    3

## Run Fisher–Freeman–Halton exact test
p_amnion <- fisher.test(tbl_combined, simulate.p.value = FALSE, workspace = 2e8)
print(p_amnion)
# Fisher's Exact Test for Count Data
# 
# data:  tbl_combined
# p-value = 0.004728
# alternative hypothesis: two.sided

##############################################
# Chorionicity
## DCDA vs MCMA + MCDA
## Combine MCMA + MCDA
mcma_mcda <- mcma + mcda

## Build 2x3 table
tbl_combined <- rbind(
  MCMA_MCDA = mcma_mcda,
  DCDA = dcda
)
colnames(tbl_combined) <- c('full', 'sub', 'para')
print(tbl_combined)
#           full sub para
# MCMA_MCDA    8   6    1
# DCDA         1   4    2

## Run Fisher–Freeman–Halton exact test

p_chorion = fisher.test(tbl_combined, simulate.p.value = FALSE, workspace = 2e8)

print(p_chorion)
# Fisher's Exact Test for Count Data
# 
# data:  tbl_combined
# p-value = 0.1623
# alternative hypothesis: two.sided

p.adjust(c(p_chorion$p.value, p_amnion$p.value), method='bonferroni')
# > p.adjust(c(p_chorion$p.value, p_amnion$p.value), method='bonferroni')
# [1] 0.324584858 0.009456797
