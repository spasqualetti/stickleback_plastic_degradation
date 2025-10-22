# xps_permanova.R
suppressPackageStartupMessages({
  library(vegan)
  library(dplyr)
})

set.seed(42)

# ---------------------------
# Data
# ---------------------------
df <- data.frame(
  SampleID = c("241019_C1","241019_C2","241019_C3",
               "241019_195_1","241019_195_2","241019_195_3",
               "241019_191x195_1","241019_191x195_2","241019_191x195_3",
               "250418_C1","250418_C2","250418_C3",
               "250418_195_1","250418_195_2","250418_195_3",
               "250418_191x195_1","250418_191x195_2","250418_191x195_3"),
  Treatment = c(rep("Control",3),rep("KMM_195",3),rep("KMM_191x195",3),
                rep("Control",3),rep("KMM_195",3),rep("KMM_191x195",3)),
  `C-C/C-H` = c(80.71666667,81.58333333,85.17333333,61.06666667,66.30333333,81.1,83.78666667,81.73,82.86,
                79.27333333,79.32666667,83.27333333,72.62,72.94666667,78.17,76.47666667,74.29666667,74.35666667),
  `C-O`     = c(12.48333333,9.713333333,10.28666667,26.24666667,22.83666667,12.44666667,10.62333333,12.59,12.27666667,
                13.78,11.27,10.48333333,18.04666667,17.20666667,11.51333333,12.68666667,14.81666667,14.45),
  `C=O`     = c(2.086666667,2.57,2.763333333,9.863333333,8.943333333,4.153333333,3.803333333,3.87,3.01,
                4.633333333,1.97,1.84,1.313333333,7.096666667,2.143333333,1.3,1.81,2.593333333),
  `O-C=O`   = c(4.713333333,6.13,1.78,2.82,1.92,2.296666667,1.79,1.803333333,1.853333333,
                2.316666667,7.433333333,4.403333333,8.02,2.746666667,8.173333333,9.54,9.08,8.606666667),
  check.names = FALSE
)

rownames(df) <- df$SampleID
df$Treatment <- factor(df$Treatment, levels = c("Control","KMM_195","KMM_191x195"))

# Features (√ transform is common before Bray–Curtis)
features <- sqrt(df[, c("C-C/C-H","C-O","C=O","O-C=O")])
rownames(features) <- rownames(df)   # keep alignment
X <- features                         # <- define X

# ---------------------------
# Pairwise PERMANOVA
# ---------------------------
pairwise_permanova_all <- function(X, grouping,
                                   method = "bray",
                                   nperm = 9999,
                                   seed = 42,
                                   p_adjust = "BH") {
  set.seed(seed)
  grp <- droplevels(factor(grouping))
  levs <- levels(grp)
  pairs <- combn(levs, 2, simplify = FALSE)
  
  res <- lapply(pairs, function(pr) {
    keep <- grp %in% pr
    Xi   <- X[keep, , drop = FALSE]
    gi   <- droplevels(grp[keep])
    d    <- vegdist(Xi, method = method)
    dat  <- data.frame(grp = gi, row.names = rownames(Xi))
    a    <- adonis2(d ~ grp, data = dat, permutations = nperm)
    data.frame(Comparison = paste(pr, collapse = " vs "),
               F  = a$F[1],
               R2 = a$R2[1],
               p  = a$`Pr(>F)`[1],
               stringsAsFactors = FALSE)
  }) |> dplyr::bind_rows()
  
  res$p_adj <- p.adjust(res$p, method = p_adjust)
  res
}

pw <- pairwise_permanova_all(X, df$Treatment)
print(pw)
