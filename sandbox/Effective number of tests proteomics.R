# Step 1: Extract the relevant columns (exposures and outcomes) for PCA
exposures <- data[c(pfas_names_all, "score")]
outcomes <- data[c(prot_names)]

# Step 2: Perform PCA for exposures
pca_exposures <- prcomp(exposures)

# Step 4: Perform PCA for outcomes
pca_outcomes <- prcomp(outcomes)

# Step 5: Calculate eigenvalues for exposures
eigenvalues_exposures <- pca_exposures$sdev^2

# Step 6: Calculate eigenvalues for outcomes
eigenvalues_outcomes <- pca_outcomes$sdev^2

# Step 7: Determine the effective number of tests using the Kaiser-Guttman rule
M_eff_exposures <- sum(eigenvalues_exposures > 1)
M_eff_outcomes <- sum(eigenvalues_outcomes > 1)

# Step 8: Sum to get total number of effective tests
(M_eff_total <- M_eff_exposures+M_eff_outcomes)
# 42 tests
# #New p value
0.05/M_eff_total
0.05/M_eff_outcomes
