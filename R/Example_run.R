# ---------- Step 3: Run it all ----------
# Load data and prepare
physeq_path <- system.file("extdata", "longdataexample.RDS", package = "eDNAModel")
physeq <- readRDS(physeq_path)

out1 <- prepare_long_data(physeq_obj = physeq)
long_df <- out1$long_df

# Add total reads for offset
long_df <- long_df %>%
  group_by(Sample) %>%
  mutate(total_reads = sum(y)) %>%
  ungroup()
long_df$log_total_reads <- log(long_df$total_reads)

# Subset 50 OTUs
otu_levels <- levels(long_df$i)[1:50]
otu_subset <- subset(long_df, i %in% otu_levels)
otu_subset <- droplevels(otu_subset)

# Fit model
out <- simulate_glm_burnin_iterations(
  data_glm = otu_subset,
  poisson_formula = y ~ Site + i + (1 | Sample) + (1 | Replicate) + offset(log_total_reads),
  binomial_formula = z_sim ~ Site + i,
  num_iterations   = 100,
  burn_in          = 50
)


# ---------- Step 4: Extract predictions ----------
# Occupancy predictions (psi)
psi_mat <- do.call(cbind, lapply(out$binomial_models, function(model) {
  predict(model, type = "response", newdata = otu_subset)
}))

# Abundance predictions (lambda)
lambda_mat <- do.call(cbind, lapply(out$poisson_models, function(model) {
  predict(model, type = "response", newdata = otu_subset)
}))

# Detection probability
p_detect_mat <- 1 - exp(-lambda_mat)

# ---------- Step 5: Summarize per OTU ----------
otu_summary <- data.frame(
  OTU           = otu_subset$i,
  psi_mean      = rowMeans(psi_mat),
  psi_se        = apply(psi_mat, 1, sd),
  psi_lwr       = apply(psi_mat, 1, quantile, probs = 0.025),
  psi_upr       = apply(psi_mat, 1, quantile, probs = 0.975),
  lambda_mean   = rowMeans(lambda_mat),
  lambda_se     = apply(lambda_mat, 1, sd),
  p_detect_mean = rowMeans(p_detect_mat),
  p_detect_se   = apply(p_detect_mat, 1, sd),
  p_detect_lwr  = apply(p_detect_mat, 1, quantile, probs = 0.025),
  p_detect_upr  = apply(p_detect_mat, 1, quantile, probs = 0.975)
) %>%
  group_by(OTU) %>%
  summarise(
    psi_mean        = mean(psi_mean),
    psi_se          = mean(psi_se),
    psi_lwr         = mean(psi_lwr),
    psi_upr         = mean(psi_upr),
    lambda_mean     = mean(lambda_mean),
    lambda_se       = mean(lambda_se),
    p_detect_mean   = mean(p_detect_mean),
    p_detect_se     = mean(p_detect_se),
    p_detect_lwr    = mean(p_detect_lwr),
    p_detect_upr    = mean(p_detect_upr),
    .groups = "drop"
  )

# ---------- View Results ----------
print(otu_summary)
