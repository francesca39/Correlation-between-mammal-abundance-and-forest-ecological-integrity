# Load required libraries
library(dplyr)
library(readxl)
library(glmmTMB)
library(MuMIn)

# Load the dataset from Excel
file_path <- "Desktop/species/mountainhare3.xlsx"
data <- read_xlsx(file_path)

# Keep only observations with up to 5 snow days
# (Note: if snow_days can be 0, log(0) below will be -Inf and later removed by na.omit)
data <- filter(data, snow_days <= 5)

# Create/cast variables ---------------------------------------------

# Hierarchical identifiers as factors (nested random effects will use these)
data$triangles     <- as.factor(data$triangle_id)
data$year_f        <- as.factor(data$year)
data$riistakeskus  <- as.factor(data$riistakeskus)
data$rhy_id        <- as.factor(data$rhy_id)

# Center continuous covariates to improve interpretability/convergence
data$year.c            <- data$year - mean(data$year, na.rm = TRUE)
data$latitude.100km    <- (data$triangle_y - mean(data$triangle_y, na.rm = TRUE)) / 100000
data$longitude.100km   <- (data$triangle_x - mean(data$triangle_x, na.rm = TRUE)) / 100000
data$FLII              <- data$indexvalue - mean(data$indexvalue, na.rm = TRUE)
data$forest_cover      <- data$forest_pro - mean(data$forest_pro, na.rm = TRUE)
data$canopy            <- data$canopy - mean(data$canopy, na.rm = TRUE)

# Log transforms
# Caution: if snow_days == 0, log(0) = -Inf
data$ln.snow.days  <- log(data$snow_days)

# Offset: log of surveyed area/effort (divide by 10,000 to rescale units if needed)
data$ln.covered.di <- log(data$covered_di / 10000)

# Remove rows with any missing or infinite values
data <- na.omit(data)

# Fit Negative Binomial GLMM (nbinom2 with log link)
# Fixed effects: FLII, forest cover, canopy, latitude, year (factor), log snow days
# Random effects: triangles nested in rhy_id nested in riistakeskus
# Offset: log(surveyed area/effort)
model_nb <- glmmTMB(
  n_tracks ~ FLII + forest_cover + canopy + latitude.100km + year_f + ln.snow.days +
    (1 | riistakeskus / rhy_id / triangles),
  offset = ln.covered.di,
  family = nbinom2(link = "log"),
  data = data
)

# Print model summary
print(summary(model_nb))

# Check whether the offset variable is present in the model frame
if ("ln.covered.di" %in% colnames(model_nb$frame)) {
  print("ln.covered.di is present in the model's data frame.")
} else {
  print("ln.covered.di is NOT present in the model's data frame.")
}

# Custom function to compute R2 on the link scale (marginal and conditional)
# - incl_DSV: include distribution-specific variance term
# - incl_offset: if FALSE, subtract the offset from fixed-effect predictions
# - offset_var: name of the offset column
# - exclude: character vector of fixed-effect columns to zero-out when predicting
glmmTMB.R2 <- function(m, incl_DSV = FALSE, incl_offset = FALSE, offset_var = NULL, exclude = NULL) {
  # Use the model frame as prediction data
  new.data.m <- m$frame

  # If requested, ensure the offset column is present in newdata
  if (!is.null(offset_var) & incl_offset == FALSE) {
    new.data.m[[offset_var]] <- data[[offset_var]]
  }

  # Optionally zero-out selected fixed-effect columns (to exclude their contribution)
  if (!is.null(exclude)) {
    new.data.m[, match(exclude, colnames(new.data.m))] <- 0
  }

  # Predicted values on link (log) scale with fixed effects only
  ln.yhat <- predict(m, type = "link", re.form = ~0, newdata = new.data.m)
  n <- length(ln.yhat)

  # Optionally remove the offset contribution from fixed-effect predictions
  if (!is.null(offset_var) & incl_offset == FALSE) {
    ln.yhat <- ln.yhat - new.data.m[[offset_var]]
  }

  # Variance components -------------------------------------------------

  # Fixed-effects variance (use population variance, n denominator)
  s2_f <- var(ln.yhat) * (n - 1) / n

  # Random-effects variances (sum across all random terms)
  s2_a <- unlist(VarCorr(m)$cond)

  # Residual/dispersion component on link scale (approximation for NB)
  # Note: this uses trigamma(sigma(m)); alternative formulations exist in the literature.
  s2_e <- trigamma(sigma(m))

  # Optional distribution-specific variance term (DSV) on link scale
  # Here approximated from the mean of ln.yhat
  s2_d <- log(1 / exp(mean(ln.yhat)) + 1)

  # Total variance and R2 metrics
  if (incl_DSV == TRUE) {
    s2_tot <- s2_f + sum(s2_a) + s2_e + s2_d
    R2_marg <- s2_f / s2_tot
    R2_cond <- (s2_f + sum(s2_a)) / s2_tot
  } else {
    s2_tot <- s2_f + sum(s2_a) + s2_e
    R2_marg <- s2_f / s2_tot
    R2_cond <- (s2_f + sum(s2_a)) / s2_tot
  }

  # Return named vector with marginal and conditional R2
  R2 <- c(R2_marg = unname(R2_marg), R2_cond = unname(R2_cond))
  return(R2)
}

# Compute custom R2 (excluding ln.snow.days from fixed-effects contribution and
# subtracting the offset from predictions)
r2_custom <- glmmTMB.R2(
  model_nb,
  incl_offset = FALSE,
  offset_var = "ln.covered.di",
  exclude = "ln.snow.days"
)
print(r2_custom)

# Compute marginal/conditional R2 via MuMIn (Nakagawa & Schielzeth)
cat("R2 MuMIn:\n")
print(r.squaredGLMM(model_nb))

