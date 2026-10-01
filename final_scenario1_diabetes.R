################################################################################
# FINAL PROJECT SCENARIO 1: UNDERSTANDING DIABETES PROGRESSION
# Complete analysis script
#
# Upload this file AND diabetes.tab.txt to Posit Cloud (or run anywhere),
# set the working directory to the folder containing them, and run top-to-bottom.
################################################################################

# ---------------------------------------------------------------------------- #
# 0. SETUP
# ---------------------------------------------------------------------------- #

needed_packages <- c("ggplot2", "dplyr", "broom")
for (pkg in needed_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
library(ggplot2)
library(dplyr)
library(broom)

options(digits = 4)

# All plots are also saved as separate PNG files in a 'figures' folder
# next to this script, so they can be dropped into the report.
if (!dir.exists("figures")) dir.create("figures")


# ============================================================================ #
# QUESTION 1. UNDERSTAND AND PREPARE THE DATA
# ============================================================================ #

# One row = one patient (n = 442). Outcome: 'progression', a continuous measure
# of disease progression one year after baseline. Predictors: age, sex, bmi,
# bp (blood pressure) and six serum measurements s1-s6.

data_file <- "diabetes.tab.txt"
if (!file.exists(data_file)) {
  message("diabetes.tab.txt not found in working directory - please choose it.")
  data_file <- file.choose()
}

diabetes <- read.table(data_file, header = TRUE, sep = "\t")
names(diabetes) <- tolower(names(diabetes))
if ("y" %in% names(diabetes)) diabetes <- rename(diabetes, progression = y)

# 'sex' is coded 1/2 but is a CATEGORY, not a number -> convert to factor.
diabetes$sex <- factor(diabetes$sex)

dim(diabetes)          # 442 rows x 11 columns
names(diabetes)
str(diabetes)

colSums(is.na(diabetes))   # missing values per variable -> all zero
sum(is.na(diabetes))       # total missing -> 0
sum(duplicated(diabetes))  # duplicated rows -> 0

summary(diabetes)
# Data quality: no missing values, no duplicates, no impossible values
# (all measurements lie in plausible clinical ranges).
# Responsibility note: these are real patient records, so results describe
# associations in THIS sample - they must not be used to label or predict
# individuals, and interpretation must avoid over-stating what observational
# data can show.


# ============================================================================ #
# QUESTION 2. EXPLORE THE DATA
# ============================================================================ #

describe_numeric <- function(x) {
  c(
    n      = sum(!is.na(x)),
    mean   = mean(x, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    sd     = sd(x, na.rm = TRUE),
    min    = min(x, na.rm = TRUE),
    q1     = quantile(x, 0.25, na.rm = TRUE),
    q3     = quantile(x, 0.75, na.rm = TRUE),
    max    = max(x, na.rm = TRUE),
    iqr    = IQR(x, na.rm = TRUE)
  )
}

describe_numeric(diabetes$progression)  # centre ~139-140, roughly symmetric, wide spread
describe_numeric(diabetes$bmi)
describe_numeric(diabetes$age)

# Distribution of the outcome.
p1 <- ggplot(diabetes, aes(x = progression)) +
  geom_histogram(bins = 25, colour = "white") +
  labs(title = "Distribution of diabetes progression",
       x = "Progression measure", y = "Count") +
  theme_minimal()
print(p1)
ggsave("figures/fig1_progression_hist.png", p1, width = 6, height = 4, dpi = 300)

# Distribution of BMI.
p2 <- ggplot(diabetes, aes(x = bmi)) +
  geom_histogram(bins = 25, colour = "white") +
  labs(title = "Distribution of BMI", x = "BMI", y = "Count") +
  theme_minimal()
print(p2)
ggsave("figures/fig2_bmi_hist.png", p2, width = 6, height = 4, dpi = 300)

# Extra visualisation: progression by sex (boxplot).
p3 <- ggplot(diabetes, aes(x = sex, y = progression)) +
  geom_boxplot() +
  labs(title = "Diabetes progression by sex",
       x = "Sex (coded 1 or 2 in the source data)", y = "Progression measure") +
  theme_minimal()
print(p3)
ggsave("figures/fig3_sex_boxplot.png", p3, width = 5, height = 4, dpi = 300)

# Comment: progression is approximately symmetric (mean ~ median ~140) with no
# extreme outliers; BMI is slightly right-skewed; the two sex groups have
# broadly similar medians.


# ============================================================================ #
# QUESTION 3. INVESTIGATE RELATIONSHIPS
# ============================================================================ #

# sex is a factor now, so it is excluded from the Pearson correlations.
numeric_data <- select(diabetes, where(is.numeric))

cor_matrix <- cor(numeric_data, use = "complete.obs", method = "pearson")
round(cor_matrix, 2)

progression_correlations <- sort(cor_matrix[, "progression"], decreasing = TRUE)
print(progression_correlations)
# Strongest positive associations with progression: bmi (~0.59), s5 (~0.57),
# bp (~0.44), s4 (~0.43). s3 has a moderate NEGATIVE association (~-0.39).
# Age (~0.19), s1 (~0.21) and s2 (~0.17) are weakly related.

# Base-R heatmap: save with the png() device instead of ggsave().
png("figures/fig4_corr_heatmap.png", width = 1600, height = 1400, res = 220)
heatmap(cor_matrix, Rowv = NA, Colv = NA, scale = "none", symm = TRUE,
        margins = c(7, 7), main = "Correlation heatmap")
dev.off()
# Among predictors, s1 & s2 are almost perfectly correlated (~0.9) - they carry
# nearly the same information, a multicollinearity risk if both enter a model.
# s4 & s5, and s5 & s6, are also strongly related (~0.66, ~0.62).

# Scatterplot 1: BMI vs progression (strongest single predictor).
p4 <- ggplot(diabetes, aes(x = bmi, y = progression)) +
  geom_point(alpha = 0.65) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(title = "Progression versus BMI", x = "BMI", y = "Progression measure") +
  theme_minimal()
print(p4)
ggsave("figures/fig5_bmi_scatter.png", p4, width = 6, height = 4, dpi = 300)

# Scatterplot 2: blood pressure vs progression.
p5 <- ggplot(diabetes, aes(x = bp, y = progression)) +
  geom_point(alpha = 0.65) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(title = "Progression versus blood pressure",
       x = "Blood pressure", y = "Progression measure") +
  theme_minimal()
print(p5)
ggsave("figures/fig6_bp_scatter.png", p5, width = 6, height = 4, dpi = 300)

# Both relationships are positive, roughly linear, with moderate scatter.


# ============================================================================ #
# QUESTION 4. BUILD AND INTERPRET REGRESSION MODELS
# ============================================================================ #

# ---- Simple linear regression: progression ~ bmi -----------------------------
simple_model <- lm(progression ~ bmi, data = diabetes)
summary(simple_model)
tidy(simple_model, conf.int = TRUE)
glance(simple_model)
coef(simple_model)

# Fitted equation:  progression_hat = -117.77 + 10.23 * bmi
# Slope: each +1 BMI unit is associated with ~10.2 points more progression,
# on average (95% CI 8.9 to 11.6).
# Intercept (-117.8) is the predicted progression at bmi = 0, which is
# biologically meaningless - it only anchors the line.
# R-squared ~ 0.34: BMI alone explains about 34% of the variation in
# progression; ~66% is unexplained.


# ---- Multiple regression: bmi + bp + s5 --------------------------------------
# Chosen from the EDA: the three strongest, well-measured predictors; s1/s2
# were avoided because they are almost perfectly correlated with each other.
multiple_model <- lm(progression ~ bmi + bp + s5, data = diabetes)

summary(multiple_model)
tidy(multiple_model, conf.int = TRUE)
glance(multiple_model)

model_comparison <- bind_rows(
  glance(simple_model)    %>% mutate(model = "Simple: BMI"),
  glance(multiple_model)  %>% mutate(model = "Multiple: BMI + BP + S5")
) %>%
  select(model, r.squared, adj.r.squared, sigma)
print(model_comparison)

# Interpretation (holding the other predictors constant):
# - bmi: +1 unit -> ~ +6.5 points progression, holding bp and s5 fixed.
# - bp:  +1 unit -> ~ +0.9 points progression, holding bmi and s5 fixed.
# - s5:  +1 unit -> ~ +49.6 points progression, holding bmi and bp fixed.
# All three remain statistically significant; adjusted R-squared rises from
# ~0.34 (simple) to ~0.48 (multiple), so the extra predictors add real
# explanatory power - but most variation (~52%) is still unexplained, and a
# bigger R-squared does not by itself make the model "better" (simplicity,
# uncertainty and purpose matter too).


# ============================================================================ #
# QUESTION 5. EVALUATE THE ANALYSIS
# ============================================================================ #

# ---- Residuals versus fitted values ------------------------------------------
diagnostic_data <- augment(multiple_model)

p6 <- ggplot(diagnostic_data, aes(x = .fitted, y = .resid)) +
  geom_point(alpha = 0.65) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_smooth(se = FALSE) +
  labs(title = "Residuals versus fitted values",
       x = "Fitted values", y = "Residuals") +
  theme_minimal()
print(p6)
ggsave("figures/fig7_residuals_fitted.png", p6, width = 6, height = 4, dpi = 300)
# Residuals scatter randomly around zero with no strong curve or funnel shape:
# the linearity and constant-variance assumptions look reasonably satisfied.

# ---- Multicollinearity: VIF --------------------------------------------------
manual_vif <- function(model) {
  X <- model.matrix(model)[, -1, drop = FALSE]  # remove intercept
  vif_values <- sapply(seq_len(ncol(X)), function(i) {
    target <- X[, i]
    others <- X[, -i, drop = FALSE]
    if (ncol(others) == 0) return(1)
    r2 <- summary(lm(target ~ others))$r.squared
    1 / (1 - r2)
  })
  setNames(vif_values, colnames(X))
}

manual_vif(multiple_model)
# All VIFs are small (well under 5): bmi, bp and s5 are only moderately
# correlated with each other, so multicollinearity is not a serious problem
# for THIS model - even though s1/s2 (excluded) would have been.

# ---- Statistical uncertainty -------------------------------------------------
coef_table <- tidy(multiple_model, conf.int = TRUE)
print(coef_table)
confint(multiple_model, level = 0.95)

# Example discussion points:
# - bmi: estimate ~6.5, 95% CI [5.1, 7.9], p < 0.001 -> clearly different
#   from zero, and the effect is large enough to matter clinically.
# - bp: estimate ~0.90, 95% CI [0.48, 1.33], p < 0.001 -> statistically
#   significant; practical size is moderate (about 0.9 points per bp unit).
# Statistical significance (p < 0.05) is not the same as practical importance:
# with n = 442 even small effects can be significant, so the CI and the units
# of the estimate must be read together.
# Weaknesses: observational data, possible confounders not measured, ~52% of
# variation unexplained, and serum variables are composite/standardised
# measures that are hard to interpret clinically.


# ============================================================================ #
# QUESTION 6. DRAW APPROPRIATE CONCLUSIONS
# ============================================================================ #

# No new code. Conclusions supported by the evidence:
# - In THIS observational sample, higher BMI, higher blood pressure and higher
#   s5 are each associated with greater diabetes progression one year later.
# - The associations persist after adjusting for the other predictors and are
#   statistically clear.
#
# Claims that would GO BEYOND the evidence:
# - "BMI causes faster diabetes progression" - the data are observational;
#   confounders (diet, activity, genetics, medication, disease severity at
#   baseline) are not controlled, so no causal claim is justified.
# - "Reducing BMI will reduce progression" - an intervention claim, which would
#   require an experiment (randomised trial).
# - "These results apply to all diabetic patients" - the sample is one
#   population (standardised baseline measures); generalisation to other
#   populations is uncertain.
# - Prediction for a single individual - the model explains only ~48% of the
#   variation, far too little for individual-level prediction.


################################################################################
# END OF SCRIPT
################################################################################
