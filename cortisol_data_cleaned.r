# load packages
library(ggplot2)
library(dplyr)
library(svglite)
library(emmeans)  
library(car)
library(multcomp)
library(glmmTMB)
library(DHARMa)
library(boot)
library(tidyr)

# set working directory
setwd("/Users/marcbauhus/Library/Mobile Documents/com~apple~CloudDocs/PhD/Astyanax/field stuff/cortisol measurements/manuscript/Environmental Biology of Fishes/revision/data")

# import cortisol data - cortisol_data_wild.csv (model1, 2), cortisol_data_labvswild.csv (model3)
cdata <- read.csv("cortisol_data_wild.csv", header = TRUE, sep = ",", dec = ".")
head(cdata)
cdata

# check if cortisol is normally distributed, if yes - symmetrical bell shape without tail
# histogram
ggplot(cdata, aes(x = cortisol)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black", alpha = 0.7) +
  geom_density(color = "red", size = 1) +
  labs(title = "Histogram of Cortisol Levels", x = "Cortisol (pg/mg)", y = "Density") +
  theme_minimal(
) # doesn't look like normal distribution!!!

# shapiro-wilk test, if p < 0.05 data is not normally distributed
shapiro.test(cdata$cortisol) # in this case not normally distributed!!!



### Fit the generalized-linear models
# Prepare data and fit Models 1 and 2

# Convert categorical variables to factors BEFORE modelling
cdata$habitat <- factor(cdata$habitat)
cdata$population <- factor(cdata$population)
cdata$year <- factor(cdata$year)
cdata$sex <- factor(cdata$sex)

# Centre handling time (and standard length) around the mean of the analysis dataset
cdata$time_capture_dissection_c <- cdata$time_capture_dissection_hours - mean(cdata$time_capture_dissection_hours)

### Model 1: habitat effect, accounting for population
model1 <- glmmTMB(cortisol ~ habitat + time_capture_dissection_c + k_index + sex +  (1 | population),
  data = cdata,
  family = Gamma(link = "log")
)

# # alternative model 1, using standard_length instead of k_index as covariate
# model1 <- glmmTMB(cortisol ~ habitat + time_capture_dissection_c + standard_length + sex +  (1 | population),
#   data = cdata,
#   family = Gamma(link = "log")
# )

# model summary
summary(model1)

# save output
capture.output(
  summary(model1),
  file = "summary_model1_mixed.txt"
)

# Inspect the estimated population variance
VarCorr(model1)

# Simulation-based residual diagnostics
svglite("model1.svg", width = 8, height = 10)
res_model1 <- simulateResiduals(model1)
plot(res_model1)
dev.off()

testUniformity(res_model1)
testDispersion(res_model1)

# check if residuals are normally-distributed
# pull residuals from the model (type in model of interest)
resid_standard <- residuals(model1)
resid_standard

# shapiro-wilk test, if p < 0.05 residuals are not normally distributed
shapiro.test(resid_standard) # residuals not normally distributed!!!


## pairwise comparisons - type in model of interest
# get population estimated marginal means to conduct pairwise comparisons
population_emmeans <- emmeans(model1, pairwise ~ habitat, adjust = "tukey")
summary(population_emmeans)
summary(population_emmeans, infer = TRUE)

# Get the compact letter display (CLD)
cld_results <- cld(population_emmeans, Letters = letters)

# View the CLD results
print(cld_results)

write.csv(cld_results, "clds_model1.csv", row.names = FALSE)

# Extract the pairwise comparison table
pairwise_results <- as.data.frame(population_emmeans$contrasts)

# Add significance markers based on p-value thresholds
pairwise_results$significance <- ifelse(pairwise_results$p.value < 0.001, "***",
                                ifelse(pairwise_results$p.value < 0.01, "**",
                                ifelse(pairwise_results$p.value < 0.05, "*", "ns")))

# View the pairwise comparison results to check significance
head(pairwise_results)

# Write the results to a CSV file
write.csv(pairwise_results, "pairwise_comparisons_model1.csv", row.names = FALSE)


## calculate EMMs and plot data for model 1 
# 1. Calculate EMMs for habitat
emmeans_table <- emmeans(model1, ~ habitat)
head(emmeans_table)

# 2. Back-transform emmeans from log to original scale
emmeans_table <- as.data.frame(emmeans_table)  # Convert to data frame
emmeans_table$exp_emmean <- exp(emmeans_table$emmean)
emmeans_table$exp_lower.CL <- exp(emmeans_table$asymp.LCL)
emmeans_table$exp_upper.CL <- exp(emmeans_table$asymp.UCL)
emmeans_table$exp_SE <- exp(emmeans_table$SE)

# Optional: Rename habitat column for consistent plotting
names(emmeans_table)[names(emmeans_table) == "habitat"] <- "population"
head(emmeans_table)
# 3. Write EMMs to CSV
write.csv(emmeans_table, "emmeans_table_transformed_model1.csv", row.names = FALSE)

# 4. Filter NAs in cortisol (if not done already)
cdata <- cdata %>% filter(is.finite(cortisol))

# Set the order of habitat levels
cdata$habitat <- factor(cdata$habitat, levels = c("river", "lake", "cave"))
emmeans_table$population <- factor(emmeans_table$population, levels = c("river", "lake", "cave"))

# 5. Violin plot of raw data + back-transformed EMMs
svglite("cortisol_data_raw_back-transformed-emm_model1.svg", width = 7, height = 10)
ggplot() +
  geom_violin(data = cdata, aes(x = habitat, y = cortisol, fill = habitat), 
              trim = FALSE, alpha = 0.5, color = "black", width = 0.7) +
  geom_jitter(data = cdata, aes(x = habitat, y = cortisol, fill = habitat), 
              size = 1.7, alpha = 0.5, position = position_jitter(width = 0.13)) +
  geom_point(data = emmeans_table, aes(x = population, y = exp_emmean), 
             color = "black", size = 5, shape = 16, stroke = 2) +
  geom_errorbar(data = emmeans_table, 
                aes(x = population, ymin = exp_lower.CL, ymax = exp_upper.CL), 
                width = 0, color = "black", linewidth = 2) +
  labs(title = "Cortisol Levels by Habitat",
       x = "Habitat", y = "Cortisol (pg/mg scales)") +
  scale_fill_manual(values = c("cave" = "black", "lake" = "white", "river" = "gray")) +
  scale_y_continuous(breaks = seq(0, max(cdata$cortisol, na.rm = TRUE), by = 20)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 20, angle = 45, hjust = 1),
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 22),
    legend.position = "none",
    plot.title = element_text(size = 20, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90")
  )
dev.off()


## plot correlation between cortisol and body condition with model fitted trend line that is based on predictions for model 1
# Create a new dataset for prediction over k_index, averaging over population
new_data <- data.frame(
  k_index = seq(min(cdata$k_index, na.rm = TRUE), 
                max(cdata$k_index, na.rm = TRUE), 
                length.out = 100),
  sex = factor("m", levels = levels(cdata$sex)),
  time_capture_dissection_c = 0,
  population = factor(
    levels(cdata$population)[1],
    levels = levels(cdata$population)
  )
)

# Get model matrix for population (to compute average effect)
habitat_levels <- unique(cdata$habitat)
hab_effects <- sapply(habitat_levels, function(hab) {
  new_data$habitat <- hab
  predict(model1, newdata = new_data, type = "link")  # Get log-scale predictions
})

# Compute the average log-scale prediction across populations
new_data$log_predicted_cortisol <- rowMeans(hab_effects)

# Convert back to the response scale
new_data$predicted_cortisol <- exp(new_data$log_predicted_cortisol)

# Get standard errors for confidence intervals (approximated by variance across population effects)
log_se <- apply(hab_effects, 1, sd)  # Standard deviation across populations
new_data$lower_ci <- exp(new_data$log_predicted_cortisol - 1.96 * log_se)
new_data$upper_ci <- exp(new_data$log_predicted_cortisol + 1.96 * log_se)

# Plot the data and the single model-based trend line
svglite("cortisol_vs_k-index_model1.svg", width = 8, height = 10)
ggplot(cdata, aes(x = k_index, y = cortisol)) +
  geom_point(aes(fill = habitat), shape = 21, color = "black", stroke = 1.2, size = 4, alpha = 0.9) +  # Raw data points
  scale_fill_manual(values = c("cave" = "black", "river" = "gray", "lake" = "white")) +
  geom_ribbon(data = new_data, aes(x = k_index, ymin = lower_ci, ymax = upper_ci), 
              fill = "blue", alpha = 0.1, inherit.aes = FALSE) +  # Shaded confidence interval
  geom_line(data = new_data, aes(x = k_index, y = predicted_cortisol), 
            color = "blue", linewidth = 1) +  # Model-based trend line
  labs(x = "Body Condition (k-index)", 
       y = "Cortisol (pg/mg scales)", 
       title = "Cortisol vs. Body Condition (Single Model-Based Trend Line)") +
  theme_minimal() +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20, face = "bold"),
        legend.text = element_text(size = 18), 
        legend.title = element_text(size = 20, face = "bold"))
dev.off()

## plot correlation between cortisol and handling time with model fitted trend line that is based on predictions for model 1

# Dataset used for plotting; should match the model-fitting data
plot_data <- cdata

# Generate a sequence of centred handling times
time_seq_c <- seq(
  min(plot_data$time_capture_dissection_c, na.rm = TRUE),
  max(plot_data$time_capture_dissection_c, na.rm = TRUE),
  length.out = 100
)

# Hold K-index at its mean and sex at male.
# Average equally over habitat levels.
emm_time <- emmeans(
  model1,
  ~ time_capture_dissection_c,
  at = list(
    time_capture_dissection_c = time_seq_c,
    k_index = mean(plot_data$k_index, na.rm = TRUE),
    sex = "m"
  ),
  weights = "equal",
  type = "response"
)

# Convert predictions to a data frame
new_data <- as.data.frame(emm_time)

new_data$predicted_cortisol <- new_data$response

# Identify confidence interval columns
lower_col <- intersect(
  c("asymp.LCL", "lower.CL"), names(new_data)
)[1]

upper_col <- intersect(
  c("asymp.UCL", "upper.CL"), names(new_data)
)[1]

new_data$lower_ci <- new_data[[lower_col]]
new_data$upper_ci <- new_data[[upper_col]]

# Convert centred time back to actual hours for plotting
time_center_constant <- mean(
  plot_data$time_capture_dissection_hours -
    plot_data$time_capture_dissection_c,
  na.rm = TRUE
)

new_data$time_capture_dissection_hours <-
  new_data$time_capture_dissection_c + time_center_constant

# Plot the data and model-based trend line
svglite("cortisol_vs_handling_time_model1.svg", width = 8, height = 10)
ggplot(plot_data, aes(x = time_capture_dissection_hours, y = cortisol)) +
  geom_point(aes(fill = habitat), shape = 21, color = "black", stroke = 1.2, size = 4, alpha = 0.9) +
  scale_fill_manual(values = c("cave" = "black", "river" = "gray", "lake" = "white")) +
  geom_ribbon(data = new_data, aes(x = time_capture_dissection_hours, ymin = lower_ci, ymax = upper_ci), fill = "blue", alpha = 0.1, inherit.aes = FALSE) +
  geom_line(data = new_data, aes( x = time_capture_dissection_hours, y = predicted_cortisol), color = "blue", linewidth = 1, inherit.aes = FALSE) +
  labs(x = "Time between capture and dissection (hours)", y = "Cortisol (pg/mg scales)", title = "Cortisol vs. Capture-to-Dissection Interval") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 20, face = "bold")
  )
dev.off()


## calculate EMMs and plot sex differences for model 1
# Calculate EMMs for sex
emm_sex <- emmeans(model1, ~ sex, type = "response", weights = "flat")  # type = "response" gives results on the original scale

# Convert to data frame for plotting
emmeans_table <- as.data.frame(emm_sex)

emmeans_table$exp_emmean <- emmeans_table$response
emmeans_table$exp_lower.CL <- emmeans_table$asymp.LCL
emmeans_table$exp_upper.CL <- emmeans_table$asymp.UCL

emmeans_table

# Write the table to a CSV file
write.csv(emmeans_table, "emmeans_table_sex_model1.csv", row.names = FALSE)

cdata_clean <- subset(cdata, !is.na(sex))

# plot
svglite("cortisol_sex_model1.svg", width = 8, height = 10)
ggplot() +
  geom_violin(data = cdata_clean, aes(x = sex, y = cortisol, fill = sex),
              trim = FALSE, alpha = 0.5, color = "black", width = 1.0) +
  geom_jitter(data = cdata_clean, aes(x = sex, y = cortisol, fill = sex),
              size = 1.7, alpha = 0.5, position = position_jitter(width = 0.13)) +
  geom_point(data = emmeans_table, aes(x = sex, y = exp_emmean),
             color = "black", size = 5, shape = 16, stroke = 2) +
  geom_errorbar(data = emmeans_table,
                aes(x = sex, ymin = exp_lower.CL, ymax = exp_upper.CL),
                width = 0, color = "black", linewidth = 2) +
  labs(title = "Cortisol Levels by Sex",
       x = "Sex", y = "Cortisol (pg/mg scales)") +
  scale_fill_manual(values = c("f" = "red", "m" = "lightblue")) +
  scale_y_continuous(breaks = seq(0, max(cdata$cortisol, na.rm = TRUE), by = 20)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 20, angle = 45, hjust = 1),
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 22),
    legend.position = "none",
    plot.title = element_text(size = 20, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90")
  )
dev.off()



### model 2 - checks for population x year effects on cortisol levels
model2 <- glm(cortisol ~ population * year + k_index + sex + time_capture_dissection_c, data=cdata, family = Gamma(link = "log"), na.action = na.exclude) 

# # alternative model 2, using standard_length instead of k_index as covariate
# model2 <- glm(cortisol ~ population * year + standard_length + sex + time_capture_dissection_c, data=cdata, family = Gamma(link = "log"), na.action = na.exclude) 

# Summary of the model 
summary(model2)

capture.output(
  summary(model2),
  file = "summary_model2.txt"
)

# check if residuals are normally-distributed
# pull residuals from the model 
resid_standard <- residuals(model2)
resid_standard

# shapiro-wilk test, if p < 0.05 residuals are not normally distributed
shapiro.test(resid_standard) # residuals not normally distributed!!!

# test for pseudo R2 (1-(residual deviance/null deviance)) - variation explained
1-(35.522/66.821) # model2 = 0.47

# model diagnostics (type in model of interest)
svglite("model2.svg", width = 6, height = 10)
par(mfrow = c(2, 2)) # Set up plotting area for multiple plots
plot(model2)
dev.off()


## pairwise comparisons 
# get population estimated marginal means to conduct pairwise comparisons
population_emmeans <- emmeans(model2, pairwise ~ population * year, adjust = "tukey")
summary(population_emmeans)
summary(population_emmeans, infer = TRUE)

# Get the compact letter display (CLD)
cld_results <- cld(population_emmeans, Letters = letters)

# View the CLD results
print(cld_results)

write.csv(cld_results, "clds_model2.csv", row.names = FALSE)

# Extract the pairwise comparison table
pairwise_results <- as.data.frame(population_emmeans$contrasts)

# Add significance markers based on p-value thresholds
pairwise_results$significance <- ifelse(pairwise_results$p.value < 0.001, "***",
                                ifelse(pairwise_results$p.value < 0.01, "**",
                                ifelse(pairwise_results$p.value < 0.05, "*", "ns")))

# View the pairwise comparison results to check significance
head(pairwise_results)

# Write the results to a CSV file
write.csv(pairwise_results, "pairwise_comparisons_model2.csv", row.names = FALSE)


## compute and back-transform EMMs by population * year, model 2, plot data
# Calculate estimated marginal means for each population × year combination
emmeans_table2 <- emmeans(model2, ~ population * year)

# Back-transform from log to original scale
emmeans_table2 <- as.data.frame(emmeans_table2)
emmeans_table2$exp_emmean <- exp(emmeans_table2$emmean)
emmeans_table2$exp_lower.CL <- exp(emmeans_table2$lower.CL)
emmeans_table2$exp_upper.CL <- exp(emmeans_table2$upper.CL)
emmeans_table2$exp_SE <- exp(emmeans_table2$SE)

emmeans_table2
emmeans_table2 <- emmeans_table2 %>% filter(is.finite(exp_emmean))

# Write the table to a CSV file
write.csv(emmeans_table2, "emmeans_table_transformed_model2.csv", row.names = FALSE)

# Create a unified column for plotting on the x-axis
emmeans_table2$population_year <- paste0(emmeans_table2$population, "_", substr(emmeans_table2$year, 3, 4))

# Make sure cortisol is finite
cdata <- cdata %>% filter(is.finite(cortisol))

# Create matching population_year column in cdata
cdata$population_year <- paste0(cdata$population, "_", substr(cdata$year, 3, 4))

# set plotting order of population_years
population_order <- c(  
  "Rio_Choy_23", "Rio_Choy_24", "Presa_El_Oyul_23", "Presa_El_Oyul_24", "Pachón_23", "Los_Sabinos_23", "Tinaja_23", "Tinaja_24", "Rio_Subterráneo_23", "Rio_Subterráneo_24" 
)

# Save SVG output
svglite("cortisol_data_raw_back-transformed-emm_model2.svg", width = 8, height = 10)
ggplot() +
  geom_violin(data = cdata, aes(x = factor(population_year, levels = population_order), y = cortisol, fill = population_year),
              trim = FALSE, alpha = 0.5, color = "black", width = 0.7) +
  geom_jitter(data = cdata, aes(x = factor(population_year, levels = population_order), y = cortisol, fill = population),
              size = 1.7, alpha = 0.5, position = position_jitter(width = 0.13)) +
  geom_point(data = emmeans_table2, aes(x = population_year, y = exp_emmean),
             color = "black", size = 4, shape = 16, stroke = 2) +
  geom_errorbar(data = emmeans_table2,
                aes(x = population_year, ymin = exp_lower.CL, ymax = exp_upper.CL),
                width = 0, color = "black", linewidth = 1.5) +
  labs(title = "Cortisol Levels by Wild Population and Sampling Year",
       x = "Population-Year", y = "Cortisol (pg/mg scales)") +
  scale_fill_manual(values = c(
    "Rio_Choy_23" = "#ff0000ff", "Rio_Choy_24" = "#ff0000ff",
    "Presa_El_Oyul_23" = "#6c006cff", "Presa_El_Oyul_24" = "#6c006cff",
    "Pachón_23" = "#131300ff", "Los_Sabinos_23" = "#ffcc00ff",
    "Tinaja_23" = "#ff2e8eff", "Tinaja_24" = "#ff2e8eff",
    "Rio_Subterráneo_23" = "#0000ffff", "Rio_Subterráneo_24" = "#0000ffff"
  )) +
  scale_y_continuous(breaks = seq(0, max(cdata$cortisol, na.rm = TRUE), by = 20)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 16, angle = 45, hjust = 1),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.position = "none",
    plot.title = element_text(size = 20, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90")
  )
dev.off()


## plot correlation between cortisol and body condition with model fitted trend line that is based on predictions for model 2
cdata$pop_year <- interaction(cdata$population, cdata$year, sep = "_")

# Define the range of k_index
new_data_base <- data.frame(
  k_index = seq(min(cdata$k_index, na.rm = TRUE),
                max(cdata$k_index, na.rm = TRUE),
                length.out = 100),
  sex = factor("m", levels = levels(cdata$sex)),
  time_capture_dissection_c = 0
)

# Get all unique population-year combinations
pop_year_levels <- as.character(unique(cdata$pop_year))

# Prepare a matrix to hold predictions
log_preds_matrix <- matrix(NA, nrow = nrow(new_data_base), ncol = length(pop_year_levels))
colnames(log_preds_matrix) <- pop_year_levels

# Loop through each population-year, get predictions on log scale
for (i in seq_along(pop_year_levels)) {
  pop_year <- pop_year_levels[i]
  
  # Extract year (last 2 digits) and population (everything else before last "_")
  yr <- sub(".*_", "", pop_year)
  pop <- sub("_(\\d+)$", "", pop_year)
  
  new_data_temp <- new_data_base
  new_data_temp$population <- pop
  new_data_temp$year <- yr
  
  log_preds_matrix[, i] <- predict(model2, newdata = new_data_temp, type = "link")
}

# Average predictions on the log scale
new_data <- new_data_base
new_data$log_predicted_cortisol <- rowMeans(log_preds_matrix)

# Back-transform to response scale
new_data$predicted_cortisol <- exp(new_data$log_predicted_cortisol)

# Approximate SE across groups for CIs
log_se <- apply(log_preds_matrix, 1, sd)
new_data$lower_ci <- exp(new_data$log_predicted_cortisol - 1.96 * log_se)
new_data$upper_ci <- exp(new_data$log_predicted_cortisol + 1.96 * log_se)

# Make sure pop_year exists
cdata$pop_year <- interaction(cdata$population, cdata$year, sep = "_", drop = TRUE)
cdata$pop_year <- as.character(cdata$pop_year)  # ensure it's character, not factor
cdata$year <- as.character(cdata$year)

# plot
svglite("cortisol_vs_k-index_model2.svg", width = 8, height = 10)
ggplot(cdata, aes(x = k_index, y = cortisol)) +
  geom_point(aes(fill = pop_year, shape = year), color = "black", stroke = 1.2, size = 4, alpha = 0.9) +
  scale_fill_manual(values = c(
    "Rio_Choy_2023" = "#ff0000ff", "Rio_Choy_2024" = "#ff0000ff",
    "Presa_El_Oyul_2023" = "#6c006cff", "Presa_El_Oyul_2024" = "#6c006cff",
    "Pachón_2023" = "#131300ff", "Los_Sabinos_2023" = "#ffcc00ff",
    "Tinaja_2023" = "#ff2e8eff", "Tinaja_2024" = "#ff2e8eff",
    "Rio_Subterráneo_2023" = "#0000ffff", "Rio_Subterráneo_2024" = "#0000ffff"
  )) +
  scale_shape_manual(values = c("2023" = 21, "2024" = 24)) +  # 21 = circle, 24 = triangle
  geom_ribbon(data = new_data, aes(x = k_index, ymin = lower_ci, ymax = upper_ci), 
              fill = "blue", alpha = 0.1, inherit.aes = FALSE) +
  geom_line(data = new_data, aes(x = k_index, y = predicted_cortisol), 
            color = "blue", linewidth = 1) +
  labs(x = "Body Condition (k-index)", 
       y = "Cortisol (pg/mg scales)", 
       title = "Cortisol vs. Body Condition (Single Model-Based Trend Line)") +
  theme_minimal() +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20, face = "bold"),
        legend.text = element_text(size = 18), 
        legend.title = element_text(size = 20, face = "bold"))

dev.off()


## plot correlation between cortisol and handling time with model fitted trend line that is based on predictions for model 2

# Use the same dataset that was used to fit model2
plot_data <- cdata

# Create population-year groups for plotting
plot_data$pop_year <- interaction(
  plot_data$population,
  plot_data$year,
  sep = "_",
  drop = TRUE
)

# Mean used to centre handling time during model fitting
time_mean <- mean(
  plot_data$time_capture_dissection_hours,
  na.rm = TRUE
)

# Generate a sequence of actual handling times
time_seq <- seq(
  min(plot_data$time_capture_dissection_hours, na.rm = TRUE),
  max(plot_data$time_capture_dissection_hours, na.rm = TRUE),
  length.out = 100
)

# Convert to the centred scale used by model2
time_seq_c <- time_seq - time_mean

# Calculate model-based predictions
emm_time <- emmeans(
  model2,
  ~ time_capture_dissection_c,
  at = list(
    time_capture_dissection_c = time_seq_c,
    k_index = mean(plot_data$k_index, na.rm = TRUE),
    sex = "m"
  ),
  weights = "flat",
  type = "response"
)

new_data <- as.data.frame(emm_time)

# Extract predicted means and confidence limits
new_data$predicted_cortisol <- new_data$response

lower_col <- intersect(
  c("asymp.LCL", "lower.CL"), names(new_data)
)[1]

upper_col <- intersect(
  c("asymp.UCL", "upper.CL"), names(new_data)
)[1]

new_data$lower_ci <- new_data[[lower_col]]
new_data$upper_ci <- new_data[[upper_col]]

# Convert centred time back to actual hours
new_data$time_capture_dissection_hours <-
  new_data$time_capture_dissection_c + time_mean

# Plot
svglite("cortisol_vs_handling_time_model2.svg", width = 8, height = 10)
ggplot(plot_data, aes(x = time_capture_dissection_hours, y = cortisol)) +
  geom_ribbon(data = new_data, aes(x = time_capture_dissection_hours, ymin = lower_ci, ymax = upper_ci), fill = "blue", alpha = 0.1, inherit.aes = FALSE) +
  geom_line(data = new_data, aes(x = time_capture_dissection_hours, y = predicted_cortisol), color = "blue", linewidth = 1, inherit.aes = FALSE) +
  geom_point(aes(fill = pop_year, shape = factor(year)), color = "black", stroke = 1.2, size = 4, alpha = 0.9) +
  scale_fill_manual(values = c(
    "Rio_Choy_2023" = "#ff0000ff",
    "Rio_Choy_2024" = "#ff0000ff",
    "Presa_El_Oyul_2023" = "#6c006cff",
    "Presa_El_Oyul_2024" = "#6c006cff",
    "Pachón_2023" = "#131300ff",
    "Los_Sabinos_2023" = "#ffcc00ff",
    "Tinaja_2023" = "#ff2e8eff",
    "Tinaja_2024" = "#ff2e8eff",
    "Rio_Subterráneo_2023" = "#0000ffff",
    "Rio_Subterráneo_2024" = "#0000ffff"
  )) +
  scale_shape_manual(values = c(
    "2023" = 21,
    "2024" = 24
  )) +
  labs(
    x = "Time between capture and dissection (hours)",
    y = "Cortisol (pg/mg scales)",
    title = "Cortisol vs. Capture-to-Dissection Interval",
    fill = "Population-year",
    shape = "Year"
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 20, face = "bold")
  )
dev.off()


## calculate EMMs and plot sex differences for model 2
# Calculate EMMs for sex
emm_sex <- emmeans(model2, ~ sex, type = "response", weights = "flat",)  # type = "response" gives results on the original scale

# Convert to data frame for plotting
emmeans_table <- as.data.frame(emm_sex)

emmeans_table$exp_emmean <- emmeans_table$response
emmeans_table$exp_lower.CL <- emmeans_table$asymp.LCL
emmeans_table$exp_upper.CL <- emmeans_table$asymp.UCL

emmeans_table
# Write the table to a CSV file
write.csv(emmeans_table, "emmeans_table_sex_model2.csv", row.names = FALSE)


cdata_clean <- subset(cdata, !is.na(sex))

# plot
svglite("cortisol_sex_model2.svg", width = 8, height = 10)
ggplot() +
  geom_violin(data = cdata_clean, aes(x = sex, y = cortisol, fill = sex),
              trim = FALSE, alpha = 0.5, color = "black", width = 0.7) +
  geom_jitter(data = cdata_clean, aes(x = sex, y = cortisol, fill = sex),
              size = 1.7, alpha = 0.5, position = position_jitter(width = 0.13)) +
  geom_point(data = emmeans_table, aes(x = sex, y = exp_emmean),
             color = "black", size = 5, shape = 16, stroke = 2) +
  geom_errorbar(data = emmeans_table,
                aes(x = sex, ymin = lower.CL, ymax = upper.CL),
                width = 0, color = "black", linewidth = 2) +
  labs(title = "Cortisol Levels by Sex",
       x = "Sex", y = "Cortisol (pg/mg scales)") +
  scale_fill_manual(values = c("f" = "red", "m" = "lightblue")) +
  scale_y_continuous(breaks = seq(0, max(cdata$cortisol, na.rm = TRUE), by = 20)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 20, angle = 45, hjust = 1),
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 22),
    legend.position = "none",
    plot.title = element_text(size = 20, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90")
  )
dev.off()



### model 3 - checks for population x origin effects on cortisol levels 

## load and prepare data
# import cortisol data - cortisol_data_labvswild.csv (model3)
cdata2 <- read.csv("cortisol_data_labvswild.csv", header = TRUE, sep = ",", dec = ".")
head(cdata2)
cdata2

# check if cortisol is normally distributed, if yes - symmetrical bell shape without tail
# histogram
ggplot(cdata2, aes(x = cortisol)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black", alpha = 0.7) +
  geom_density(color = "red", size = 1) +
  labs(title = "Histogram of Cortisol Levels", x = "Cortisol (pg/mg)", y = "Density") +
  theme_minimal(
) # doesn't look like normal distribution!!!

# shapiro-wilk test, if p < 0.05 data is not normally distributed
shapiro.test(cdata2$cortisol) # in this case not normally distributed!!!

# Convert categorical variables to factors BEFORE modelling
cdata2$origin <- factor(cdata2$origin)
cdata2$sex <- factor(cdata2$sex)


## fit model 3
model3 <- glm(cortisol ~ population * origin + k_index + sex, data=cdata2, family = Gamma(link = "log"), na.action = na.exclude) 

# # alternative model 3, using standard_length instead of k_index as covariate
# model3 <- glm(cortisol ~ population * origin + standard_length + sex, data=cdata2, family = Gamma(link = "log"), na.action = na.exclude) 

# Summary of the model 
summary(model3)

# save output
capture.output(
  summary(model3),
  file = "summary_model3.txt"
)

# check if residuals are normally-distributed
# pull residuals from the model 
resid_standard <- residuals(model3)
resid_standard

# shapiro-wilk test, if p < 0.05 residuals are not normally distributed
shapiro.test(resid_standard) # residuals not normally distributed!!!

# test for pseudo R2 (1-(residual deviance/null deviance)) - variation explained
1-(14.476/20.103) # model3 = 0.28

# model diagnostics (type in model of interest)
svglite("model3.svg", width = 6, height = 10)
par(mfrow = c(2, 2)) # Set up plotting area for multiple plots
plot(model3)
dev.off()


## pairwise comparisons
# Calculate estimated marginal means for each population x origin combination
emmeans_table_3 <- emmeans(model3, ~ population * origin)

# pairwise comparisons between all combinations
pairwise_comparisons <- pairs(emmeans_table_3, adjust = "tukey")

# Extract the pairwise comparison table from interaction EMMs
pairwise_results <- as.data.frame(pairwise_comparisons)

# Add significance markers based on p-value thresholds
pairwise_results$significance <- ifelse(pairwise_results$p.value < 0.001, "***",
                                 ifelse(pairwise_results$p.value < 0.01, "**",
                                 ifelse(pairwise_results$p.value < 0.05, "*", "ns")))

# View the pairwise comparison results to check significance
head(pairwise_results)

# Write the results to a CSV file
write.csv(pairwise_results, "pairwise_comparisons_model3_significance.csv", row.names = FALSE)

# Get the compact letter display (CLD)
cld_results <- cld(emmeans_table_3, Letters = letters)

# View the CLD results
print(cld_results)
write.csv(cld_results, "clds_model3.csv", row.names = FALSE)

# Back-transform from log to original scale
emmeans_table_3 <- as.data.frame(emmeans_table_3)
emmeans_table_3$exp_emmean <- exp(emmeans_table_3$emmean)
emmeans_table_3$exp_lower.CL <- exp(emmeans_table_3$lower.CL)
emmeans_table_3$exp_upper.CL <- exp(emmeans_table_3$upper.CL)
emmeans_table_3$exp_SE <- exp(emmeans_table_3$SE)

emmeans_table_3
emmeans_table_3 <- emmeans_table_3 %>% filter(is.finite(exp_emmean))

# Create a unified column for plotting on the x-axis
emmeans_table_3$population_origin <- paste0(emmeans_table_3$population, "_", substr(emmeans_table_3$origin, 1, 4))
emmeans_table_3$population_origin

# Write the table to a CSV file
write.csv(emmeans_table_3, "emmeans_table_transformed_model3.csv", row.names = FALSE)

# Make sure cortisol is finite
cdata2 <- cdata2 %>% filter(is.finite(cortisol))

# Create matching population_year column in cdata
cdata2$population_origin <- paste0(cdata2$population, "_", substr(cdata2$origin, 1, 4))
cdata2$population_origin

population_order <- c("Rio_Choy_wild", "Rio_Choy_lab ", "Pachón_wild", "Pachón_lab ")

# Save SVG plot
svglite("cortisol_plot_model3_population_origin.svg", width = 8, height = 10)
ggplot() +
  geom_violin(data = cdata2, aes(x = factor(population_origin, levels=population_order), y = cortisol, fill = population_origin),
              trim = FALSE, alpha = 0.5, color = "black", width = 0.7) +
  geom_jitter(data = cdata2, aes(x = factor(population_origin, levels=population_order), y = cortisol, fill = population),
              size = 1.7, alpha = 0.5, position = position_jitter(width = 0.13)) +
  geom_point(data = emmeans_table_3, aes(x = population_origin, y = exp_emmean),
             color = "black", size = 4, shape = 16, stroke = 2) +
  geom_errorbar(data = emmeans_table_3,
                aes(x = population_origin, ymin = exp_lower.CL, ymax = exp_upper.CL),
                width = 0, color = "black", linewidth = 1.5) +
  labs(title = "Cortisol Levels by Population and Origin",
       x = "Population-Origin", y = "Cortisol (pg/mg scales)") +
  scale_fill_manual(values = c(
    "Rio_Choy_wild" = "#ff0000ff",
    "Rio_Choy_lab " = "#ff0000ff",
    "Pachón_wild" = "#131300ff",
    "Pachón_lab " = "#131300ff"
  )) +
  scale_y_continuous(breaks = seq(0, max(cdata2$cortisol, na.rm = TRUE) + 10, by = 20)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 16, angle = 45, hjust = 1),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.position = "none",
    plot.title = element_text(size = 20, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90")
  )
dev.off()


## plot correlation between cortisol and body condition with model fitted trend line that is based on predictions for model 2
# Clean population and origin
cdata2$population <- factor(trimws(as.character(cdata2$population)))
cdata2$origin <- factor(trimws(as.character(cdata2$origin)))

# Create combined factor: pop_origin
cdata2$pop_origin <- interaction(cdata2$population, cdata2$origin, sep = "_")
cdata2$pop_origin <- factor(trimws(as.character(cdata2$pop_origin)))

# Create shape by origin
cdata2$origin_shape <- ifelse(trimws(as.character(cdata2$origin)) == "lab", "lab", "wild")

# Define color map
my_colors_po <- c(
  "Rio_Choy_wild" = "#ff0000ff",
  "Rio_Choy_lab"  = "#ff0000ff",
  "Pachón_wild"   = "#131300ff",
  "Pachón_lab"    = "#131300ff"
)
names(my_colors_po) <- trimws(names(my_colors_po))

# Define shapes: wild = circle (21), lab = triangle (24)
shape_map <- c("wild" = 21, "lab" = 24)

# Prediction: Averaged line across all pop x origin 
# Create base data
k_seq <- seq(min(cdata2$k_index, na.rm = TRUE), max(cdata2$k_index, na.rm = TRUE), length.out = 100)
new_data_base <- data.frame(
  k_index = k_seq,
  sex = "m"
)

# Get levels
pop_levels <- levels(cdata2$population)
origin_levels <- levels(cdata2$origin)

# Matrix to store predictions
pred_matrix <- matrix(NA, nrow = length(k_seq), ncol = length(pop_levels) * length(origin_levels))
col_idx <- 1

for (pop in pop_levels) {
  for (org in origin_levels) {
    temp_data <- new_data_base
    temp_data$population <- factor(pop, levels = pop_levels)
    temp_data$origin <- factor(org, levels = origin_levels)
    
    pred_matrix[, col_idx] <- predict(model16, newdata = temp_data, type = "link")
    col_idx <- col_idx + 1
  }
}

# Average predictions and get CI
log_pred <- rowMeans(pred_matrix)
log_se <- apply(pred_matrix, 1, sd)

new_data_avg <- data.frame(
  k_index = k_seq,
  predicted = exp(log_pred),
  lower_ci = exp(log_pred - 1.96 * log_se),
  upper_ci = exp(log_pred + 1.96 * log_se)
)

# plot 
svglite("cortisol_vs_k-index_model3.svg", width = 8, height = 10)
ggplot(cdata2, aes(x = k_index, y = cortisol)) +
  # Data points
  geom_point(aes(fill = pop_origin, shape = origin_shape),
             color = "black", stroke = 1.2, size = 4, alpha = 0.9) +
  scale_fill_manual(values = my_colors_po) +
  scale_shape_manual(values = shape_map) +

  # Prediction ribbon + average line
  geom_ribbon(data = new_data_avg, aes(x = k_index, ymin = lower_ci, ymax = upper_ci),
              fill = "blue", alpha = 0.1, inherit.aes = FALSE) +
  geom_line(data = new_data_avg, aes(x = k_index, y = predicted),
            color = "blue", linewidth = 1.2, linetype = "solid", inherit.aes = FALSE) +

  # Labels and theme
  labs(
    x = "Body Condition (k-index)", 
    y = "Cortisol (pg/mg scales)", 
    title = "Cortisol vs. Body Condition (Model 3: Averaged Across Population × Origin)",
    fill = "Population × Origin",
    shape = "Origin"
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 18), 
    legend.title = element_text(size = 20, face = "bold"),
    plot.title = element_text(size = 20, face = "bold")
  )
dev.off()


## calculate EMMs and plot sex differences for model 3
emm_sex <- emmeans(model3, ~ sex, type = "response", weights = "flat")  # type = "response" gives results on the original scale

# Convert to data frame for plotting
emmeans_table <- as.data.frame(emm_sex)

emmeans_table$exp_emmean <- emmeans_table$response
emmeans_table$exp_lower.CL <- emmeans_table$asymp.LCL
emmeans_table$exp_upper.CL <- emmeans_table$asymp.UCL

emmeans_table
# Write the table to a CSV file
write.csv(emmeans_table, "emmeans_table_sex_model3.csv", row.names = FALSE)

cdata_clean <- subset(cdata2, !is.na(sex))

# plot
svglite("cortisol_sex_model3.svg", width = 8, height = 10)
ggplot() +
  geom_violin(data = cdata_clean, aes(x = sex, y = cortisol, fill = sex),
              trim = FALSE, alpha = 0.5, color = "black", width = 0.7) +
  geom_jitter(data = cdata_clean, aes(x = sex, y = cortisol, fill = sex),
              size = 1.7, alpha = 0.5, position = position_jitter(width = 0.13)) +
  geom_point(data = emmeans_table, aes(x = sex, y = exp_emmean),
             color = "black", size = 5, shape = 16, stroke = 2) +
  geom_errorbar(data = emmeans_table,
                aes(x = sex, ymin = lower.CL, ymax = upper.CL),
                width = 0, color = "black", linewidth = 2) +
  labs(title = "Cortisol Levels by Sex",
       x = "Sex", y = "Cortisol (pg/mg scales)") +
  scale_fill_manual(values = c("f" = "red", "m" = "lightblue")) +
  scale_y_continuous(breaks = seq(0, max(cdata$cortisol, na.rm = TRUE), by = 20)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 20, angle = 45, hjust = 1),
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 22),
    legend.position = "none",
    plot.title = element_text(size = 20, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90")
  )
dev.off()



### bootstrapping to compare individual variation among populations

# import cortisol data - cortisol_data.csv
cdata3 <- read.csv("cortisol_data.csv", header = TRUE, sep = ",", dec = ".")
head(cdata3)
cdata3

# Exclude missing group labels and invalid cortisol measurements
cdata3_clean <- cdata3 %>%
  filter(
    !is.na(population_year),
    is.finite(cortisol)
  )

# Function to calculate variance
var_func <- function(data, indices) {
  d <- data[indices]
  return(var(d))
}

# List to store results
boot_results <- list()

# Get unique populations
populations <- unique(cdata3_clean$population_year)

# Perform bootstrapping for each population
for (pop in populations) {
  pop_data <- cdata3_clean$cortisol[cdata3_clean$population_year == pop]
  boot_results[[pop]] <- boot(pop_data, var_func, R = 1000)
}

# Extract confidence intervals
ci_results <- lapply(boot_results, function(x) boot.ci(x, type = "perc"))


# Initialize an empty data frame
ci_table <- data.frame(Population = character(),
                       Lower_CI = numeric(),
                       Upper_CI = numeric(),
                       stringsAsFactors = FALSE)

# Store results in the data frame
for (pop in populations) {
  ci <- ci_results[[pop]]
  
  # Check if boot.ci was successful (avoid errors for small sample sizes)
  if (!is.null(ci) && !is.null(ci$perc)) {
    ci_table <- rbind(ci_table, 
                      data.frame(Population = pop,
                                 Lower_CI = ci$perc[4],  # Lower bound of 95% CI
                                 Upper_CI = ci$perc[5])) # Upper bound of 95% CI
  } else {
    ci_table <- rbind(ci_table, 
                      data.frame(Population = pop,
                                 Lower_CI = NA,
                                 Upper_CI = NA))
  }
}

# Print the table
print(ci_table)

# Save the table as a CSV file
write.csv(ci_table, "bootstrap_CI_results.csv", row.names = FALSE)


population_order <- c(  
  "Rio_Subterraneo_24", "Rio_Subterraneo_23", "Tinaja_24", "Tinaja_23", "Los_Sabinos_23", "Pachon_lab", "Pachon_23",
  "Presa_El_Oyul_24", "Presa_El_Oyul_23", "Rio_Choy_lab", "Rio_Choy_24", "Rio_Choy_23"
)

# Convert Population column to a factor with the specified order
ci_table$Population <- factor(ci_table$Population, levels = population_order)

# Create the plot
svglite("cortisol_variance.svg", width = 7, height = 10)
ggplot(ci_table, aes(x = Population, y = (Lower_CI + Upper_CI) / 2)) + 
  geom_point(size = 3) +  
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.2) +
  theme_minimal() +
  labs(title = "Bootstrapped Confidence Intervals of Variance",
       x = "Population",
       y = "Variance (with 95% CI)") +
  coord_flip() +  # Flip for better readability
  theme(
    axis.text = element_text(size = 14),  # Increase axis text size
    axis.title = element_text(size = 16),  # Increase axis titles size
    plot.title = element_text(size = 16, face = "bold")  # Increase title size and make it bold
  ) 
dev.off()



### exploratory analysis to test for effect of mean_abiotic_PCs on cortisol levels

# Prepare data and identify environmental sampling units
cdata_env <- cdata %>%
  mutate(
    environmental_unit = interaction(
      population, year, pool, drop = TRUE
    ),
    population = factor(population),
    year = factor(year),
    sex = factor(sex),
    time_capture_dissection_c =
      time_capture_dissection_hours -
      mean(time_capture_dissection_hours, na.rm = TRUE)
  ) %>%
  group_by(environmental_unit) %>%
  fill(mean_abiotic_PC1, mean_abiotic_PC2, .direction = "downup") %>%
  ungroup()

# Create datasets for each PC model
cdata_PC1 <- cdata_env %>%
  filter(
    is.finite(cortisol), cortisol > 0,
    is.finite(mean_abiotic_PC1),
    is.finite(k_index),
    is.finite(time_capture_dissection_c),
    !is.na(sex), !is.na(year),
    !is.na(population),
    !is.na(environmental_unit)
  ) %>%
  droplevels()

cdata_PC2 <- cdata_env %>%
  filter(
    is.finite(cortisol), cortisol > 0,
    is.finite(mean_abiotic_PC2),
    is.finite(k_index),
    is.finite(time_capture_dissection_c),
    !is.na(sex), !is.na(year),
    !is.na(population),
    !is.na(environmental_unit)
  ) %>%
  droplevels()

## fit models using these datasets
# model_PC1
model_PC1 <- glmmTMB(cortisol ~ mean_abiotic_PC1 + year + k_index + sex + time_capture_dissection_c + (1 | population) + (1 | environmental_unit),
  data = cdata_PC1,
  family = Gamma(link = "log")
)

summary(model_PC1)

# save output
capture.output(
  summary(model_PC1),
  file = "summary_model_PC1.txt"
)

svglite("model_PC1.svg", width = 8, height = 10)
res_model_PC1 <- simulateResiduals(model_PC1)
plot(res_model_PC1)
dev.off()

# model_PC2
model_PC2 <- glmmTMB(cortisol ~ mean_abiotic_PC2 + year + k_index + sex + time_capture_dissection_c + (1 | population) + (1 | environmental_unit),
  data = cdata_PC2,
  family = Gamma(link = "log")
)

summary(model_PC2)

# save output
capture.output(
  summary(model_PC2),
  file = "summary_model_PC2.txt"
)

svglite("model_PC2.svg", width = 8, height = 10)
res_model_PC2 <- simulateResiduals(model_PC2)
plot(res_model_PC2)
dev.off()

# Helper function for model-based prediction curves
make_curve <- function(model, dat, pc_var, pc_label) {
  x_seq <- seq(
    min(dat[[pc_var]], na.rm = TRUE),
    max(dat[[pc_var]], na.rm = TRUE),
    length.out = 100
  )
  # Predict across PC values, holding continuous covariates
  # at their means and averaging equally over year and sex
  at_values <- list(
    k_index = mean(dat$k_index, na.rm = TRUE),
    time_capture_dissection_c =
      mean(dat$time_capture_dissection_c, na.rm = TRUE)
  )
  at_values[[pc_var]] <- x_seq
  pred <- as.data.frame(
    emmeans(
      model,
      specs = as.formula(paste("~", pc_var)),
      at = at_values,
      weights = "equal",
      type = "response"
    )
  )
  # Identify estimated means and confidence limits
  fit_col <- if ("response" %in% names(pred)) {
    "response"
  } else {
    "emmean"
  }
  lower_col <- intersect(
    c("asymp.LCL", "lower.CL"), names(pred)
  )[1]
  upper_col <- intersect(
    c("asymp.UCL", "upper.CL"), names(pred)
  )[1]
  data.frame(
    axis = pc_label,
    x = pred[[pc_var]],
    fit = pred[[fit_col]],
    lower = pred[[lower_col]],
    upper = pred[[upper_col]]
  )
}

# Use the same filtered datasets used to fit each model
# cdata_PC1 and cdata_PC2 should contain the respective complete-case observations from your model workflow.

predictions <- bind_rows(
  make_curve(
    model_PC1, cdata_PC1,
    "mean_abiotic_PC1", "PC1"
  ),
  make_curve(
    model_PC2, cdata_PC2,
    "mean_abiotic_PC2", "PC2"
  )
)

# Summarise cortisol for each environmental sampling unit
plot_units <- cdata_env %>%
  group_by(
    environmental_unit, population, year, pool
  ) %>%
  summarise(
    n_fish = sum(is.finite(cortisol) & cortisol > 0),
    mean_cortisol = mean(
      cortisol[cortisol > 0 & is.finite(cortisol)],
      na.rm = TRUE
    ),
    sd_cortisol = sd(
      cortisol[cortisol > 0 & is.finite(cortisol)],
      na.rm = TRUE
    ),
    PC1 = mean(mean_abiotic_PC1, na.rm = TRUE),
    PC2 = mean(mean_abiotic_PC2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    se_cortisol = sd_cortisol / sqrt(n_fish),
    ci_halfwidth = ifelse(
      n_fish > 1,
      qt(0.975, df = n_fish - 1) * se_cortisol,
      NA_real_
    ),
    lower_mean = mean_cortisol - ci_halfwidth,
    upper_mean = mean_cortisol + ci_halfwidth
  ) %>%
  pivot_longer(
    cols = c(PC1, PC2),
    names_to = "axis",
    values_to = "x"
  ) %>%
  filter(
    is.finite(mean_cortisol),
    is.finite(x)
  ) %>%
  mutate(
    axis = factor(axis, levels = c("PC1", "PC2")),
    year = factor(year)
  )

p_abiotic <- ggplot() +
  geom_ribbon(data = predictions, aes(x = x, ymin = lower, ymax = upper), fill = "grey70", alpha = 0.4) + # Confidence intervals from the fitted GLMMs
  geom_line(data = predictions, aes(x = x, y = fit), color = "black", linewidth = 1) + # Model-adjusted trend lines
  geom_errorbar(data = plot_units, aes(x = x, ymin = lower_mean, ymax = upper_mean, color = population), width = 0) + # Observed means for each environmental sampling unit
  geom_point(data = plot_units, aes(x = x, y = mean_cortisol, color = population, shape = year), size = 4.5) +
  scale_color_manual(values = c("Rio_Choy" = "#ff0000ff", "Presa_El_Oyul" = "#6c006cff", "Pachón" = "#131300ff", "Los_Sabinos" = "#ffcc00ff", "Tinaja" = "#ff2e8eff", "Rio_Subterráneo" = "#0000ffff")) +
  facet_wrap(~ axis, scales = "free_x") +
  labs(x = "Environmental principal component score", y = "Mean scale cortisol (pg/mg)", color = "Population", shape = "Year") +
  theme_classic(base_size = 18)
p_abiotic

ggsave(
  "supplementary_abiotic_PC_cortisol.svg",
  p_abiotic,
  width = 9,
  height = 10
)


