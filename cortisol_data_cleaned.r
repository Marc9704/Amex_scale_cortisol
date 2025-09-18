
### cortisol data

library(ggplot2)
library(dplyr)
library(svglite)

# set working directory

setwd("")

# import cortisol data - cortisol_data_wild.csv (model1, 2), cortisol_data_labvswild.csv (model3)

cdata <- read.csv("cortisol_data_wild.csv", header = TRUE, sep = ",", dec = ".")

head(cdata)
cdata


### statistical analysis

# check if cortisol is normally distributed, if yes - symmetrical bell shape without tail

# histogram
ggplot(cdata, aes(x = cortisol)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black", alpha = 0.7) +
  geom_density(color = "red", size = 1) +
  labs(title = "Histogram of Cortisol Levels", x = "Cortisol (pg/mg)", y = "Density") +
  theme_minimal() # doesn't look like normal distribution!!!

# shapiro-wilk test, if p < 0.05 data is not normally distributed

shapiro.test(cdata$cortisol) # in this case not normally distributed!!!


library(emmeans)  # For pairwise comparisons
library(car)
library(multcomp)

# Fit the generalized-linear models

# model 1 - checks for overall habitat effect on cortisol levels
model1 <- glm(cortisol ~ habitat + k_index + sex, data=cdata, family = Gamma(link = "log"), na.action = na.exclude) 

# set population and year as categorical factor
cdata$year <- as.factor(cdata$year)
cdata$population <- as.factor(cdata$population)

#check if it worked
levels(cdata$year)
levels(cdata$population)

table(cdata$population, cdata$year)

# model 2 - checks for population x year effects on cortisol levels
model2 <- glm(cortisol ~ population * year + k_index + sex, data=cdata, family = Gamma(link = "log"), na.action = na.exclude) # second model used! pseudo R2=0.4590323, AIC=1293.5


cdata$origin

# model 3 - checks for population x origin effects on cortisol levels 
model3 <- glm(cortisol ~ population * origin + k_index + sex, data=cdata, family = Gamma(link = "log"), na.action = na.exclude) # third model used! pseudo R2=0.2799085, AIC=488.14



## additional models

# model4 <- glm(cortisol ~ population, data = cdata, family = Gamma(link = "log")) 

# model5 <- glm(cortisol ~ population + k_index, data=cdata, family = Gamma(link = "log")) 

# model6 <- glm(cortisol ~ k_index * population, data=cdata, family = Gamma(link = "log")) 

# model7 <- glm(cortisol ~ population + sex, data=cdata, family = Gamma(link = "log")) 

# model8 <- glm(cortisol ~ population + k_index + sex, data=cdata, family = Gamma(link = "log")) 

# model9 <- glm(cortisol ~ population + mass, data=cdata, family = Gamma(link = "log"))

# model10 <- glm(cortisol ~ population * mass, data=cdata, family = Gamma(link = "log"))

# model11 <- glm(cortisol ~ population_year + k_index + standard_length_c + sex, data=cdata, family = Gamma(link = "log"), na.action = na.exclude) 

# model12 <- glm(cortisol ~ population + k_index + standard_length_c, data=cdata, family = Gamma(link = "log"), na.action = na.exclude) 

# model13 <- glm(cortisol ~ population_year + mass + standard_length_c + sex, data=cdata, family = Gamma(link = "log"), na.action = na.exclude)

# model14 <- glm(cortisol ~ population_year + k_index + habitat + standard_length_c + sex, data=cdata, family = Gamma(link = "log"), na.action = na.exclude)

# model15 <- glm(cortisol ~ habitat + mass + sex, data=cdata, family = Gamma(link = "log"), na.action = na.exclude) 

# model16 <- glm(cortisol ~ habitat * k_index + sex, data=cdata, family = Gamma(link = "log"), na.action = na.exclude) 

# model17 <- glm(cortisol ~ habitat * sex + k_index, data=cdata, family = Gamma(link = "log"), na.action = na.exclude) 

# model18 <- glm(cortisol ~ population * year * k_index + sex, data=cdata, family = Gamma(link = "log"), na.action = na.exclude) 

# model19 <- glm(cortisol ~ population * year * sex + k_index, data=cdata, family = Gamma(link = "log"), na.action = na.exclude) 

# model20 <- glm(cortisol ~ population * origin * k_index + sex, data=cdata, family = Gamma(link = "log"), na.action = na.exclude) 

# model21 <- glm(cortisol ~ habitat/population * year + k_index + sex, data = cdata, family = Gamma(link = "log"), na.action = na.exclude) 

# model22 <- glm(cortisol ~ population_year + k_index + sex, data=cdata, family = Gamma(link = "log"), na.action = na.exclude) 


# check if residuals are normally-distributed

# pull residuals from the model (type in model of interest)
resid_standard <- residuals(model1)

resid_standard

# shapiro-wilk test, if p < 0.05 residuals are not normally distributed

shapiro.test(resid_standard) # residuals not normally distributed!!!


# Summary of the model (type in model of interest)
summary(model1)


# test for pseudo R2 (1-(residual deviance/null deviance)) - variation explained
1-(36.148/66.821)


# model diagnostics (type in model of interest)

svglite("model1.svg", width = 6, height = 10)

par(mfrow = c(2, 2)) # Set up plotting area for multiple plots
plot(model1)

dev.off()

### pairwise comparisons - type in model of interest

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
write.csv(pairwise_results, "pairwise_comparisons_model15_significance.csv", row.names = FALSE)


### bootstrapping to compare individual variation among populations

library(boot)

# Function to calculate variance
var_func <- function(data, indices) {
  d <- data[indices]
  return(var(d))
}

# List to store results
boot_results <- list()

# Get unique populations
populations <- unique(cdata$population)

# Perform bootstrapping for each population
for (pop in populations) {
  pop_data <- cdata$cortisol[cdata$population == pop]
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



### calculate EMMs and plot data for model 1 

# 1. Calculate EMMs for habitat
emmeans_table <- emmeans(model1, ~ habitat)

# 2. Back-transform emmeans from log to original scale
emmeans_table <- as.data.frame(emmeans_table)  # Convert to data frame
emmeans_table$exp_emmean <- exp(emmeans_table$emmean)
emmeans_table$exp_lower.CL <- exp(emmeans_table$lower.CL)
emmeans_table$exp_upper.CL <- exp(emmeans_table$upper.CL)
emmeans_table$exp_SE <- exp(emmeans_table$SE)


# Optional: Rename habitat column for consistent plotting
names(emmeans_table)[names(emmeans_table) == "habitat"] <- "population"

# 3. Write EMMs to CSV
write.csv(emmeans_table, "emmeans_table_transformed_model1.csv", row.names = FALSE)

# 4. Filter NAs in cortisol (if not done already)
cdata <- cdata %>% filter(is.finite(cortisol))

# Set the order of habitat levels
cdata$habitat <- factor(cdata$habitat, levels = c("surface", "cave"))
emmeans_table$population <- factor(emmeans_table$population, levels = c("surface", "cave"))


# 5. Violin plot of raw data + back-transformed EMMs
svglite("cortisol_data_raw_back-transformed-emm_model1.svg", width = 7, height = 10)

ggplot() +
  geom_violin(data = cdata, aes(x = habitat, y = cortisol, fill = habitat), 
              trim = FALSE, alpha = 0.5, color = "black", width = 1.0) +
  geom_jitter(data = cdata, aes(x = habitat, y = cortisol, fill = habitat), 
              size = 1.7, alpha = 0.5, position = position_jitter(width = 0.13)) +
  geom_point(data = emmeans_table, aes(x = population, y = exp_emmean), 
             color = "black", size = 5, shape = 16, stroke = 2) +
  geom_errorbar(data = emmeans_table, 
                aes(x = population, ymin = exp_lower.CL, ymax = exp_upper.CL), 
                width = 0, color = "black", linewidth = 2) +
  labs(title = "Cortisol Levels by Habitat",
       x = "Habitat", y = "Cortisol (pg/mg scales)") +
  scale_fill_manual(values = c("cave" = "gray", "surface" = "white")) +
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



### compute and back-transform EMMs by population * year, model 2, plot data


# Calculate estimated marginal means for each population × year combination
emmeans_table <- emmeans(model2, ~ population * year)

# Back-transform from log to original scale
emmeans_table <- as.data.frame(emmeans_table)
emmeans_table$exp_emmean <- exp(emmeans_table$emmean)
emmeans_table$exp_lower.CL <- exp(emmeans_table$lower.CL)
emmeans_table$exp_upper.CL <- exp(emmeans_table$upper.CL)
emmeans_table$exp_SE <- exp(emmeans_table$SE)

emmeans_table
emmeans_table <- emmeans_table %>% filter(is.finite(exp_emmean))

# Write the table to a CSV file
write.csv(emmeans_table, "emmeans_table_transformed_model2.csv", row.names = FALSE)


# Create a unified column for plotting on the x-axis
emmeans_table$population_year <- paste0(emmeans_table$population, "_", substr(emmeans_table$year, 3, 4))

# Make sure cortisol is finite
cdata <- cdata %>% filter(is.finite(cortisol))

# Create matching population_year column in cdata
cdata$population_year <- paste0(cdata$population, "_", substr(cdata$year, 3, 4))

population_order <- c(  
  "Rio_Choy_23", "Rio_Choy_24", "Presa_El_Oyul_23", "Presa_El_Oyul_24", "Pachón_23", "Los_Sabinos_23", "Tinaja_23", "Tinaja_24", "Rio_Subterráneo_23", "Rio_Subterráneo_24" 
)

# Save SVG output
svglite("cortisol_data_raw_back-transformed-emm_model2.svg", width = 8, height = 10)

ggplot() +
  geom_violin(data = cdata, aes(x = factor(population_year, levels = population_order), y = cortisol, fill = population_year),
              trim = FALSE, alpha = 0.5, color = "black", width = 1.0) +
  geom_jitter(data = cdata, aes(x = factor(population_year, levels = population_order), y = cortisol, fill = population),
              size = 1.7, alpha = 0.5, position = position_jitter(width = 0.13)) +
  geom_point(data = emmeans_table, aes(x = population_year, y = exp_emmean),
             color = "black", size = 4, shape = 16, stroke = 2) +
  geom_errorbar(data = emmeans_table,
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



### model 3

# Calculate estimated marginal means for each population x origin combination
emmeans_table <- emmeans(model3, ~ population * origin)

# pairwise comparisons between all combinations
pairwise_comparisons <- pairs(emmeans_table, adjust = "tukey")

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
cld_results <- cld(emmeans_table, Letters = letters)

# View the CLD results
print(cld_results)

write.csv(cld_results, "clds_model3.csv", row.names = FALSE)


# Back-transform from log to original scale
emmeans_table <- as.data.frame(emmeans_table)
emmeans_table$exp_emmean <- exp(emmeans_table$emmean)
emmeans_table$exp_lower.CL <- exp(emmeans_table$lower.CL)
emmeans_table$exp_upper.CL <- exp(emmeans_table$upper.CL)
emmeans_table$exp_SE <- exp(emmeans_table$SE)

emmeans_table
emmeans_table <- emmeans_table %>% filter(is.finite(exp_emmean))

# Create a unified column for plotting on the x-axis
emmeans_table$population_origin <- paste0(emmeans_table$population, "_", substr(emmeans_table$origin, 1, 4))
emmeans_table$population_origin

# Write the table to a CSV file
write.csv(emmeans_table, "emmeans_table_transformed_model16.csv", row.names = FALSE)

# Make sure cortisol is finite
cdata <- cdata %>% filter(is.finite(cortisol))

# Create matching population_year column in cdata
cdata$population_origin <- paste0(cdata$population, "_", substr(cdata$origin, 1, 4))
cdata$population_origin

population_order <- c("Rio_Choy_wild", "Rio_Choy_lab ", "Pachón_wild", "Pachón_lab ")


# Save SVG plot
svglite("cortisol_plot_model3_population_origin.svg", width = 8, height = 10)

ggplot() +
  geom_violin(data = cdata, aes(x = factor(population_origin, levels=population_order), y = cortisol, fill = population_origin),
              trim = FALSE, alpha = 0.5, color = "black", width = 1.0) +
  geom_jitter(data = cdata, aes(x = factor(population_origin, levels=population_order), y = cortisol, fill = population),
              size = 1.7, alpha = 0.5, position = position_jitter(width = 0.13)) +
  geom_point(data = emmeans_table, aes(x = population_origin, y = exp_emmean),
             color = "black", size = 4, shape = 16, stroke = 2) +
  geom_errorbar(data = emmeans_table,
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
  scale_y_continuous(breaks = seq(0, max(cdata$cortisol, na.rm = TRUE) + 10, by = 20)) +
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





### plot correlation between cortisol and body condition with model fitted trend line that is based on predictions for model 1


# Create a new dataset for prediction over k_index, averaging over population
new_data <- data.frame(
  k_index = seq(min(cdata$k_index, na.rm = TRUE), 
                max(cdata$k_index, na.rm = TRUE), 
                length.out = 100),
  sex = "m"
)

# Get model matrix for population (to compute average effect)
habitat_levels <- unique(cdata$habitat)
hab_effects <- sapply(habitat_levels, function(hab) {
  new_data$habitat <- hab
  predict(model13, newdata = new_data, type = "link")  # Get log-scale predictions
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
  scale_fill_manual(values = c("cave" = "gray", "surface" = "white")) +
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


# same for model 2

cdata$pop_year <- interaction(cdata$population, cdata$year, sep = "_")

# Define the range of k_index
new_data_base <- data.frame(
  k_index = seq(min(cdata$k_index, na.rm = TRUE),
                max(cdata$k_index, na.rm = TRUE),
                length.out = 100),
  sex = "m"  # Set to a reference level or average as needed
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
  
  log_preds_matrix[, i] <- predict(model14, newdata = new_data_temp, type = "link")
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

setdiff(unique(cdata$pop_year), names(my_colors))

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



# same for model 3

# Clean population and origin
cdata$population <- factor(trimws(as.character(cdata$population)))
cdata$origin <- factor(trimws(as.character(cdata$origin)))

# Create combined factor: pop_origin
cdata$pop_origin <- interaction(cdata$population, cdata$origin, sep = "_")
cdata$pop_origin <- factor(trimws(as.character(cdata$pop_origin)))

# Create shape by origin
cdata$origin_shape <- ifelse(trimws(as.character(cdata$origin)) == "lab", "lab", "wild")

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
k_seq <- seq(min(cdata$k_index, na.rm = TRUE), max(cdata$k_index, na.rm = TRUE), length.out = 100)
new_data_base <- data.frame(
  k_index = k_seq,
  sex = "m"
)

# Get levels
pop_levels <- levels(cdata$population)
origin_levels <- levels(cdata$origin)

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

# Final Plot 

svglite("cortisol_vs_k-index_model3.svg", width = 8, height = 10)

ggplot(cdata, aes(x = k_index, y = cortisol)) +
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
    title = "Cortisol vs. Body Condition (Model 16: Averaged Across Population × Origin)",
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



### calculate EMMs and plot sex differences for each model

# model 1

# Calculate EMMs for sex
emm_sex <- emmeans(model1, ~ sex, type = "response")  # type = "response" gives results on the original scale

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


# model 2

# Calculate EMMs for sex
emm_sex <- emmeans(model2, ~ sex, type = "response")  # type = "response" gives results on the original scale

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
              trim = FALSE, alpha = 0.5, color = "black", width = 1.0) +
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


# model 3

# Calculate EMMs for sex
emm_sex <- emmeans(model3, ~ sex, type = "response")  # type = "response" gives results on the original scale

# Convert to data frame for plotting
emmeans_table <- as.data.frame(emm_sex)

emmeans_table$exp_emmean <- emmeans_table$response
emmeans_table$exp_lower.CL <- emmeans_table$asymp.LCL
emmeans_table$exp_upper.CL <- emmeans_table$asymp.UCL

emmeans_table
# Write the table to a CSV file
write.csv(emmeans_table, "emmeans_table_sex_model3.csv", row.names = FALSE)


cdata_clean <- subset(cdata, !is.na(sex))

# plot

svglite("cortisol_sex_model3.svg", width = 8, height = 10)

ggplot() +
  geom_violin(data = cdata_clean, aes(x = sex, y = cortisol, fill = sex),
              trim = FALSE, alpha = 0.5, color = "black", width = 1.0) +
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



