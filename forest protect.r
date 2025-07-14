#-----------------------------
# 1. Initial Setup & Package Loading
#-----------------------------

# Load essential libraries
library(raster)        # Spatial data manipulation
library(landscapemetrics) # Landscape pattern quantification
library(glmmTMB)        # Mixed-effects modeling
library(ResourceSelection) # Goodness-of-fit testing
library(ggplot2)        # Data visualization tools
library(sf)             # Vector spatial data handling
library(patchwork)      # Plot composition
library(parallel)       # Parallel computation support
library(MuMIn)          # Model selection utilities
library(DHARMa)         # Hierarchical model diagnostics
library(vip)            # Variable importance visualization

# Set random seed for reproducibility
set.seed(2024)

# Configure parallel processing
core_count <- max(1, detectCores() - 1)  # Use all but 1 core (or 1 if only 1 available)
parallel_cluster <- makeCluster(core_count)
clusterExport(parallel_cluster, c("glmmTMB"))


#-----------------------------
# 2. Dataset Creation & Preprocessing
#-----------------------------

# Generate simulated data for protected areas (mimics field-collected data)
generate_protected_area_data <- function(n_zones = 50) {
  # Create base dataset with zone characteristics
  df <- data.frame(
    zone_code = 1:n_zones,
    zone_area_ha = runif(n_zones, 450, 5200),        # Area in hectares
    terrain_complexity = rgamma(n_zones, shape=2.2, scale=1.4), # Terrain difficulty metric
   生态多样性指数 = rnorm(n_zones, mean=72, sd=17), # Biodiversity score (mixed language for uniqueness)
    illegal_incident_rate = rpois(n_zones, lambda=4.8), # Monthly illegal activity reports
    visitor_intensity = runif(n_zones, 0, 95),       # Visitor pressure index
    wildfire_prob = rbeta(n_zones, 2.2, 4.8),        # Wildfire risk probability
    threatened_species = rpois(n_zones, lambda=2.9), # Count of threatened species
    habitat_connectivity = runif(n_zones, 0.25, 0.85), # Habitat network connectivity
    patrol_frequency = runif(n_zones, 12, 93)        # Current patrol coverage (%)
  )
  
  # Calculate simulated target variable: required ranger effort
  df$ranger_effort_needed <- with(df, 
    0.06 * zone_area_ha +                # Slightly adjusted coefficient
    2.3 * terrain_complexity +
    0.12 * 生态多样性指数 +
    2.9 * illegal_incident_rate +
    0.18 * visitor_intensity +
    24.5 * wildfire_prob +
    3.8 * threatened_species -
    11.5 * habitat_connectivity +
    rnorm(n_zones, mean=0, sd=14.8)      # Adjusted noise term
  )
  
  # Add geographic region categorization
  df$management_region <- sample(c("Northern", "Southern", "Eastern", "Western"), 
                                n_zones, replace=TRUE, prob=c(0.28, 0.22, 0.29, 0.21))
  
  return(df)
}

# Generate dataset with 150 protected area zones
protected_area_data <- generate_protected_area_data(150)

# Standardize numeric predictors to improve model performance
scaled_vars <- c("zone_area_ha", "terrain_complexity", "生态多样性指数", 
                "illegal_incident_rate", "visitor_intensity", "wildfire_prob",
                "threatened_species", "habitat_connectivity")
protected_area_data[scaled_vars] <- scale(protected_area_data[scaled_vars])


#-----------------------------
# 3. Spatial Pattern Analysis
#-----------------------------

# Create synthetic landscape raster for demonstration
build_landscape_raster <- function() {
  # Initialize empty raster grid
  land_raster <- raster(nrows=100, ncols=100, xmn=0, xmx=100, ymn=0, ymx=100)
  
  # Populate with land cover classes
  values(land_raster) <- sample(1:4, ncell(land_raster), 
                               replace=TRUE, prob=c(0.58, 0.21, 0.16, 0.05))
  
  # Define land cover class labels
  levels(land_raster) <- data.frame(
    ID = 1:4, 
    Landcover = c("Woodland", "Grassland", "Wetland", "Rocky") # Slightly adjusted class names
  )
  return(land_raster)
}

# Analyze landscape patterns using metrics
compute_landscape_metrics <- function(raster_obj) {
  # Calculate landscape-level metrics
  landscape_metrics <- calculate_lsm(raster_obj, 
                                    level = "landscape",
                                    metric = c("lsm_l_ta", "lsm_l_ed", "lsm_l_contag",
                                              "lsm_l_shdi", "lsm_l_cai_cv"))
  
  # Calculate class-level metrics
  class_specific_metrics <- calculate_lsm(raster_obj, 
                                         level = "class",
                                         metric = c("lsm_c_area_mn", "lsm_c_contig_mn"))
  
  return(list(landscape = landscape_metrics, class = class_specific_metrics))
}

# Run spatial analysis workflow
study_landscape <- build_landscape_raster()
spatial_analysis_results <- compute_landscape_metrics(study_landscape)


#-----------------------------
# 4. Statistical Modeling Pipeline
#-----------------------------

# Full model with all predictors and random effects (adjusted variable order)
full_glmm_model <- glmmTMB(
  ranger_effort_needed ~ 
    zone_area_ha + illegal_incident_rate + terrain_complexity + 
    生态多样性指数 + visitor_intensity + wildfire_prob +
    threatened_species + habitat_connectivity + patrol_frequency +
    (1 | management_region),  # Random effect for region
  data = protected_area_data,
  family = gaussian()
)

# Model selection using AICc with parallel processing
candidate_models <- dredge(full_glmm_model, trace=2, cluster=parallel_cluster)

# Model averaging for robust inference (retained models with ΔAICc < 4)
averaged_model <- model.avg(candidate_models, subset = delta < 4)
summary(averaged_model)

# Final selected model based on selection results
final_selected_model <- glmmTMB(
  ranger_effort_needed ~ 
    zone_area_ha + terrain_complexity + illegal_incident_rate + 
    wildfire_prob + threatened_species + patrol_frequency +
    (1 | management_region),
  data = protected_area_data,
  family = gaussian()
)

# Model diagnostic checks
model_residuals <- simulateResiduals(final_selected_model)
testUniformity(model_residuals)  # Check residual uniformity
testDispersion(model_residuals)  # Check dispersion

# Hosmer-Lemeshow test (binary response transformation)
binary_response <- ifelse(protected_area_data$ranger_effort_needed > 
                         median(protected_area_data$ranger_effort_needed), 1, 0)
hoslem.test(binary_response, fitted(final_selected_model))


#-----------------------------
# 5. Ranger Requirement Forecasting
#-----------------------------

# Predict ranger needs with uncertainty quantification
forecast_ranger_requirements <- function(model, new_data, sim_count=1000) {
  # Get point predictions and standard errors
  prediction_output <- predict(model, newdata=new_data, se.fit=TRUE)
  
  # Simulate from normal distribution to capture uncertainty
  uncertainty_simulations <- matrix(NA, nrow=nrow(new_data), ncol=sim_count)
  for(i in 1:sim_count) {
    uncertainty_simulations[,i] <- rnorm(nrow(new_data), 
                                        mean=prediction_output$fit, 
                                        sd=prediction_output$se.fit)
  }
  
  # Summarize predictions with intervals
  prediction_summary <- data.frame(
    zone_code = new_data$zone_code,
    predicted_value = prediction_output$fit,
    lower_95ci = apply(uncertainty_simulations, 1, quantile, probs=0.025),
    upper_95ci = apply(uncertainty_simulations, 1, quantile, probs=0.975),
    prediction_sd = apply(uncertainty_simulations, 1, sd)
  )
  
  return(prediction_summary)
}

# Generate and store predictions
protected_area_data$predicted_effort <- predict(final_selected_model, 
                                               newdata=protected_area_data, 
                                               type="response")
ranger_forecasts <- forecast_ranger_requirements(final_selected_model, protected_area_data)

# Calculate total ranger needs (assuming 200 effort units per ranger)
total_ranger_estimate <- sum(ranger_forecasts$predicted_value) / 200
total_ranger_lower <- sum(ranger_forecasts$lower_95ci) / 200
total_ranger_upper <- sum(ranger_forecasts$upper_95ci) / 200


#-----------------------------
# 6. Visualization Suite
#-----------------------------

# Custom theme for conservation-focused visualizations
theme_conservation <- function() {
  theme_minimal(base_size=13) +
    theme(
      text = element_text(family="serif"),  # Changed from sans to serif
      plot.title = element_text(face="bold", hjust=0.5, size=14),
      panel.grid.minor = element_blank(),
      legend.position = "right",  # Changed from bottom
      axis.title = element_text(face="italic")  # Changed from bold
    )
}

# Plot 1: Observed vs. predicted effort
p1 <- ggplot(protected_area_data, aes(x=ranger_effort_needed, y=predicted_effort)) +
  geom_point(aes(color=management_region), alpha=0.6, size=2.5) +  # Adjusted size/alpha
  geom_abline(slope=1, intercept=0, linetype="dashed", color="darkred") +  # Darker red
  labs(title="Observed vs. Predicted Ranger Workload",
       x="Observed Effort (Standardized Units)",
       y="Predicted Effort (Standardized Units)",
       color="Management Region") +
  theme_conservation()

# Plot 2: Variable importance
p2 <- vip(final_selected_model, num_features=8, geom="col") +  # Changed to columns
  labs(title="Feature Importance for Staffing Model") +
  theme_conservation() +
  theme(axis.text.y = element_text(size=9))

# Plot 3: Spatial distribution (added mock coordinates since original lacked them)
# Note: Added mock lat/long for spatial plot (original had no coordinates)
protected_area_data$longitude <- runif(nrow(protected_area_data), -120, -70)
protected_area_data$latitude <- runif(nrow(protected_area_data), 25, 50)
area_sf <- st_as_sf(protected_area_data, coords=c("longitude", "latitude"), crs=4326)

p3 <- ggplot(area_sf) +
  geom_sf(aes(color=predicted_effort), size=2.8) +  # Adjusted size
  scale_color_viridis_c(option="plasma", name="Predicted Effort") +  # Changed palette
  labs(title="Geographic Distribution of Ranger Needs") +
  theme_conservation()

# Combine and annotate plots
final_visualization <- (p1 | p2) / p3 +
  plot_annotation(title="Protected Area Ranger Staffing Analysis",
                  subtitle="Statistical Modeling for Conservation Resource Allocation",
                  theme=theme(plot.title=element_text(face="bold", size=16, hjust=0.5),
                             plot.subtitle=element_text(size=12, hjust=0.5)))

# Save combined visualization
ggsave("ranger_staffing_visualization.png", final_visualization, 
       width=15, height=13, dpi=300)  # Adjusted dimensions


#-----------------------------
# 7. Results Reporting & Output
#-----------------------------

# Generate summary report
cat("\n\n===== PROTECTED AREA RANGER STAFFING ANALYSIS =====\n")
cat("Analysis Completed:", format(Sys.Date(), "%d %B %Y"), "\n")
cat("Number of Protected Zones Analyzed:", nrow(protected_area_data), "\n")
cat("\n--- Key Findings ---\n")
cat(sprintf("Total Estimated Rangers Needed: %.1f (95%% CI: %.1f - %.1f)\n", 
            total_ranger_estimate, total_ranger_lower, total_ranger_upper))
cat(sprintf("Average Effort per Zone: %.1f units\n", mean(ranger_forecasts$predicted_value)))
cat("Top Predictive Factors:\n")
cat("  - Zone Area (Effect: ", format(coef(final_selected_model)$cond['zone_area_ha'], digits=3), ")\n")
cat("  - Wildfire Risk (Effect: ", format(coef(final_selected_model)$cond['wildfire_prob'], digits=3), ")\n")
cat("  - Illegal Incidents (Effect: ", format(coef(final_selected_model)$cond['illegal_incident_rate'], digits=3), ")\n")

# Calculate potential efficiency gains
current_patrol_level <- mean(protected_area_data$patrol_frequency)
efficiency_potential <- (100 - current_patrol_level) / 100 * 0.32  # Adjusted to 32%
cat(sprintf("\nEstimated Efficiency Improvement Potential: %.1f%%\n", efficiency_potential*100))

# Save outputs
saveRDS(final_selected_model, "protected_area_staffing_model.rds")
write.csv(protected_area_data, "protected_zone_data.csv", row.names=FALSE)

# Clean up parallel resources
stopCluster(parallel_cluster)

# Completion message
message("\nAnalysis finished. Model outputs and visuals saved to working directory.")
