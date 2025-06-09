# Install required packages (only once)
# install.packages(c("randomForest", "caret", "ggplot2", "mgcv", "Metrics", "tidyverse"))

# Load libraries
library(randomForest)
library(caret)
library(ggplot2)
library(mgcv)
library(Metrics)
library(tidyverse)

# 🔹 Define plotting function
plot_obs_vs_pred <- function(actual, predicted, model_name) {
  df_plot <- data.frame(Actual = actual, Predicted = predicted)
  
  ggplot(df_plot, aes(x = Actual, y = Predicted)) +
    geom_point(color = "#0072B2", alpha = 0.6, size = 2) +
    geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "red", size = 1) +
    coord_equal() +
    labs(
      title = paste(model_name, "Model: Observed vs Predicted AGC"),
      x = "Observed AGC",
      y = "Predicted AGC"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.title = element_text(size = 16, face = "bold"),
      axis.text = element_text(size = 14),
      axis.line = element_line(color = "black", size = 0.8),
      panel.grid = element_blank()
    )
}



# 🔹 Import and clean data
df <- read.csv("NDVI_EVI_SAVI_LAI_AGC_Sample_globbiomass.csv")
names(df)[names(df) == "b1"] <- "AGC"

df_clean <- df %>%
  select(AGC, NDVI, EVI, SAVI, LAI) %>%
  drop_na()

# 🔹 Split data
set.seed(123)
train_index <- createDataPartition(df_clean$AGC, p = 0.8, list = FALSE)
train_data <- df_clean[train_index, ]
train_data
test_data  <- df_clean[-train_index, ]
test_data

##############################################
# ✅ 1. Random Forest Model
##############################################
rf_model <- randomForest(AGC ~ NDVI + EVI + SAVI + LAI,
                         data = train_data,
                         ntree = 500,
                         importance = TRUE)

pred_rf <- predict(rf_model, newdata = test_data)
cat("Random Forest Performance:\n")
print(postResample(pred_rf, test_data$AGC))
plot_obs_vs_pred(test_data$AGC, pred_rf, "Random Forest")

##############################################
# ✅ 2. Linear Regression Model
##############################################
lm_model <- lm(AGC ~ NDVI + EVI + SAVI + LAI, data = df_clean)
pred_lm <- predict(lm_model, newdata = test_data)

cat("\nLinear Regression Performance:\n")
print(postResample(pred_lm, test_data$AGC))
plot_obs_vs_pred(test_data$AGC, pred_lm, "Linear Regression")

##############################################
# ✅ 3. Polynomial Regression Model
##############################################
poly_model <- lm(AGC ~ poly(NDVI, 2) + poly(EVI, 2) +
                   poly(SAVI, 2) + poly(LAI, 2), data = df_clean)

pred_poly <- predict(poly_model, newdata = test_data)

cat("\nPolynomial Regression Performance:\n")
print(postResample(pred_poly, test_data$AGC))
plot_obs_vs_pred(test_data$AGC, pred_poly, "Polynomial Regression")

##############################################
# ✅ 4. Generalized Additive Model (GAM)
##############################################
gam_model <- gam(AGC ~ s(NDVI) + s(EVI) + s(SAVI) + s(LAI), data = df_clean)

pred_gam <- predict(gam_model, newdata = test_data)

cat("\nGAM Model Performance:\n")
print(postResample(pred_gam, test_data$AGC))
plot_obs_vs_pred(test_data$AGC, pred_gam, "GAM")

#####################################################################
### AGC  calculation for each year using NDVI, EVI, SAVI, and LAI ###
#####################################################################
install.packages("raster")
install.packages("terra")
library(raster)
library(sp)
library(terra)
library(randomForest)

# Load the multi-band image
veg_stack_2003 <- stack("2003_VIs/2003_VI.tif")

# Rename to match model training
names(veg_stack_2003) <- c("NDVI", "EVI", "SAVI", "LAI")


# Predict AGC using the trained Random Forest model
agc_raster_2003 <- raster::predict(veg_stack_2003, rf_model, progress = "text")

# Save the AGC raster
writeRaster(agc_raster_2003, "AGC_prediction_2003.tif", format = "GTiff", overwrite = TRUE)

# Optional: Plot the result
plot(agc_raster_2003, main = "Predicted AGC for 2003")


##################################################################
###                Uncertainty analysis       ###################
##################################################################

### model based anucertainty
## Crewate raster layer of quantiles
# install.packages("quantregForest")
library(quantregForest)

# Convert raster stack to dataframe
predictors_df <- as.data.frame(veg_stack_2003, na.rm = TRUE)

# Remove rows with NA
predictors_df <- na.omit(predictors_df)

chunk_size <- 10000  # or smaller if needed
n <- nrow(predictors_df)
chunks <- split(predictors_df, ceiling(seq_len(n) / chunk_size))

q_preds_list <- lapply(chunks, function(chunk) {
  predict(qrf_model, newdata = chunk, what = c(0.05, 0.5, 0.95))
})

# Combine the results
q_preds <- do.call(rbind, q_preds_list)

# Convert raster stack to data frame with all cell values
predictors_df_all <- as.data.frame(veg_stack_2003)

# Identify rows (cells) that have complete (non-NA) data across all layers
cell_ids <- which(complete.cases(predictors_df_all))

# Subset only those rows (valid data for prediction)
predictors_df <- predictors_df_all[cell_ids, ]

# Create empty rasters
r_low <- r_median <- r_high <- raster(veg_stack_2003[[1]])
r_low[] <- r_median[] <- r_high[] <- NA  # Initialize with NA

# Fill only valid cells using cell_ids
r_low[cell_ids]    <- q_preds[, 1]  # 5th percentile
r_median[cell_ids] <- q_preds[, 2]  # median
r_high[cell_ids]   <- q_preds[, 3]  # 95th percentile

qrf_stack <- stack(r_low, r_median, r_high)
names(qrf_stack) <- c("Q05", "Q50", "Q95")

# Plot
plot(qrf_stack)

# Save
writeRaster(qrf_stack, filename = "QRF_predictions_2003.tif", format = "GTiff", overwrite = TRUE)

##### Analyse or plot prediction as plots
# Convert to data frame with column names
q_preds_df <- as.data.frame(q_preds)
colnames(q_preds_df) <- c("Q05", "Q50", "Q95")

# Histogram of predicted medians
hist(q_preds_df$Q50, main = "Predicted Median", xlab = "Median Value")

# Compare spread
boxplot(q_preds_df, main = "Quantile Spread")

# Save predictions
write.csv(q_preds_df, "qrf_predictions.csv", row.names = FALSE)

r_uncertainty <- r_high - r_low
plot(r_uncertainty, main = "Prediction Uncertainty (Q95 - Q05)")

summary(values(r_uncertainty))
hist(values(r_uncertainty), main = "Histogram of AGC Prediction Uncertainty", xlab = "Uncertainty")

