# Add this to your R code:
.libPaths(c("/projappl/project_2009874/project_rpackages_4.3.2", .libPaths()))
libpath <- .libPaths()[1]
install.packages("performance")
install.packages("MuMIn")

# This command can be used to check that the folder is now visible:
.libPaths() # It should be first on the list
# Load necessary libraries
library(dplyr)
library(performance)
library(glmmTMB)
library(DHARMa)
library(MuMIn)
library(tools)
library(readxl)

# Define the dataset filenames with full paths
file_path <- "Desktop/species/"  # Adjust this path as needed
datasets <- c("mountainhare3.xlsx", "pinemarten3.xlsx", "forestreindeer3.xlsx", "lynx3.xlsx", "wolf3.xlsx", "wolverine3.xlsx", "otter3.xlsx", "moose3.xlsx","brownhare3.xlsx", "ermine3.xlsx", "roedeer3.xlsx", "mink3.xlsx", "redfox3.xlsx", "redsquirrel3.xlsx", "racoondog3.xlsx", "wtaileddeer3.xlsx", "weasel3.xlsx")
datasets <- paste0(file_path, datasets)

# Initialize a list to store models and their R^2 values
models <- list()
r2_values <- list()
collinearity_checks <- list()

# Function to process each dataset
process_dataset <- function(filename) {
  # Extract the species name from the filename
  species <- file_path_sans_ext(basename(filename))
  
  # Read the dataset
  data <- read_excel(filename)
  
  # Filter the data
  data <- data[which(data$snow_days <= 5), ]
  # Define factor variables
  data$triangles <- as.factor(data$triangle_id)
  data$year_f <- as.factor(data$year)
  data$riistakeskus <- as.factor(data$riistakeskus)
  data$rhy_id <- as.factor(data$rhy_id)
  
  # Scale numerical variables
  data$year.c <- data$year - mean(data$year)
  data <- na.omit(data)
  data$latitude.100km <- (data$triangle_y - mean(data$triangle_y)) / 100000
  data$longitude.100km <- (data$triangle_x - mean(data$triangle_x)) / 100000
  data$FLII <- data$indexvalue - mean(data$indexvalue)
  data()$forest_cover <- data$forest_pro - mean(data$forest_pro)
  data$canopy <- data()$canopy - mean(data$canopy)
#NB

model<- glmmTMB(n_tracks ~ FLII + forest_cover + canopy + latitude.100km + year_f+log(snow_days)+(1|riistakeskus/rhy_id/triangles), 
                   offset=log(covered_di/10000),
                   family =nbinom2(link = "log"),
                   data = data) 
summary(model_nb)


# Store the model in the list with the species name
models[[species]] <- model

# Print the summary of the model with species name
cat("Summary of the model for", species, ":\n")
print(summary(model))

# Calculate R^2 for the model
r2_nb <- tryCatch({
  r.squaredGLMM(model)
}, error = function(e) {
  cat("Error calculating R^2 for", species, ":", e$message, "\n")
  NULL
})
if (!is.null(r2_nb)) {
  r2_values[[species]] <- r2_nb
  cat("R^2 for the model for", species, ":\n")
  print(r2_nb)
}

# Check for collinearity
collinearity_check <- tryCatch({
  check_collinearity(model)
}, error = function(e) {
  cat("Error checking collinearity for", species, ":", e$message, "\n")
  NULL
})
if (!is.null(collinearity_check)) {
  collinearity_checks[[species]] <- collinearity_check
  cat("Collinearity check for the model for", species, ":\n")
  print(collinearity_check)
}
}

# Apply the function to each dataset one by one
process_dataset(datasets[1])  # brownhare3.xlsx
process_dataset(datasets[2])  # ermine3.xlsx
process_dataset(datasets[3])  # roedeer3.xlsx
process_dataset(datasets[4])  # mink3.xlsx
process_dataset(datasets[5])  # redfox3.xlsx
process_dataset(datasets[6])  # redsquirrel3.xlsx
process_dataset(datasets[7])  # racoondog3.xlsx
process_dataset(datasets[8])  # wtaileddeer3.xlsx
process_dataset(datasets[9])  # weasel3.xlsx
process_dataset(datasets[1])  # mountainhare3.xlsx
process_dataset(datasets[2])  # pinemarten3.xlsx
process_dataset(datasets[3])  # forestreindeer3.xlsx
process_dataset(datasets[4])  # lynx3.xlsx
process_dataset(datasets[5])  # wolf3.xlsx
process_dataset(datasets[6])  # wolverine3.xlsx
process_dataset(datasets[7])  # otter3.xlsx
process_dataset(datasets[8])  # moose3.xlsx

# Print all R^2 values and collinearity checks at the end
cat("All R^2 values:\n")
print(r2_values)
cat("All collinearity checks:\n")
print(collinearity_checks)

# The models list now contains models for each species with corresponding names
models
