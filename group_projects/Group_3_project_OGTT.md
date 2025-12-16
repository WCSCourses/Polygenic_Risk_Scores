
# Group 3: Can PRSmix Predict if You Will Say No To the Doughnut? A Genetic Crystall Ball for Type 2 Diabetes in African Populations (Using the Oral Glucose Tolerance Test aka, The Sugar Rush Test

## Author: Himran Moundi, Andrea Sene, Olga Tendo, Pheziwe Mshunqwane, Nto Johnson Nto

## 🧾 Introduction

The aim of this analysis is to assess the predictive power of polygenic risk scores (PRS) and traditional covariates (such as age, sex, and ancestry principal components) for 2-hour post-load blood glucose levels—a key marker for glucose regulation and Type 2 diabetes risk.

To achieve this, we apply **Elastic Net regression**, a regularized linear modeling technique that combines both Lasso and Ridge penalties. Elastic Net is particularly suitable for datasets with many correlated predictors, such as multiple PRS.

The workflow includes:
- Loading and preprocessing the data
- Log-transforming the outcome variable (2hrglc) to normalize it
- Standardizing the predictors
- Splitting the data into training and testing sets
- Using cross-validation to tune the Elastic Net's alpha and lambda parameters
- Evaluating model performance on the test set using R² and RMSE
- Comparing the incremental R² of each individual PRS to the Elastic Net model
- Visualizing feature importance and predictive accuracy

This approach helps identify which genetic and demographic factors are most predictive of glucose response, and whether a combined model like Elastic Net outperforms individual PRS in explaining variation in glucose levels.

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```


## 📦 Load Required Packages

```{r}
# Load all necessary packages for modeling, visualization, and data manipulation
library(glmnet)   # For fitting Elastic Net
library(caret)    # For training/tuning Elastic Net via cross-validation
library(dplyr)    # For data wrangling
library(ggplot2)  # For plotting
library(readr)    # For reading data
library(vip)      # For visualizing feature importance
```

## 📂 Load and Prepare Dataset

```{r}
# Load dataset
data <- read.table("2hrglcdta.txt", header = TRUE)

# Rename column for easier use
names(data)[names(data) == "X2hrglc"] <- "2hrglc"
```

# Convert outcome to numeric and log-transform for normalization
```{r}
# Convert the glucose outcome to numeric
data$`2hrglc` <- as.numeric(data$`2hrglc`)

# Log-transform the outcome to reduce skewness
data$log_2hrglc <- log(data$`2hrglc`)
```

## 📊 Visualize Distribution of Outcome
```{r}
# Histogram of original 2hr glucose
ggplot(data, aes(x = `2hrglc`)) +
  geom_histogram(bins = 30, fill = "grey") +
  ggtitle("Original 2hr Glucose Distribution") +
  theme_minimal()

# Histogram of log-transformed glucose
ggplot(data, aes(x = log_2hrglc)) +
  geom_histogram(bins = 30, fill = "skyblue") +
  ggtitle("Log-transformed 2hr Glucose Distribution") +
  theme_minimal()
```

![Alt text](https://github.com/WCSCourses/PRS_2025/blob/main/group_projects/group3_images/image5.png)

## 🔍 Define PRS and Covariates
```{r}
# Select all PRS columns using pattern matching
prs_cols <- grep("^PGS", names(data), value = TRUE)

# Select covariates: age, sex, and PCs
covars <- data %>% select(age, sex, PC1:PC10)

# Convert sex to numeric (e.g., 0 and 1)
covars$sex <- as.numeric(as.factor(covars$sex)) - 1
```

# Create model matrix (X) and outcome (y)
```{r}
# Combine covariates and PRS into predictor matrix X
X <- cbind(as.matrix(covars), as.matrix(data[, prs_cols]))
y <- data$log_2hrglc  # Use log-transformed outcome
```

## ✂️ Split Data into Training and Testing
```{r}
# For reproducibility
set.seed(42)

# Random 80% train / 20% test split
n <- nrow(data)
train_index <- sample(1:n, size = 0.8 * n, replace = FALSE)

# Create train/test sets
X_train <- X[train_index, ]
y_train <- y[train_index]
X_test <- X[-train_index, ]
y_test <- y[-train_index]
```

## ⚙️ Preprocess: Standardize Numeric Predictors
```{r}
# Identify numeric columns to standardize
num_cols <- colnames(X_train)

# Fit preprocessing model (center and scale)
pp <- preProcess(as.data.frame(X_train)[, num_cols], method = c("center", "scale"))

# Apply preprocessing to train and test sets
X_train_scaled <- predict(pp, as.data.frame(X_train))
X_test_scaled  <- predict(pp, as.data.frame(X_test))
```

## 🔁 Train Elastic Net with caret (Optimize Alpha & Lambda)
```{r}
# Prepare data for caret
train_data <- X_train_scaled
train_data$y <- y_train

# Set up 5-fold cross-validation
set.seed(123)
ctrl <- trainControl(method = "cv", number = 5, verboseIter = TRUE)

# Train Elastic Net model, letting caret tune alpha and lambda
elastic_model <- train(
  y ~ .,
  data = train_data,
  method = "glmnet",
  trControl = ctrl,
  tuneLength = 10
)
```


```{r}
# Print best model results
print(elastic_model)
plot(elastic_model)
```

```{r}
# Extract best alpha and lambda values
best_alpha <- elastic_model$bestTune$alpha
best_lambda <- elastic_model$bestTune$lambda
cat("Best alpha:", best_alpha, "\n")
cat("Best lambda:", best_lambda, "\n")
```

## 🧠 Predict and Evaluate on Test Set
```{r}
# Predict on test set using best model
y_pred <- predict(elastic_model, newdata = X_test_scaled)

# Calculate R² and RMSE
SSE <- sum((y_test - y_pred)^2)
SST <- sum((y_test - mean(y_test))^2)
R2 <- 1 - (SSE / SST)
RMSE <- sqrt(mean((y_test - y_pred)^2))

# Output performance metrics
cat("Test R²:", round(R2, 4), "\n")
cat("Test RMSE:", round(RMSE, 4), "\n")
```

## 📈 Observed vs Predicted Plot
```{r}
# Create data frame for plotting
df_plot <- data.frame(True = y_test, Predicted = y_pred)

# Plot true vs predicted values
ggplot(df_plot, aes(x = True, y = Predicted)) +
  geom_point(color = "#2c7bb6") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(title = "Elastic Net Prediction: Log 2hr Glucose",
       x = "Observed",
       y = "Predicted")
```



## 📜 Model Coefficients and Variable Importance
```{r}
# Extract final model coefficients at best lambda
final_model <- elastic_model$finalModel
coefs <- coef(final_model, s = best_lambda)

# Convert to data frame and filter out zero coefficients
coefs_df <- as.data.frame(as.matrix(coefs))
coefs_df <- coefs_df[coefs_df != 0, , drop = FALSE]
coefs_df

# Plot top 20 important features
vip(final_model, num_features = 32)
```

## 📊 Compare Incremental R² of PRS vs Elastic Net
```{r}
# Function to compute incremental R² of predictor beyond covariates
compute_r2 <- function(predictor, data, outcome = "2hrglc") {
    outcome <- paste0("`", outcome, "`")
    formula_null <- as.formula(paste(outcome, "~ age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"))
    formula_full <- as.formula(paste(outcome, "~", predictor, "+ age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"))
    
    null_model <- lm(formula_null, data = data)
    full_model <- lm(formula_full, data = data)
    
    r2_null <- summary(null_model)$r.squared
    r2_full <- summary(full_model)$r.squared
    return(r2_full - r2_null)
  }
```


```{r}
# Compute R² improvement for each PRS
prs_names <- grep("^PGS", names(data), value = TRUE)
  prs_r2 <- sapply(prs_names, function(p) compute_r2(p, data = data, outcome = "log_2hrglc"))
  prs_r2
```


```{r}
# Add Elastic Net R² for comparison
prs_r2["ElasticNet"] <- R2
prs_r2
```


```{r}
# Sort and plot R² comparisons
prs_r2_sorted <- sort(prs_r2, decreasing = TRUE)
  barplot(prs_r2_sorted,
          las = 2,
          col = ifelse(names(prs_r2_sorted) == "ElasticNet", "red", "steelblue"),
          main = "Incremental R²: Individual PRS vs Elastic Net",
          ylab = "Incremental R²",
          cex.names = 0.7)
  legend("topright", legend = c("Individual PRS", "Elastic Net"), fill = c("steelblue", "red"), bty = "n")
```

![Alt text](https://github.com/WCSCourses/PRS_2025/blob/main/group_projects/group3_images/image5.png)
