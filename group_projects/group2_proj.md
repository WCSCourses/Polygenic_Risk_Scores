---
title: "Leveraging pre-existing PGS to improve prediction of insulin resistance in African populations"
author: "Group2"
date: "`r Sys.Date()`"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```
### Loading the neccesary Libraries

```{r load_libraries, message=FALSE}
library(glmnet)
library(dplyr)
library(caret)
library(ggplot2)
library(gridExtra)
library(knitr)
library(tibble)
```


### Loading the dataset

```{r load_dataset, message=FALSE, warning=FALSE}
homair_df <- read.table("../Downloads/homairdta.txt", sep = " ", header = TRUE)
```

## Dataset Descriptive statistics
### Demographic and Metabolic Analysis





#### 1. Basic Sex Distribution
```{r sex_barplot_base}
# Simple barplot using base R graphics
homair_df$sex<- factor(homair_df$sex, levels = c(1,2), labels = c("Male", "Female"))
barplot(table(homair_df$sex),
        main = "Sex Distribution",
        col = c("lightpink", "lightblue"),
        ylab = "Frequency")
```

### 2. Enhanced Visualizations (ggplot2)

#### Age Distribution
```{r age_distribution}
p1 <- ggplot(homair_df, aes(x = age)) + 
  geom_histogram(fill = "steelblue", bins = 30, color = "white") +
  labs(title = "Age Distribution", x = "Age (years)", y = "Count") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
```

#### Sex Distribution (Enhanced)
```{r sex_distribution}

p2 <- ggplot(homair_df, aes(x = sex, fill = sex)) + 
  geom_bar() +
  scale_fill_manual(values = c("Male" = "#3498db", "Female" = "#e74c3c")) +
  labs(title = "Sex Distribution", x = "Sex", y = "Count") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none")
```

#### HOMA-IR by Sex
```{r homair_distribution}
p3 <- ggplot(homair_df, aes(x = sex, y = HOMAIR, fill = sex)) + 
  geom_boxplot() +
  scale_fill_manual(values = c("Male" = "#3498db", "Female" = "#e74c3c")) +
  labs(x = "Sex", y = "HOMA-IR") +
  theme_minimal() +
  ggtitle("HOMA-IR by Sex") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
```

### 3. Combined Visualization
```{r combined_plots, fig.width=10, fig.height=6}
grid.arrange(p1, p2, p3, ncol = 2,
             top = "Demographic and Metabolic Characteristics")
```

### 4. Frequency Table
```{r frequency_table}
sex_table <- table(homair_df$sex)
print(sex_table)
```




###
### Define PGS and covariate columns

```{r define_columns, message=FALSE, warning=FALSE }
pgs_cols <- grep("^PGS00", colnames(homair_df), value = TRUE)
covariate_cols <- c("age", "sex", paste0("PC", 1:10))
```


### Extract predictors and outcome

``` {r split_dataset, message=FALSE, warning=FALSE }
X_all <- homair_df[, c(covariate_cols, pgs_cols)]
y_all <- homair_df$HOMAIR

set.seed(123)
train_index <- createDataPartition(y_all, p = 0.8, list = FALSE)
X_train_raw <- X_all[train_index, ]
X_test_raw  <- X_all[-train_index, ]
y_train <- y_all[train_index]
y_test  <- y_all[-train_index]

```

### Data Normalisation the data

``` {r normalisation, message=FALSE, warning=FALSE}
# Create model matrix to handle factors like 'sex'
X_train_mat <- model.matrix(~ ., data = X_train_raw)[, -1]  # remove intercept
X_test_mat  <- model.matrix(~ ., data = X_test_raw)[, -1]

# Standardize using training set statistics
means <- apply(X_train_mat, 2, mean)
sds   <- apply(X_train_mat, 2, sd)
X_train_std <- scale(X_train_mat, center = means, scale = sds)

```
### Fit the Elastic Net model
``` {r model_fitting, message=FALSE, warning=FALSE }
cv_fit <- cv.glmnet(X_train_std, y_train, alpha = 0.5, nfolds = 5, standardize = FALSE)
coefs_std <- as.vector(coef(cv_fit, s = "lambda.min"))[-1]
adjusted_weights <- coefs_std / sds
combined_pred_test <- X_test_mat %*% adjusted_weights
```


### Evaluate model
```{r evaluate-model, echo=TRUE}
r2 <- cor(combined_pred_test, y_test)^2
mse <- mean((combined_pred_test - y_test)^2)

cat("Elastic Net PRSmix with covariates:\n")
cat("R-squared:", round(r2, 4), "\n")
cat("MSE:", round(mse, 4), "\n")

pred_vs_actual_df <- data.frame(
  Actual = y_test,
  Predicted = as.vector(combined_pred_test)
)
```
### Compute R² for each individual PGS predictor
```{r}

pgs_r2 <- sapply(pgs_cols, function(pgs_name) {
  single_pgs <- homair_df[-train_index, pgs_name]
  cor(single_pgs, y_test)^2
})

```
### Compare Combined Model to Individual PGS Predictors


``` {r Combining_results, message=FALSE, warning=FALSE}
pgs_r2_df <- data.frame(
  Model = c(pgs_cols, "Combined_EN"),
  R2 = c(pgs_r2, r2)
) %>%
  arrange(desc(R2))
ggplot(pgs_r2_df, aes(x = R2, y = reorder(Model, R2), fill = Model == "Combined_EN")) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("gray70", "orange")) +
  coord_flip() +
  labs(
    title = "Individual PGS vs Combined Elastic Net Model",
    y = "PGS Models",
    x = "Prediction Accuracy (R²)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

```