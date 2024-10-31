library(tidyverse)
library(infer)
library(randomForest)
library(tidymodels)
library(modelr)
library(yardstick)
library(pROC)
library(xgboost)
library(vip)
setwd('../pstat197a/module1-group9/scripts')
####################################################################
# Question 1
var_names <- read_csv('../data/biomarker-raw.csv', 
                      col_names = F, 
                      n_max = 2, 
                      show_col_types = FALSE,
                      col_select = -(1:2)) %>%
  t() %>%
  as_tibble() %>%
  rename(name = V1, 
         abbreviation = V2) %>%
  na.omit()

biomarker_not_transformed <- read_csv('../data/biomarker-raw.csv', 
                                      skip = 2,
                                      col_select = -2L,
                                      show_col_types = FALSE,
                                      col_names = c('group', 
                                                    'empty',
                                                    pull(var_names, abbreviation),
                                                    'ados'),
                                      na = c('-', '')) %>%
  filter(!is.na(group)) %>%
  # reorder columns
  select(group, ados, everything())


# Select a sample of proteins
set.seed(123)
sample_proteins <- sample(colnames(biomarker_not_transformed)[3:ncol(biomarker_not_transformed)], 10)

# Handle missing values 
biomarker_not_transformed <- na.omit(biomarker_not_transformed)  

# Melt the data for easier plotting
biomarker_long <- biomarker_not_transformed %>%
  pivot_longer(cols = all_of(sample_proteins),
               names_to = "protein",
               values_to = "value")

# Create a histogram facet plot
ggplot(biomarker_long, aes(x = value)) +
  geom_histogram(fill = "steelblue", color = "black") +
  facet_wrap(~ protein, scales = "free_x") +
  labs(title = "Distribution of Selected Proteins",
       x = "Protein Value", y = "Frequency")

####################################################################
# Question 2.

# get names
var_names <- read_csv('data/biomarker-raw.csv', 
                      col_names = F, 
                      n_max = 2, 
                      col_select = -(1:2)) %>%
  t() %>%
  as_tibble() %>%
  rename(name = V1, 
         abbreviation = V2) %>%
  na.omit()

# read in data (without trimming)
bio_clean <- read_csv('data/biomarker-raw.csv', 
                      skip = 2,
                      col_select = -2L,
                      col_names = c('group', 
                                    'empty',
                                    pull(var_names, abbreviation),
                                    'ados'),
                      na = c('-', '')) %>%
  filter(!is.na(group)) %>%
  # log transform, center and scale
  mutate(across(.cols = -c(group, ados), 
                ~ (scale(log10(.x))[, 1]))) %>%
  # reorder columns
  select(group, ados, everything())

# computer number of outliers for each individual
bio_outliers <- bio_clean %>%
  mutate(
    id = row_number(),
    .before = group
  ) %>%
  pivot_longer(
    -c(id, group, ados),
    names_to = "protein",
    values_to = "level"
  ) %>%
  group_by(
    id
  ) %>%
  summarize(
    group = first(group),
    ados = first(ados),
    n.outlier = sum(abs(level) > 3)
  )

# extract top 10 observations with largest number of outliers
bio_outliers %>%
  arrange(
    desc(n.outlier)
  ) %>%
  head(n = 10)

# generate histogram of number of outliers
bio_outliers %>%
  ggplot() +
  geom_histogram(
    aes(x = n.outlier)
  ) +
  labs(
    title = "Histogram of Number of Outliers",
    x = "Number of Outliers",
    y = "Count"
  ) +
  theme_minimal()

# generate histogram of number of outliers by group
bio_outliers %>%
  ggplot() +
  geom_histogram(
    aes(
      x = n.outlier,
      after_stat(density),
      fill = group
    ),
  ) +
  facet_wrap(vars(group)) +
  labs(
    title = "Histogram of Number of Outliers by Group",
    x = "Number of Outliers",
    y = "Frequency",
    fill = "Group"
  ) +
  theme_minimal()

####################################################################
# Question 3A


load('../data/biomarker-clean.RData')
library(randomForest)

# Split the data into training and testing sets
set.seed(101422)
biomarker_split <- initial_split(biomarker_clean, prop = 0.8)
train_data <- training(biomarker_split)
test_data <- testing(biomarker_split)

## Multiple Testing
test_fn <- function(.df){
  t_test(.df,
         formula = level ~ group,
         order = c('ASD', 'TD'),
         alternative = 'two-sided',
         var.equal = F)
}

ttests_out <- train_data %>%
  # drop ADOS score
  select(-ados) %>%
  # arrange in long format
  pivot_longer(-group,
               names_to = 'protein',
               values_to = 'level') %>%
  # nest by protein
  nest(data = c(level, group)) %>%
  # compute t tests
  mutate(ttest = map(data, test_fn)) %>%
  unnest(ttest) %>%
  # sort by p-value
  arrange(p_value) %>%
  # multiple testing correction
  mutate(m = n(),
         hm = log(m) + 1/(2*m) - digamma(1),
         rank = row_number(),
         p.adj = m*hm*p_value/rank)

# select significant proteins based on adjusted p-value
proteins_s1 <- ttests_out %>%
  slice_min(p.adj, n = 10) %>%
  pull(protein)

## Random Forest 
# store predictors and response separately
predictors <- train_data %>% 
  select(-c(group, ados))

response <- train_data %>% pull(group) %>% factor()

# fit RF
set.seed(123)
rf_out <- randomForest(x = predictors,
                       y = response,
                       ntree = 1000,
                       importance = T)

# check errors 
rf_out$confusion

# compute importance scores
proteins_s2 <- rf_out$importance %>% 
  as_tibble() %>%
  mutate(protein = rownames(rf_out$importance)) %>%
  slice_max(MeanDecreaseGini, n = 10) %>%
  pull(protein)


## Logistic Regression
# select subset of interest
proteins_sstar <- intersect(proteins_s1, proteins_s2)

biomarker_sstar <- train_data %>%
  select(group, any_of(proteins_sstar)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)

biomarker_sstar_test <- test_data %>%
  select(group, any_of(proteins_sstar)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)

# fit logistic regression model to training set
fit <- glm(class ~ ., 
           data = biomarker_sstar, 
           family = 'binomial')

# evaluate errors on test set
class_metrics <- metric_set(sensitivity, 
                            specificity, 
                            accuracy,
                            roc_auc)

biomarker_sstar_test %>%
  add_predictions(fit, type = 'response') %>%
  mutate(est = as.factor(pred > 0.5), tr_c = as.factor(class)) %>%
  class_metrics(estimate = est,
                truth = tr_c, pred,
                event_level = 'second')

####################################################################
# Question 3b

## MULTIPLE TESTING
####################

# function to compute tests
test_fn <- function(.df){
  t_test(.df, 
         formula = level ~ group,
         order = c('ASD', 'TD'),
         alternative = 'two-sided',
         var.equal = F)
}

ttests_out <- biomarker_clean %>%
  # drop ADOS score
  select(-ados) %>%
  # arrange in long format
  pivot_longer(-group, 
               names_to = 'protein', 
               values_to = 'level') %>%
  # nest by protein
  nest(data = c(level, group)) %>% 
  # compute t tests
  mutate(ttest = map(data, test_fn)) %>%
  unnest(ttest) %>%
  # sort by p-value
  arrange(p_value) %>%
  # multiple testing correction
  mutate(m = n(),
         hm = log(m) + 1/(2*m) - digamma(1),
         rank = row_number(),
         p.adj = m*hm*p_value/rank)

# select significant proteins, n = 10
proteins_s1_10 <- ttests_out %>%
  slice_min(p.adj, n = 10) %>%
  pull(protein)

# select significant proteins, n = 15
proteins_s1_15 <- ttests_out %>%
  slice_min(p.adj, n = 15) %>%
  pull(protein)

# select significant proteins, n = 20
proteins_s1_20 <- ttests_out %>%
  slice_min(p.adj, n = 20) %>%
  pull(protein)

# select significant proteins, n = 25
proteins_s1_25 <- ttests_out %>%
  slice_min(p.adj, n = 25) %>%
  pull(protein)

## RANDOM FOREST
##################

# store predictors and response separately
predictors <- biomarker_clean %>%
  select(-c(group, ados))

response <- biomarker_clean %>% pull(group) %>% factor()

# fit RF
set.seed(101422)
rf_out <- randomForest(x = predictors, 
                       y = response, 
                       ntree = 1000, 
                       importance = T)

# check errors
rf_out$confusion

# compute importance scores, n = 10
proteins_s2_10 <- rf_out$importance %>% 
  as_tibble() %>%
  mutate(protein = rownames(rf_out$importance)) %>%
  slice_max(MeanDecreaseGini, n = 10) %>%
  pull(protein)

# protein selection, n = 15
proteins_s2_15 <- rf_out$importance %>%
  as_tibble() %>%
  mutate(protein = rownames(rf_out$importance)) %>%
  slice_max(MeanDecreaseGini, n = 15) %>%
  pull(protein)

# protein selection, n = 20
proteins_s2_20 <- rf_out$importance %>%
  as_tibble() %>%
  mutate(protein = rownames(rf_out$importance)) %>%
  slice_max(MeanDecreaseGini, n = 20) %>%
  pull(protein)

# protein selection, n = 25
proteins_s2_25 <- rf_out$importance %>%
  as_tibble() %>%
  mutate(protein = rownames(rf_out$importance)) %>%
  slice_max(MeanDecreaseGini, n = 25) %>%
  pull(protein)

## LOGISTIC REGRESSION
#######################

## n = 10

# select subset of interest
proteins_sstar <- intersect(proteins_s1_10, proteins_s2_10)

biomarker_sstar <- biomarker_clean %>%
  select(group, any_of(proteins_sstar)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)

# partition into training and test set
set.seed(101422)
biomarker_split <- biomarker_sstar %>%
  initial_split(prop = 0.8)

# fit logistic regression model to training set
fit <- glm(class ~ ., 
           data = training(biomarker_split), 
           family = 'binomial')

# evaluate errors on test set
class_metrics <- metric_set(sensitivity, 
                            specificity, 
                            accuracy,
                            roc_auc)

metric_10 <- testing(biomarker_split) %>%
  add_predictions(fit, type = 'response') %>%
  mutate(est = as.factor(pred > 0.5), tr_c = as.factor(class)) %>%
  class_metrics(estimate = est,
                truth = tr_c, pred,
                event_level = 'second')

## n = 15

# select subset of interest
proteins_sstar <- intersect(proteins_s1_15, proteins_s2_15)

biomarker_sstar <- biomarker_clean %>%
  select(group, any_of(proteins_sstar)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)

# partition into training and test set
set.seed(101422)
biomarker_split <- biomarker_sstar %>%
  initial_split(prop = 0.8)

# fit logistic regression model to training set
fit <- glm(class ~ ., 
           data = training(biomarker_split), 
           family = 'binomial')

# evaluate errors on test set
class_metrics <- metric_set(sensitivity, 
                            specificity, 
                            accuracy,
                            roc_auc)

metric_15 <- testing(biomarker_split) %>%
  add_predictions(fit, type = 'response') %>%
  mutate(est = as.factor(pred > 0.5), tr_c = as.factor(class)) %>%
  class_metrics(estimate = est,
                truth = tr_c, pred,
                event_level = 'second')

## n = 20

# select subset of interest
proteins_sstar <- intersect(proteins_s1_20, proteins_s2_20)

biomarker_sstar <- biomarker_clean %>%
  select(group, any_of(proteins_sstar)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)

# partition into training and test set
set.seed(101422)
biomarker_split <- biomarker_sstar %>%
  initial_split(prop = 0.8)

# fit logistic regression model to training set
fit <- glm(class ~ ., 
           data = training(biomarker_split), 
           family = 'binomial')

# evaluate errors on test set
class_metrics <- metric_set(sensitivity, 
                            specificity, 
                            accuracy,
                            roc_auc)

metric_20 <- testing(biomarker_split) %>%
  add_predictions(fit, type = 'response') %>%
  mutate(est = as.factor(pred > 0.5), tr_c = as.factor(class)) %>%
  class_metrics(estimate = est,
                truth = tr_c, pred,
                event_level = 'second')

## n = 25

# select subset of interest
proteins_sstar <- intersect(proteins_s1_25, proteins_s2_25)

biomarker_sstar <- biomarker_clean %>%
  select(group, any_of(proteins_sstar)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)

# partition into training and test set
set.seed(101422)
biomarker_split <- biomarker_sstar %>%
  initial_split(prop = 0.8)

# fit logistic regression model to training set
fit <- glm(class ~ ., 
           data = training(biomarker_split), 
           family = 'binomial')

# evaluate errors on test set
class_metrics <- metric_set(sensitivity, 
                            specificity, 
                            accuracy,
                            roc_auc)

metric_25 <- testing(biomarker_split) %>%
  add_predictions(fit, type = 'response') %>%
  mutate(est = as.factor(pred > 0.5), tr_c = as.factor(class)) %>%
  class_metrics(estimate = est,
                truth = tr_c, pred,
                event_level = 'second')

# add column for n in error dataframes
metric_10 <- metric_10 %>%
  mutate(
    n = 10
  )

metric_15 <- metric_15 %>%
  mutate(
    n = 15
  )

metric_20 <- metric_20 %>%
  mutate(
    n = 20
  )

metric_25 <- metric_25 %>%
  mutate(
    n = 25
  )

# create single dataframe for all model errors
metrics <- rbind(metric_10, metric_15, metric_20, metric_25)

# plot model performance against number of top proteins selected
metrics %>%
  ggplot(
    aes(
      x = n,
      y = .estimate
    )
  ) +
  geom_point() +
  geom_line() +
  facet_wrap(
    vars(.metric)
  ) +
  labs(
    title = "Model Performance Metrics by Number of Top Proteins Selected",
    x = "Number of Top Proteins Selected",
    y = "Estimate"
  ) +
  theme_minimal()

####################################################################
# Question 3c
## MULTIPLE TESTING
####################

# function to compute tests
test_fn <- function(.df){
  t_test(.df, 
         formula = level ~ group,
         order = c('ASD', 'TD'),
         alternative = 'two-sided',
         var.equal = F)
}

ttests_out <- biomarker_clean %>%
  # drop ADOS score
  select(-ados) %>%
  # arrange in long format
  pivot_longer(-group, 
               names_to = 'protein', 
               values_to = 'level') %>%
  # nest by protein
  nest(data = c(level, group)) %>% 
  # compute t tests
  mutate(ttest = map(data, test_fn)) %>%
  unnest(ttest) %>%
  # sort by p-value
  arrange(p_value) %>%
  # multiple testing correction
  mutate(m = n(),
         hm = log(m) + 1/(2*m) - digamma(1),
         rank = row_number(),
         p.adj = m*hm*p_value/rank)

# select significant proteins
proteins_s1 <- ttests_out %>%
  slice_min(p.adj, n = 10) %>%
  pull(protein)

## RANDOM FOREST
##################

# store predictors and response separately
predictors <- biomarker_clean %>%
  select(-c(group, ados))

response <- biomarker_clean %>% pull(group) %>% factor()

# fit RF
set.seed(101422)
rf_out <- randomForest(x = predictors, 
                       y = response, 
                       ntree = 1000, 
                       importance = T)

# check errors
rf_out$confusion

# compute importance scores
proteins_s2 <- rf_out$importance %>% 
  as_tibble() %>%
  mutate(protein = rownames(rf_out$importance)) %>%
  slice_max(MeanDecreaseGini, n = 10) %>%
  pull(protein)

## LOGISTIC REGRESSION
#######################

# select subset of interest by fuzzy intersecting
proteins_sstar <- unique(proteins_s1, proteins_s2)

biomarker_sstar <- biomarker_clean %>%
  select(group, any_of(proteins_sstar)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)

# partition into training and test set
set.seed(101422)
biomarker_split <- biomarker_sstar %>%
  initial_split(prop = 0.8)

# fit logistic regression model to training set
fit <- glm(class ~ ., 
           data = training(biomarker_split), 
           family = 'binomial')

# evaluate errors on test set
class_metrics <- metric_set(sensitivity, 
                            specificity, 
                            accuracy,
                            roc_auc)

testing(biomarker_split) %>%
  add_predictions(fit, type = 'response') %>%
  mutate(est = as.factor(pred > 0.5), tr_c = as.factor(class)) %>%
  class_metrics(estimate = est,
                truth = tr_c, pred,
                event_level = 'second')


#####################################################################
# Question 4

## BOOSTED TREES
#############################
set.seed(456)
# modifying group into 1 and 0 
biomarker_bt <- biomarker_clean %>% 
  mutate(ASD = as.factor(group == 'ASD'))

# splitting data into training and testing
biomarker_split_bt <- initial_split(biomarker_bt, 0.8)
biomarker_train_bt <- training(biomarker_split_bt)
biomarker_test_bt <- testing(biomarker_split_bt)

biomarker_recipe <- recipe(ASD ~ ., data = biomarker_train_bt) %>%
  step_rm(c("ados", "group")) %>% 
  step_center() %>% 
  step_scale() 

biomarker_folds <- vfold_cv(biomarker_train_bt, v = 10)
# creating model 
bt_mod <- boost_tree(mtry = tune(), 
                     trees = tune(), 
                     learn_rate = tune()) %>%
  set_engine("xgboost") %>% 
  set_mode("classification")


# boosted trees workflow
bt_workflow <- workflow() %>%
  add_recipe(biomarker_recipe) %>%
  add_model(bt_mod)

# Define a grid for Boosted Trees
bt_grid <- grid_regular(mtry(range = c(4,12)), 
                        trees(range = c(350,500)), 
                        learn_rate(range = c(-6,-3)), 
                        levels = 5)

# Tune Boosted Trees
bt_results <- bt_workflow %>%
  tune_grid(resamples = biomarker_folds, 
            grid = bt_grid)

show_best(bt_results, metric = "roc_auc")

best_bt <- select_by_one_std_err(bt_results, 
                                 metric = "roc_auc",
                                 mtry,
                                 trees,
                                 learn_rate)

final_bt_model <- finalize_workflow(bt_workflow, best_bt)
final_bt_model <- fit(final_bt_model, biomarker_train_bt)
final_bt_model %>% extract_fit_parsnip() %>% 
  vip() +
  theme_minimal()


bt_panel <- final_bt_model %>% 
  extract_fit_parsnip() %>% 
  vip::vi() %>%               
  arrange(desc(Importance)) %>%  
  pull(Variable) %>% 
  head(10)

# check the top 10 proteins
bt_panel

# use only top 10 proteins in a new model and fit to training data 
panel_recipe <- biomarker_recipe <- recipe(ASD ~ ., data = biomarker_train_bt) %>%
  step_rm(c("ados", "group")) %>% 
  update_role(all_of(bt_panel), new_role = "predictor") %>%
  step_center() %>% 
  step_scale() 

panel_workflow <- workflow() %>%
  add_recipe(panel_recipe) %>%
  add_model(bt_mod)

panel_model <- finalize_workflow(panel_workflow, best_bt)
panel_model <- fit(panel_model, biomarker_train_bt)

# calculate sensitivity specificity accuracy to benchmark
multi_metric <- metric_set(sensitivity, specificity, accuracy)
augment(panel_model, new_data = biomarker_test_bt) %>% 
  multi_metric(truth = ASD, estimate = .pred_class)

