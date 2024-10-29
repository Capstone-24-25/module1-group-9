library(tidyverse)
library(infer)
library(randomForest)
library(tidymodels)
library(modelr)
library(yardstick)
library(pROC)
library(xgboost)
library(vip)
setwd("../module1-group-9") 
load('data/biomarker-clean.RData')
set.seed(200)
## BOOSTED TREES
#############################
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
  vip::vi() %>%                 # Get variable importance scores
  arrange(desc(Importance)) %>%  # Sort by importance               # Select top 10
  pull(Variable) %>% 
  head(10)
  
bt_panel

panel_recipe <- biomarker_recipe <- recipe(ASD ~ ., data = biomarker_train_bt) %>%
  step_rm(c("ados", "group")) %>% 
  step_select(bt_panel) %>% 
  step_center() %>% 
  step_scale() 

panel_workflow <- workflow() %>%
  add_recipe(panel_recipe) %>%
  add_model(bt_mod)

panel_model <- finalize_workflow(panel_workflow, best_bt)
panel_model <- fit(panel_model, biomarker_train_bt)

multi_metric <- metric_set(sensitivity, specificity, accuracy)
augment(panel_model, new_data = biomarker_test_bt) %>% 
  multi_metric(truth = ASD, estimate = .pred_class)

