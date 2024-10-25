library(tidyverse)

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
