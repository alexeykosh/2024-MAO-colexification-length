library(tidyverse)
library(lme4)
library(marginaleffects)
library(DHARMa)
library(lattice)


#Started running today 2024-04-08 at 12.


#0. Data processing #####

# Data
## Extracting
lexibank_all <- read.csv('data/forms_total.csv')
## Preprocessing
lexibank_all$Concepticon_ID <- as.numeric(lexibank_all$Concepticon_ID)
lexibank_all_colex <- lexibank_all %>%
  dplyr::select(Glottocode, Form,
         Concepticon_ID, Length,
         Family) %>%
  group_by(Glottocode, Form) %>%
  mutate(colexification = n_distinct(Concepticon_ID)) %>%
  group_by(Glottocode) %>%
  mutate(n_forms = n()) %>%
  distinct(Form, .keep_all = TRUE) %>%
  mutate(colexification_bin = ifelse(colexification > 1, 1, 0))
logistic_data <- lexibank_all_colex %>%
  group_by(Glottocode) %>%
  mutate(z_length = (Length - mean(Length)) / sd(Length))  %>%
  ungroup()  %>%
  mutate(z_n_forms = (n_forms - mean(n_forms)) / sd(n_forms))
## Removing non-labelled families
logistic_data <- logistic_data[-which(logistic_data$Family == ""), ]


#1. Basic analysis, as preregistered ######

#1.1. Regression models ######

# Statistical analysis
## Null model without length
model_shuffle_null <- glmer(colexification_bin ~ 1 + z_n_forms 
                            + (1 | Family / Glottocode), 
                       family = "binomial",
                       data = logistic_data)
summary(model_shuffle_null)
AIC(model_shuffle_null)

## Model with random intercepts on length
model_shuffle_no_rd_slope <- glmer(colexification_bin ~ 1 + z_length + z_n_forms + 
                                     (1 | Family / Glottocode ), family="binomial",
                                   data = logistic_data)
summary(model_shuffle_no_rd_slope)
AIC(model_shuffle_no_rd_slope)

## Model with random slopes and intercepts on length 
model_shuffle <- glmer(colexification_bin ~ 1 + z_length + z_n_forms 
                       + (1 + z_length | Family / Glottocode),
     family = "binomial",
     data = logistic_data)
summary(model_shuffle)
AIC(model_shuffle)
# confidence intervals
confint(model_shuffle)
confint(model_shuffle_no_rd_slope)

#1.2. Plots ######

# Plotting the random slopes and intercepts
## Dotplot 
# pdf('figures/random_intercepts.pdf', width=20, height=20)
#print(dotplot(ranef(model_shuffle)))
# dev.off()
random_effects <- ranef(model_shuffle, condVar = TRUE)
random_effects
## Plotting random slopes 
### Create a data frame for plotting
mean_forms <- logistic_data %>%
  group_by(Family) %>%
  summarize(mean_n = mean(n_forms))
plot_data <- data.frame(group = factor(rownames(ranef(model_shuffle)$Family)),
                        slope = ranef(model_shuffle)$Family$z_length)
plot_data <- merge(plot_data, mean_forms, by.x = "group", by.y = "Family", 
                   all.x = TRUE)
documentation = logistic_data %>% group_by(Family) %>% 
  mutate(z_n_forms = (n_forms - mean(n_forms)) / sd(n_forms))
plot_data <- merge(plot_data, documentation, by.x="group", by.y="Family", 
                   all.x = TRUE)
### Plotting
plot_data <- plot_data %>%
  arrange(slope) 
# reorder by slope
plot_data$group <- reorder(plot_data$group, plot_data$slope)
p <- ggplot(plot_data, aes(x = slope, y = group, color = mean_n)) +
  geom_point() +
  labs(y = "", x = "Random Slope (length)", color = "Mean Forms") +
  theme_bw() +
  geom_vline(xintercept = 0, color = 'red', linetype = "dashed") +
  theme(axis.text.y = element_text(size=3)) +
  scale_color_gradientn(colours = rainbow(5))
p
ggsave('figures/random_intercepts.png', p, width = 10, height = 10)

# Avg comparisons
## Computing
avg <- model_shuffle %>%
  avg_predictions(variables = "z_n_forms", 
              by = "z_n_forms", type = "response")
avg
## Plotting avg comparisons
plot_predictions(model_shuffle, 
                 condition = list("z_length", 
                                  "z_n_forms"),
                 type = "response") +
  xlab('Length (z-scored)') +
  ylab('Probability of colexification') +
  # facet_wrap(~condition2, nrow=1) +
  # theme(legend.position = 'none') +
  ylim(0, 1)



#2. Removing outlier languages ######

# Removing outliers by slope 
random_effects_d <- data.frame(random_effects$Family)
random_effects_d <- random_effects_d %>% 
  filter(z_length <= 0.25 & z_length >= -0.25)
## List of families which are in the [-0.25, 0.25] range by slope (length)
family_no <- rownames(random_effects_d)
## Filter families
logistic_data_filtered <- logistic_data %>%
  filter(Family %in% family_no)

# New regression without outliers
model_no_rs_filtered <- glmer(colexification_bin ~ 1 + z_length + z_n_forms 
                          + (1 | Family / Glottocode), 
                          family = "binomial",
                          data = logistic_data_filtered)
summary(model_no_rs_filtered)

#model with random slope
model_rs_filtered <- glmer(colexification_bin ~ 1 + z_length + z_n_forms 
                              + (1 + z_length | Family / Glottocode), 
                              family = "binomial",
                              data = logistic_data_filtered)
summary(model_rs_filtered)

#3. Adding interaction of length*n of forms ####


