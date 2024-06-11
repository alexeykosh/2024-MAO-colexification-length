library(tidyverse)
library(lme4)
library(marginaleffects)
library(DHARMa)
library(lattice)

# set ggplot theme
theme_set(theme_bw())

# run the python preprocessing script
# system2('python3 preprocessing.py')

# Inverse logit function
inverse_logit <- function(x){
  exp(x)/(1+exp(x))
}

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


# Plotting the random slopes and intercepts
## Dotplot 
# pdf('figures/random_intercepts.pdf', width=20, height=20)
# print(dotplot(ranef(model_shuffle)))
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
### Ploting
plot_data <- plot_data %>%
  arrange(slope) 
p <- ggplot(plot_data, aes(x = slope, y = group, color = mean_n)) +
  geom_point() +
  labs(y = "", x = "Random Slope (length)", color = "Mean Forms") +
  theme_bw() +
  geom_vline(xintercept = 0, color = 'red', linetype = "dashed") +
  theme(axis.text.y = element_text(size=6)) +
  scale_color_gradientn(colours = rainbow(5))
ggsave('figures/random_intercepts.pdf', p, width = 10, height = 10)

# Avg comparisons
## Computing
avg <- model_shuffle %>%
  avg_predictions(variables = "z_n_forms", 
              by = "z_n_forms", type = "response")

## Plotting avg comparisons
avg %>%
  ggplot(aes(x = z_n_forms, y=estimate))+
  geom_point()+
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high))+
  xlab('Number of forms') +
  ylab('Slope')

## Plotting avg length against slope:
plot_data %>%
  group_by(group) %>%
  summarize(mean_slope = max(slope), mean_length = mean(Length)) %>%
  ggplot(aes(x=mean_length, y=mean_slope)) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, 
                ymin = -Inf, ymax = 0), fill = "grey80", alpha = 0.1) +
  geom_point() + # annotate with geom_repel
  xlab('Average length') +
  ylab('Slope')
  

# plot_predictions(model_shuffle, 
#                  condition = list("z_length", 
#                                   "z_n_forms"),
#                  type = "response") +
#   xlab('Length (z-scored)') +
#   ylab('Probability of colexification') +
#   # facet_wrap(~condition2, nrow=1) +
#   # theme(legend.position = 'none') +
#   ylim(0, 1)


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



######## Test with small dataset (only 4 languages, 3 families) Marie #######
small_dataset <- logistic_data %>%
  filter(Glottocode %in% c("sart1249","wang1301","kolu1245","swed1254"))

# randomly sample 40 languages:

# do the subset with the 40 random lgs: 


s_model_shuffle_null <- glmer(colexification_bin ~ 1 + z_n_forms 
                            + (1 | Family / Glottocode), 
                            family = "binomial",
                            data = small_dataset)
summary(s_model_shuffle_null)
AIC(s_model_shuffle_null)
## Model with random intercepts on length
s_model_shuffle_no_rd_slope <- glmer(colexification_bin ~ 1 + z_length + z_n_forms + 
                                     (1 | Family / Glottocode ), family="binomial",
                                   data = small_dataset)
summary(s_model_shuffle_no_rd_slope)
AIC(s_model_shuffle_no_rd_slope)
## Model with random slopes and intercepts on length 
s_model_shuffle <- glmer(colexification_bin ~ 1 + z_length + z_n_forms 
                       + (1 + z_length | Family / Glottocode), 
                       family = "binomial",
                       data = small_dataset)
summary(s_model_shuffle)
AIC(s_model_shuffle)
## Plotting the random slopes and intercepts
# pdf('figures/random_intercepts.pdf', width=20, height=20)
# print(dotplot(ranef(model_shuffle)))
# dev.off()

random_effects <- ranef(s_model_shuffle, condVar = TRUE)
random_effects

### Create a data frame for plotting
s_mean_forms <- small_dataset %>%
  group_by(Family) %>%
  summarize(mean_n = mean(n_forms))
s_plot_data <- data.frame(group = factor(rownames(ranef(s_model_shuffle)$Family)),
                        intercept = ranef(s_model_shuffle)$Family$z_length)
s_plot_data <- merge(s_plot_data, s_mean_forms, by.x = "group", by.y = "Family", all.x = TRUE)
s_documentation = small_dataset %>% group_by(Family) %>% mutate(z_n_forms = (n_forms - mean(n_forms)) / sd(n_forms))
s_plot_data <- merge(s_plot_data, s_documentation, by.x="group", by.y="Family", all.x = TRUE)

### Plot using ggplot2
s_plot_data <- s_plot_data %>%
  arrange(intercept) %>%
  mutate(group = factor(group, levels = group))
ggplot(s_plot_data, aes(x = intercept, y = group, color = mean_n)) +
  geom_point() +
  labs(y = "", x = "Random Slope", color = "Mean Forms") +
  theme_bw() +
  geom_vline(xintercept = 0, color = 'red', linetype = "dashed") +
  theme(axis.text.y = element_text(size=6)) +
  scale_color_gradientn(colours = rainbow(5))
# s
# ggsave('figures/random_intercepts.pdf', p, width = 10, height = 10)

s_plot_data <- s_plot_data %>%
  arrange(intercept)

ggplot(s_plot_data, aes(x = intercept, y = group, color = z_n_forms)) +
  geom_point() +
  labs(y = "", x = "Random Slope", color = "Level of documentation") +
  theme_bw() +
  geom_vline(xintercept = 0, color = 'red', linetype = "dashed") +
  theme(axis.text.y = element_text(size=6)) +
  scale_color_gradientn(colours = rainbow(5))


s_intercept <- ranef(s_model_shuffle)$Family[, "(Intercept)"]
s_median_intercept <- median(s_intercept)
s_sd_intercept <- sd(s_intercept)
s_outliers <- levels(small_dataset$Family)[abs(s_intercept - s_median_intercept) > 0.25 * s_sd_intercept]
s_data_filtered <- small_dataset %>%
  filter(!(Family %in% s_outliers))
s_model_filtered <- glmer(colexification_bin ~ 1 + z_length + z_n_forms + (1 + z_length | Family / Glottocode), 
                        family = "binomial",
                        data = s_data_filtered)
s_model_filtered


# Avg comparisons
## Computing
s_avg <- s_model_shuffle %>%
  avg_predictions(variables = "z_n_forms", 
                  by = "z_n_forms", type = "response")
s_avg
## Plotting avg comparisons
plot_predictions(s_model_shuffle, 
                 condition = list("z_length", 
                                  "z_n_forms"),
                 type = "response") +
  xlab('Length (z-scored)') +
  ylab('Probability of colexification') +
  # facet_wrap(~condition2, nrow=1) +
  # theme(legend.position = 'none') +
  ylim(0, 1)

# Outlier removal and new analysis
## Getting the intercept
intercept <- ranef(model_shuffle)$Family[, "(Intercept)"]
# computing the median
median_intercept <- median(intercept)
sd_intercept <- sd(intercept)
outliers <- levels(logistic_data$Family)[abs(intercept - median_intercept) > 
                                           0.25 * sd_intercept]
## Filter the outlies from the initial dataset
logistic_data_filtered <- logistic_data %>%
  filter(!(Family %in% outliers))
## Fit a new model 
model_filtered <- glmer(colexification_bin ~ 1 + z_length + z_n_forms + (1 + z_length | Family / Glottocode), 
                        family = "binomial",
                        data = logistic_data_filtered)
model_filtered


###

## Summary table
type_of_model=c("null", "no random slope", "random slope", "no random slope")
outliers_present=c("yes", "yes", "yes", "no")
AIC=c(AIC(model_shuffle_null), AIC(model_shuffle_no_rd_slope),
      AIC(model_shuffle), AIC(model_no_rs_filtered))

df_summary= tibble(type_of_model, outliers_present, AIC)
df_summary
