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
  # transform(shuffle_colexification_bin = sample(colexification_bin)) %>%
  group_by(Glottocode) %>%
  mutate(z_length = (Length - mean(Length)) / sd(Length))  %>%
  ungroup()  %>%
  mutate(z_n_forms = (n_forms - mean(n_forms)) / sd(n_forms))


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
## Plotting the random slopes and intercepts
# pdf('figures/random_intercepts.pdf', width=20, height=20)
# print(dotplot(ranef(model_shuffle)))
# dev.off()
### Create a data frame for plotting
mean_forms <- logistic_data %>%
  group_by(Family) %>%
  summarize(mean_n = mean(n_forms))
plot_data <- data.frame(group = factor(rownames(ranef(model_shuffle)$Family)),
                        intercept = ranef(model_shuffle)$Family$z_length)
plot_data <- merge(plot_data, mean_forms, by.x = "group", by.y = "Family", all.x = TRUE)
### Plot using ggplot2
plot_data <- plot_data %>%
  arrange(intercept) %>%
  mutate(group = factor(group, levels = group))
p <- ggplot(plot_data, aes(x = intercept, y = group, color = mean_n)) +
  geom_point() +
  labs(y = "", x = "Random Intercept", color = "Mean Forms") +
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
 