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


#j'ai appelé les 3 modèles model_null, model_no_rd_slope, model
#et voilà les analyses que j'ai faites même si tu m'as dit qu'avec geom_smooth ça n'allait pas

library(ggplot2)

ggplot(logistic_data, aes(y = colexification_bin, x=z_length)) +  
  geom_point(size=0.2) + 
  geom_smooth(method="glm", method.args = list(family = binomial))  +
  labs(
    title = "Probability of colexification according to wordform length all data", 
    x = "Length (z-scored)",
    y = "Probability of colexification"
  )

quantiles <- quantile(logistic_data$z_n_forms, probs = seq(0, 1, by = 0.20))
logistic_data$z_n_forms_category <- cut(logistic_data$z_n_forms, breaks = quantiles, include.lowest = TRUE)
logistic_data$z_n_forms_category <- factor(logistic_data$z_n_forms_category)
levels(logistic_data$z_n_forms_category)

ggplot(logistic_data, aes(y = colexification_bin, x = z_length,color=z_n_forms_category)) +  
  geom_point(size = 0.2) + 
  geom_smooth(method = "glm", method.args = list(family = binomial)) +
  labs(
    title = "Probability of colexification according to wordform length all data", 
    x = "Length (z-scored)",
    y = "Probability of colexification"
  ) 

residuals <- resid(model)

ggplot(data.frame(residuals = residuals), aes(x = seq_along(residuals), y = residuals)) +
  geom_point(size=0.1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Observation", y = "Residuals", title = "Plot of Residuals")

residus_simules <- simulateResiduals(fittedModel = model)
plot(residus_simules, type = "qqplot")
resultat_test <- testOutliers(simulationOutput = residus_simules, type = 'bootstrap')
print(resultat_test)
testOutliers(residus_simules, alternative = c("two.sided", "greater",
                                              "less"), margin = c("both", "upper", "lower"), type = c("default",
                                                                                                      "bootstrap", "binomial"), nBoot = 100, plot = T)

install.packages("sjPlot")
install.packages("glmmTMB")
library(sjPlot)
plot_model(model, type = "re", vspacing = 10, vsize = 0.02,)

install.packages("boot")
library(boot)
help(boot)

residuals(model, type = "deviance")

data_S_T<- logistic_data %>% 
  filter(Family == "Sino-Tibetan")

ggplot(data_S_T, aes(y = colexification_bin, x = z_length)) +  
  geom_point(size = 0.2) + 
  geom_smooth(method = "glm", method.args = list(family = binomial)) +
  labs(
    title ="Probability of colexification according to wordform length for Sino_Tibetan families", 
    x = "Length (z-scored)",
    y = "Probability of colexification"
  )

data_I_E<- logistic_data %>% 
  filter(Family == "Indo-European")

ggplot(data_I_E, aes(y = colexification_bin, x = z_length)) +  
  geom_point(size = 0.2) + 
  geom_smooth(method = "glm", method.args = list(family = binomial)) +
  labs(
    title ="Probability of colexification according to wordform length for Indo-European families", 
    x = "Length (z-scored)",
    y = "Probability of colexification"
  )

data_Ijoid<- logistic_data %>% 
  filter(Family == "Ijoid")

ggplot(data_Ijoid, aes(y = colexification_bin, x = z_length)) +  
  geom_point(size = 0.2) + 
  geom_smooth(method = "glm", method.args = list(family = binomial)) +
  labs(
    title ="Probability of colexification according to wordform length for Ijoid families", 
    x = "Length (z-scored)",
    y = "Probability of colexification"
  )

plot(model)
 