library(tidyverse)
library(lme4)

# set ggplot theme
theme_set(theme_bw())

# run the python preprocessing script
# system2('python3 preprocessing.py')

# Inverse logit function
inverse_logit <- function(x){
  exp(x)/(1+exp(x))
}

lexibank_all <- read.csv('data/forms_total.csv')

lexibank_all$Concepticon_ID <- as.numeric(lexibank_all$Concepticon_ID)

lexibank_all_colex <- lexibank_all %>%
  select(Glottocode, Form, 
         Concepticon_ID, Length, 
         Family) %>%
  group_by(Glottocode, Form) %>%
  mutate(colexification = n_distinct(Concepticon_ID)) %>%
  group_by(Glottocode) %>%
  mutate(n_forms = n()) %>%
  distinct(Form, .keep_all = TRUE) %>%
  mutate(colexification_bin = ifelse(colexification > 1, 1, 0))
  

logistic_data <- lexibank_all_colex %>%
  transform(shuffle_colexification_bin = sample(colexification_bin)) %>%
  group_by(Glottocode) %>%
  mutate(z_length = (Length - mean(Length)) / sd(Length))  %>%
  ungroup()  %>%
  mutate(z_n_forms = (n_forms - mean(n_forms)) / sd(n_forms))

model_shuffle_null <- glmer(shuffle_colexification_bin ~ 1 + z_n_forms 
                            + (1 | Family / Glottocode), 
                       family = "binomial",
                       data = logistic_data)
summary(model_shuffle_null)


model_shuffle_no_rd_slope <- glmer(shuffle_colexification_bin ~ 1 + z_length + z_n_forms + 
                                     (1 | Family / Glottocode ), family="binomial",
                                   data = logistic_data)
summary(model_shuffle_no_rd_slope)

model_shuffle <- glmer(shuffle_colexification_bin ~ 1 + z_length + z_n_forms 
                       + (1 + z_length | Family / Glottocode), 
     family = "binomial",
     data = logistic_data)
summary(model_shuffle)
