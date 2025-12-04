#This code is based on Alexey's code, downloaded 2025-08-08 from his github,
#changed today by Olivier. 
#The data used is the file combined_ambiguity_freq.csv, sent to me by Alexey on 2025-07-02
#The modifications were all discussed on workflowy in July, MAO thread.


library(tidyverse)
library(lme4)
library(marginaleffects)
library(DHARMa)
library(lattice)
library(wesanderson)
library(broom.mixed)
library(ggpubr)
library(ggrepel)
library(ggridges)
library(grid)
library(scales)

# set ggplot theme
theme_set(theme_bw())
# colors 
palette_colors <- wes_palette("Zissou1", 100, 
                              type = "continuous")


# Data processing -----------------------------------------------------------

## Extracting and preprocessing data
all_dat <- read.csv('data/combined_ambiguity_freq.csv') %>%
  group_by(lang) %>%
  mutate(
    z_freq = (log(freq) - mean(log(freq))) / sd(log(freq)),
    z_bert_polysemy = (bert_polysemy - mean(bert_polysemy)) / sd(bert_polysemy),
    z_length = (word_length - mean(word_length)) / sd(word_length)
  )


# Analysis ----------------------------------------------------------------


#Analysis: BERT polysemy #####

## Null model (just frequency)
model_null <- lmer(z_bert_polysemy ~ 1 + z_freq 
                            + (1 | language_name), 
                            # family = "gaussian",
                            data = all_dat)
summary(model_null)
model_null_ <- broom.mixed::tidy(model_null, 
                                         effects = "fixed", 
                                         conf.int=TRUE) %>%
  mutate(term = recode(term,
                       "z_freq" = "Relative frequency(log)")) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))
model_null_

## Full model frequency and length with a random intercept and slope
model_full_BERT <- lmer(z_bert_polysemy ~ 1 + z_freq + z_length
                   + (1 + z_length| language_name), 
                   # family = "gaussian",
                   data = all_dat) 
summary(model_full_BERT)
model_full_BERT_ <- broom.mixed::tidy(model_full_BERT, 
                                 effects = "fixed", 
                                 conf.int=TRUE) %>%
  mutate(term = recode(term,
                       "z_freq" = "Relative frequency (log)",
                       "z_length" = "Length")) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))
model_full_BERT_

## Model comparison
round(AIC(model_full_BERT) - AIC(model_null))

## Partial correlation
ppcor::pcor.test(all_dat$z_length, all_dat$z_bert_polysemy, 
                 all_dat$z_freq, method='spearman')

languages <- unique(all_dat$lang)

#Looking at the raw correlation, language by language:
for(i in languages){
  print(i)
  set <- subset(all_dat, all_dat$lang == i)
  print(cor(set$z_bert_polysemy, set$z_length, method = "spearman"))
}


#Analysis: WORDNET polysemy #####

#Removing the languages for which there is no data on Wordnet polysemy 
#(i.e. wordnet polysemy = -1) 

raw_data <- read.csv('data/combined_ambiguity_freq.csv') 
wn_data <-  subset(raw_data, raw_data$wordnet_polysemy > -1)
#That's 7 languages:
unique(wn_data$lang)
#But Hebrew only has 0 values for WN polysemy:
unique(wn_data[wn_data$lang == "he",]$wordnet_polysemy)

#and so we remove it too: 
wn_data <-  subset(raw_data, raw_data$wordnet_polysemy > -1)
wn_data <-  wn_data[wn_data$lang != "he", ]

#Now 6 languages:
unique(wn_data$lang)

#z-scoring:
wn_data <- wn_data %>%
  group_by(lang) %>%
  mutate(
    z_freq = (log(freq) - mean(log(freq))) / sd(log(freq)),
    z_wordnet_polysemy = (wordnet_polysemy - mean(wordnet_polysemy)) / sd(wordnet_polysemy),
    z_length = (word_length - mean(word_length)) / sd(word_length)
  )

## Null model (just frequency), for these 7 languages:
model_null_wn <- lmer(z_wordnet_polysemy ~ 1 + z_freq 
                   + (1 | language_name), 
                   # family = "gaussian",
                   data = wn_data)
summary(model_null_wn)
model_null_wn_ <- broom.mixed::tidy(model_null_wn, 
                                 effects = "fixed", 
                                 conf.int=TRUE) %>%
  mutate(term = recode(term,
                       "z_freq" = "Relative frequency(log)"))  %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))
model_null_wn_


model_full_WORDNET <- lmer(z_wordnet_polysemy ~ 1 + z_freq + z_length
                        + (1 + z_length| language_name), 
                        # family = "gaussian",
                        data = wn_data)
summary(model_full_WORDNET)


model_full_WORDNET_ <- broom.mixed::tidy(model_full_WORDNET, 
                                      effects = "fixed", 
                                      conf.int=TRUE) %>%
  mutate(term = recode(term,
                       "z_freq" = "Relative frequency (log)",
                       "z_length" = "Length"))   %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))
model_full_WORDNET_

## Model comparison
round(AIC(model_full_WORDNET) - AIC(model_null))

## Partial correlation
v <- ppcor::pcor.test(wn_data$z_length, wn_data$wordnet_polysemy, 
                 wn_data$z_freq, method='spearman')


#Looking at the raw correlation, language by language:

languages <- unique(wn_data$lang)
for(i in languages){
  print(i)
  set <- subset(wn_data, wn_data$lang == i)
print(cor(set$z_wordnet_polysemy, set$z_length, method = "spearman"))
}

# Plotting ----------------------------------------------------------------

## Create a prediction grid: different z_length sequence for each language
newdat <- all_dat %>%
  group_by(language_name) %>%
  summarise(min_z = min(z_length, na.rm = TRUE),
            max_z = max(z_length, na.rm = TRUE)) %>%
  rowwise() %>%
  do({
    data.frame(
      z_length = seq(.$min_z, .$max_z, length.out = 100),
      language_name = .$language_name,
      z_freq = 0  # keep z_freq fixed at mean
    )
  }) %>%
  ungroup()
## Compute predictions
pred <- predictions(model_full_BERT, newdata = newdat)
## Plot with model predictions
ggplot() +
  geom_point(data = all_dat, 
             aes(x = z_length, 
                 y = z_bert_polysemy, 
                 color = z_freq), 
             alpha = 0.2) +
  geom_line(data = pred, 
            aes(x = z_length, 
                y = estimate, 
                group = language_name),  
            color = "black",
            linewidth = 1.2) +           
  # scale_color_gradient(low = "blue", high = "red") +
  scale_color_gradientn(
    colours = palette_colors) +
  labs(
    x = "Word length (z-scored)",
    y = "Estimated degree of polysemy (z-scored)",
    color = "Word Frequency (z-scored)\n\u2190 Low                          High \u2192"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.key.width = unit(1, "cm"),
    legend.key.height = unit(0.4, "cm"),
    legend.box = "vertical",
    legend.title = element_text(hjust = 0.5, size = 10),
    legend.text = element_text(size = 8),
    strip.text = element_text(size = 10),
    plot.margin = margin(20, 40, 20, 40) # Add extra margin space for arrows
  ) +
  guides(
    color = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5
    )
  ) +
  facet_wrap(~language_name, scales = 'free', nrow = 3)
## Save as png
ggsave('figures/word/length_degree.png', width = 10, height = 5)
## Save as pdf
ggsave('figures/vector/length_degree.pdf', width = 10, height = 5)

## Create a prediction grid for WordNet data: different z_length sequence for each language
newdat_wn <- wn_data %>%
  group_by(language_name) %>%
  summarise(min_z = min(z_length, na.rm = TRUE),
            max_z = max(z_length, na.rm = TRUE)) %>%
  rowwise() %>%
  do({
    data.frame(
      z_length = seq(.$min_z, .$max_z, length.out = 100),
      language_name = .$language_name,
      z_freq = 0  # keep z_freq fixed at mean
    )
  }) %>%
  ungroup()
## Compute predictions
pred_wn <- predictions(model_full_WORDNET, newdata = newdat_wn)
## Plot with model predictions
ggplot() +
  geom_point(data = wn_data, 
             aes(x = z_length, 
                 y = z_wordnet_polysemy, 
                 color = z_freq), 
             alpha = 0.2) +
  geom_line(data = pred_wn, 
            aes(x = z_length, 
                y = estimate, 
                group = language_name),  
            color = "black",
            linewidth = 1.2) +           
  # scale_color_gradient(low = "blue", high = "red") +
  scale_color_gradientn(
    colours = palette_colors) +
  labs(
    x = "Word length (z-scored)",
    y = "Estimated degree of polysemy \n (WordNet, z-scored)",
    color = "Word Frequency (z-scored)\n\u2190 Low                          High \u2192"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.key.width = unit(1, "cm"),
    legend.key.height = unit(0.4, "cm"),
    legend.box = "vertical",
    legend.title = element_text(hjust = 0.5, size = 10),
    legend.text = element_text(size = 8),
    strip.text = element_text(size = 10),
    plot.margin = margin(20, 40, 20, 40)
  ) +
  guides(
    color = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5
    )
  ) +
  facet_wrap(~language_name, scales = 'free', nrow = 2)
## Save as png
ggsave('figures/word/length_degree_wordnet.png', width = 10, height = 4)
## Save as pdf
ggsave('figures/vector/length_degree_wordnet.pdf', width = 10, height = 4)
