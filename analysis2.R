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

# set ggplot theme
theme_set(theme_bw())


# Data processing -----------------------------------------------------------

## Extracting and preprocessing data
all_dat <- read.csv('data/combined_ambiguity_freq.csv') %>%
  mutate(
    z_freq = (log(freq) - mean(log(freq))) / sd(log(freq)),
    z_bert_polysemy = (bert_polysemy - mean(bert_polysemy)) / sd(bert_polysemy),
    z_surprisal = (surprisal - mean(surprisal)) / sd(surprisal),
    z_length = (word_length - mean(word_length)) / sd(word_length)
  )


# Analysis ----------------------------------------------------------------

## Null model (just frequency)
model_null <- lmer(z_bert_polysemy ~ 1 + z_freq 
                            + (1 | language_name), 
                            # family = "gaussian",
                            data = all_dat)
model_null_ <- broom.mixed::tidy(model_null, 
                                         effects = "fixed", 
                                         conf.int=TRUE) %>%
  mutate(term = recode(term,
                       "z_freq" = "Relative frequency(log)"))
model_null_
## Full model frequency and length with a random intercept and slope
model_full <- lmer(z_bert_polysemy ~ 1 + z_freq + z_length
                   + (1 + z_length| language_name), 
                   # family = "gaussian",
                   data = all_dat)
model_full_ <- broom.mixed::tidy(model_full, 
                                 effects = "fixed", 
                                 conf.int=TRUE) %>%
  mutate(term = recode(term,
                       "z_freq" = "Relative frequency (log)",
                       "z_length" = "Length"))
model_full_
## Model comparison
round(AIC(model_full) - AIC(model_null))


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
pred <- predictions(model_full, newdata = newdat)
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
  scale_color_gradient(low = "blue", high = "red") +
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
