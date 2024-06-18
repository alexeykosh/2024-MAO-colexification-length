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

# set ggplot theme
theme_set(theme_bw())

# Data processing ---------------------------------------------------------

## Extracting data
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
## Z-scoring
logistic_data <- lexibank_all_colex %>%
  group_by(Glottocode) %>%
  mutate(z_length = (Length - mean(Length)) / sd(Length))  %>%
  ungroup()  %>%
  mutate(z_n_forms = (n_forms - mean(n_forms)) / sd(n_forms))
## Removing empty families
logistic_data <- logistic_data[-which(logistic_data$Family == ""), ]
## Number of entries
print(paste('Number of wordforms:', nrow(logistic_data)))
## Median n. forms
logistic_data %>%
  distinct(Glottocode, .keep_all = TRUE) %>%
  summarise(mean_n_forms = median(n_forms))
## Mean colexifying words
logistic_data %>%
  group_by(Glottocode) %>%
  summarise(mean_colex = mean(colexification_bin)) %>%
  summarise(mean(mean_colex) * 100)
## Median length 
median(logistic_data$Length)


# Preregistered analysis --------------------------------------------------

##Regression models
### Null model without length
model_shuffle_null <- glmer(colexification_bin ~ 1 + z_n_forms 
                            + (1 | Family / Glottocode), 
                       family = "binomial",
                       data = logistic_data)
model_shuffle_null_ <- broom.mixed::tidy(model_shuffle_null, 
                                        effects = "fixed", 
                                        conf.int=TRUE) %>%
  mutate(term = recode(term,
                       "z_n_forms" = "Number of forms"),
         p.value = round(p.value, 4))
model_shuffle_null_
### Model with random intercepts
model_shuffle_no_rd_slope <- glmer(colexification_bin ~ 1 + z_length + z_n_forms + 
                                     (1 | Family / Glottocode ), family="binomial",
                                   data = logistic_data)
round(AIC(model_shuffle_no_rd_slope) - AIC(model_shuffle_null))
model_shuffle_no_rd_slope_ <- broom.mixed::tidy(model_shuffle_no_rd_slope, 
                                                effects = "fixed", 
                                                conf.int=TRUE) %>%
  mutate(p.value = round(p.value, 4))
model_shuffle_no_rd_slope_
### Model with random slopes and intercepts on length 
model_shuffle <- glmer(colexification_bin ~ 1 + z_length + z_n_forms
                       + (1 + z_length | Family / Glottocode),
     family = "binomial",
     data = logistic_data)
round(AIC(model_shuffle)[1] - AIC(model_shuffle_no_rd_slope)[1])
model_shuffle_ <- broom.mixed::tidy(model_shuffle, 
                                    effects = "fixed", 
                                    conf.int=TRUE, 
                                    conf.level=0.95) %>%
  mutate(term = recode(term,
                       "z_length" = "Length",
                       "z_n_forms" = "Number of forms"))
model_shuffle_

## Plotting the random slopes for families
random_effects <- ranef(model_shuffle, condVar = TRUE) 
### Create a data frame for plotting
mean_forms <- logistic_data %>%
  distinct(Glottocode, .keep_all = TRUE) %>%  # Keep only unique combinations of Family and Language
  group_by(Family) %>%
  summarize(mean_n = mean(n_forms))
plot_data <- data.frame(group = factor(rownames(ranef(model_shuffle)$Family)),
                        slope = ranef(model_shuffle)$Family$z_length)
plot_data <- merge(plot_data, mean_forms, by.x = "group", by.y = "Family", 
                   all.x = TRUE)
documentation <- logistic_data %>% group_by(Family) %>% 
  mutate(z_n_forms = (n_forms - mean(n_forms)) / sd(n_forms))
plot_data <- merge(plot_data, documentation, by.x="group", by.y="Family", 
                   all.x = TRUE)
## Plotting random intercepts
plot_data <- plot_data %>%
  arrange(slope) 
### Create color scale
palette_colors <- wes_palette("Zissou1", 100, type = "continuous")
### Calculate the limits for the color scale
min_forms <- min(plot_data$mean_n)
max_forms <- max(plot_data$mean_n)
### Calculate the breaks at the quantiles
quantile_breaks <- quantile(plot_data$mean_n)
### Reorder by slope
plot_data$group <- reorder(plot_data$group, plot_data$slope)
### Plot random slopes
p <- ggplot(plot_data, aes(x = slope, y = group, color = mean_n)) +
  geom_point() +
  labs(y = "", x = "β-coefficient for length (deviation from the fixed effect)", 
       color = "Mean Forms") +
  theme_bw() +
  geom_vline(xintercept = 0, color = 'red', linetype = "dashed") +
  theme(axis.text.y = element_text(size=6)) +
  scale_color_gradientn(
    colours = palette_colors,  # Use the chosen wesanderson palette
    trans = "log10",  # Apply log10 transformation
    name = "Number of forms",  # Title for the colorbar
    limits = c(min_forms, max_forms),  # Set limits for the color scale
    breaks = quantile_breaks,  # Set the breaks at quantiles
    labels = round(quantile_breaks, 2)  # Optionally, round the break labels
  )
### Save figure in png
ggsave('figures/word/random_intercepts.png', p, width = 10, height = 10)
### Save figure in pdf
ggsave('figures/vector/random_intercepts.pdf', p, width = 10, height = 10)

## Avg comparisons
### Calculate mean and standard deviation of n_forms
sd_n_forms <- sd(logistic_data$n_forms)
mean_n_forms <- mean(logistic_data$n_forms)
### Function to reverse z-scoring
unzscore_n_forms <- function(z) {
  return(z * sd_n_forms + mean_n_forms)
}
### Calculate the quantiles of z_n_forms
quantiles_z_n_forms <- quantile(logistic_data$z_n_forms)
### Reverse z-scoring transformation for each quantile
unzscored_quantiles <- sapply(quantiles_z_n_forms, unzscore_n_forms)
### Plot predictions
p <- plot_predictions(model_shuffle, 
                      condition = list("z_length", 
                                       "z_n_forms" = quantiles_z_n_forms),
                      type = "response")
### Extract plot data and compute y_center
plot_data_ <- ggplot_build(p)$data[[1]] %>% 
  rowwise() %>% 
  mutate(y_center = mean(c(ymin, ymax))) 
### Create a mapping from groups to unzscored quantiles
group_quantile_map <- setNames(round(unzscored_quantiles), 
                               seq_along(unzscored_quantiles))
### Assign the unzscored quantiles to the plot data based on groups
plot_data_ <- plot_data_ %>% mutate(z_n_forms = 
                                      group_quantile_map[as.character(group)])
### Create the relation plot
relation_plot <- ggplot(plot_data_, aes(x = x, y = y_center, 
                                       color = factor(z_n_forms))) +
  geom_line(linewidth=1) +
  xlab('Length (z-score)') +
  ylab('Probability of colexification') +
  scale_color_manual(name = "Number of Forms\n(rounded quantile):",
                     values = wes_palette("Zissou1", type = "discrete")) +
  scale_fill_manual(values = wes_palette("Zissou1", type = "discrete")) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax, 
                  fill = factor(z_n_forms)), 
              alpha = 0.5, linetype='blank') +
  ylim(0, 1) +
  guides(fill="none") +
  theme(legend.key = element_blank(),
        legend.box.background = element_blank(),
        legend.key.size = unit(1, "cm")) + 
  guides(color = guide_legend(override.aes = list(fill = NA))) 
### Save for word
ggsave('figures/word/quantile_control.png', relation_plot, 
       width = 10, height = 5)
### Save in pdf
ggsave('figures/vector/quantile_control.pdf', relation_plot, 
       width = 10, height = 5)

## Plotting avg length against slope:
mean_forms <- logistic_data %>% # Keep only unique combinations of Family and Language
  group_by(Glottocode) %>%
  summarize(mean_length = mean(Length))
random_effects <- ranef(model_shuffle, condVar = TRUE)
random_effects_l <- random_effects$`Glottocode:Family`
random_effects_l$z_length <- fixef(model_shuffle)[2] +random_effects_l$z_length
random_effects_l <- cbind(language = rownames(random_effects_l), 
                          random_effects_l)
row.names(random_effects_l) <- NULL
random_effects_l <- random_effects_l %>%
  separate(language, into = c("Glottocode", "family"), sep = ":")
plot_data_sl <- merge(random_effects_l, mean_forms, by="Glottocode", 
                   all.x = TRUE)
### Create the plot
avg_length <- ggplot(data = plot_data_sl, aes(x = mean_length, y = z_length)) +
  geom_point(alpha = 0.8) +  # Default color  # Linear model smooth line without confidence interval
  xlab('Average Length') +  # Custom x-axis label
  ylab('Slope') +  # Custom y-axis label  # Use a minimal theme for a cleaner look
  stat_cor(method = "pearson", p.accuracy = 0.001, r.accuracy = 0.01) 
### Save for word
ggsave('figures/word/mean_length_language.png', avg_length, 
       width = 10, height = 5)
### Save in pdf
ggsave('figures/vector/mean_length_language.pdf', avg_length, 
       width = 10, height = 5)


## Beta coefficients distribution
### Outliers
coef_h_0 <- random_effects_l[random_effects_l$z_length >= 0,]
### Plot figure
coef_dist <- ggplot(random_effects_l, aes(x = z_length)) +
  geom_dotplot(method = 'histodot', binwidth = 0.0105, alpha=1) +
  geom_text_repel(data=coef_h_0, aes(x = z_length, y = 0, label = Glottocode), 
                  nudge_y = 0.2, 
                  size = 6,
                  color = 'blue') +
  labs(x = expression(beta * "-coefficient for length"),
       y = 'Density')
### Save for word
ggsave('figures/word/beta_distr.png', coef_dist, 
       width = 10, height = 5)
### Save in pdf
ggsave('figures/vector/beta_distr.pdf', coef_dist, 
       width = 10, height = 5)

# Plotting descriptive maps -----------------------------------------------

## Map with random slopes
### Get slopes values
random_effects <- ranef(model_shuffle, condVar = TRUE)
random_effects_l <- random_effects$`Glottocode:Family`
random_effects_l$z_length <- fixef(model_shuffle)[2] +random_effects_l$z_length
random_effects_l <- cbind(language = rownames(random_effects_l), 
                          random_effects_l)
row.names(random_effects_l) <- NULL
random_effects_l <- random_effects_l %>%
  separate(language, into = c("language", "family"), sep = ":")
random_effects_l <- random_effects_l %>% rename(Glottocode = language)
### Get latitudes and longitudes
languages_info <- read.csv('data/lexibank-lexibank-analysed-a4c0952/cldf/languages.csv') %>%
  distinct(Glottocode, .keep_all = TRUE)
### Merge the datasets
languages_info <- languages_info %>%
  inner_join(random_effects_l, by = "Glottocode")
### Number of forms
forms_n <- logistic_data %>%
  distinct(Glottocode, .keep_all = T) %>%
  select(Glottocode, n_forms)
languages_info <- languages_info %>%
  inner_join(forms_n, by = "Glottocode")
### Map data
world <- map_data("world") %>% 
  filter(lat > -62)
### Identify the top-10 lowest z_length values
top10_lowest_z_length <- languages_info %>%
  arrange(z_length) %>%
  head(10)
### Create the plot with ggrepel
full_map_results <- ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group), 
               colour = NA, fill = 'grey', alpha = 0.4) +
  geom_point(data = languages_info, aes(x = Longitude, y = Latitude, 
                                        color = z_length, size = n_forms)) +
  scale_color_gradientn(
    colors = c('#c55a01', 'white', '#0097c3'), 
    values = scales::rescale(c(min(languages_info$z_length), 
                               0, 
                               max(languages_info$z_length))),
    guide = guide_colorbar(title = expression(beta * "-coefficient for length:"))
  ) +
  scale_size_continuous(
    name = "Number of forms:"
  ) +
  geom_text_repel(data = top10_lowest_z_length, aes(x = Longitude, 
                                                    y = Latitude, label = Name),
                  size = 4, 
                  nudge_y = 10, 
                  segment.size = 0.2) +
  theme_void() +
  theme(
    legend.position = 'bottom',
    legend.direction = "horizontal",
    legend.title = element_text(vjust = 0.8)  # Center align the legend title
  )
### Save for word
ggsave('figures/word/map_results.png', full_map_results, width = 10, height = 5)
### Save in pdf
ggsave('figures/vector/map_results.pdf', full_map_results, width = 10, height = 5)

## Map with number of forms
###  Create the palette colors using wes_palette
palette_colors <- wes_palette("Zissou1", 100, type = "continuous")
###  Calculate the limits for the color scale
min_forms <- min(languages_info$n_forms)
max_forms <- max(languages_info$n_forms)
###  Calculate the breaks at the quantiles
quantile_breaks <- quantile(logistic_data$n_forms)
###  Create the full basic map
full_basic_map <- ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group), 
               colour = NA, fill = 'grey', alpha = 0.4) +
  geom_point(data = languages_info, aes(x = Longitude, y = Latitude, 
                                        color = n_forms)) +
  scale_color_gradientn(
    colours = palette_colors,  # Use the chosen wesanderson palette
    trans = "log10",  # Apply log10 transformation
    name = "Number of forms",  # Title for the colorbar
    limits = c(min_forms, max_forms),  # Set limits for the color scale
    breaks = quantile_breaks,  # Set the breaks at quantiles
    labels = round(quantile_breaks, 2)  # Optionally, round the break labels
  ) +
  theme_void() +
  theme(
    legend.position = 'bottom',
    legend.direction = "horizontal",
    legend.title = element_text(vjust = 0.8),  # Center align the legend title
    legend.key.width = unit(2, 'cm')  # Increase the width of the colorbar
  )
### Save for word
ggsave('figures/word/basic_map.png', full_basic_map, width = 10, height = 5)
### Save in pdf
ggsave('figures/vector/basic_map.pdf', full_basic_map, width = 10, height = 5)
