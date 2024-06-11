library(tidyverse)
library(lme4)
library(marginaleffects)
library(DHARMa)
library(lattice)
library(wesanderson)
library(broom.mixed)
library(ggpubr)
library(ggrepel)

# set ggplot theme
theme_set(theme_bw())

#0. Data processing #####

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
## Removing empty families
logistic_data <- logistic_data[-which(logistic_data$Family == ""), ]

yoruba_data <- logistic_data[logistic_data$Glottocode == 'yoru1245',]

#1. Basic analysis, as preregistered ######

#1.1. Regression models ######

## Null model without length
model_shuffle_null <- glmer(colexification_bin ~ 1 + z_n_forms 
                            + (1 | Family / Glottocode), 
                       family = "binomial",
                       data = logistic_data)
AIC(model_shuffle_null)
model_shuffle_null_ <- broom.mixed::tidy(model_shuffle_null, 
                                        effects = "fixed", 
                                        conf.int=TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(term = recode(term,
                       "z_n_forms" = "Number of forms"))

## Model with random intercepts
model_shuffle_no_rd_slope <- glmer(colexification_bin ~ 1 + z_length + z_n_forms + 
                                     (1 | Family / Glottocode ), family="binomial",
                                   data = logistic_data)
AIC(model_shuffle_no_rd_slope)
model_shuffle_no_rd_slope_ <- broom.mixed::tidy(model_shuffle_no_rd_slope, 
                                                effects = "fixed", 
                                                conf.int=TRUE)


## Model with random slopes and intercepts on length 
model_shuffle <- glmer(colexification_bin ~ 1 + z_length + z_n_forms
                       + (1 + z_length | Family / Glottocode),
     family = "binomial",
     data = logistic_data)
AIC(model_shuffle)
model_shuffle_ <- broom.mixed::tidy(model_shuffle, 
                                    effects = "fixed", 
                                    conf.int=TRUE, 
                                    conf.level=0.95) %>%
  filter(term != "(Intercept)") %>%
  mutate(term = recode(term,
                       "z_length" = "Length",
                       "z_n_forms" = "Number of forms"))


#1.2. Plots ######

## Plotting regression coefficients
### Plot the coefficients of model w/o length
p1 <- ggplot(data=model_shuffle_null_, aes(x = estimate, y = term)) +
  geom_point(size=5) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "", y = "") +
  xlim(-0.75, 0.75)
### Plot the coefficients of model w/ length 
p2 <- ggplot(data=model_shuffle_, aes(x = estimate, y = term)) +
  geom_point(size=5) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Estimate", y = "") +
  xlim(-0.75, 0.75)
### Combine the two figures
ggarrange(p1, p2, ncol=1, nrow=2, 
          heights = c(1, 2))

## Plotting the random slopes and intercepts
random_effects <- ranef(model_shuffle, condVar = TRUE) 
# mean(fixef(model_shuffle)[2] + random_effects$`Glottocode:Family`$z_length)
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
  labs(y = "", x = "β-coefficient for length (deviation from the mean)", color = "Mean Forms") +
  theme_bw() +
  geom_vline(xintercept = 0, color = 'red', linetype = "dashed") +
  theme(axis.text.y = element_text(size=6)) +
  scale_color_gradientn(colours = rainbow(5))
ggsave('figures/random_intercepts.png', p, width = 10, height = 10)

# Avg comparisons

# Example mean and standard deviation values
mean_n_forms <- mean(logistic_data$n_forms)
sd_n_forms <- sd(logistic_data$n_forms)

unzscore_n_forms <- function(z) {
  return(z * sd_n_forms + mean_n_forms)
}


p <- plot_predictions(model_shuffle, 
                      condition = list("z_length", 
                                       "z_n_forms"),
                      type = "response")

# Extract the plot data
plot_ <- ggplot_build(p)$data[[1]] %>% 
  rowwise() %>% 
  mutate(y_center = mean(c(ymin, ymax))) 

plot_$z_n_forms <- round(unzscore_n_forms(plot_data$group))

relation_plot <- ggplot(plot_, aes(x = x, y = y_center, 
                                       color = factor(z_n_forms))) +
  geom_line(linewidth=1) +
  xlab('Length (z-score)') +
  ylab('Probability of colexification') +
  scale_color_manual(name = " Number of Forms\n (rounded quantile):",
                     values = wes_palette("Zissou1", type = "discrete")) +
  scale_fill_manual(values = wes_palette("Zissou1", type = "discrete"))+
  geom_ribbon(aes(ymin = ymin, ymax = ymax, 
                  fill = factor(z_n_forms)), 
              alpha = 0.5, linetype='blank') +
  ylim(0, 1) +
  guides(fill="none") +
  theme(legend.key = element_blank(),
        legend.box.background = element_blank(),
        legend.key.size = unit(1, "cm")) +  # Adjust the size of the legend keys if necessary
  guides(color = guide_legend(override.aes = list(fill = NA)))  # This line removes the background

relation_plot

ggsave('figures/quantile_control.png', relation_plot, 
       width = 10, height = 5)

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
AIC(model_no_rs_filtered)

#model with random slope and intercept
model_rs_filtered <- glmer(colexification_bin ~ 1 + z_length + z_n_forms 
                              + (1 + z_length | Family / Glottocode), 
                              family = "binomial",
                              data = logistic_data_filtered)
summary(model_rs_filtered)

#3. Adding interaction of length*n of forms ####

## Plotting avg length against slope:

plot_data %>%
  group_by(group) %>%
  summarize(mean_slope = mean(fixef(model_shuffle)[2] + slope), 
            mean_length = mean(Length)) %>%
  ggplot(aes(x=mean_length, y=mean_slope)) +
  geom_smooth(method = 'lm') +
  geom_point() + # annotate with geom_repel
  xlab('Average length') +
  ylab('Slope')

data_corr <- plot_data %>%
  group_by(group) %>%
  summarize(mean_slope = mean(fixef(model_shuffle)[2] + slope),
            mean_length = mean(Length))

cor.test(data_corr$mean_slope, data_corr$mean_length, method='pearson')


#4. Plotting the map

## Get intercept values
random_effects <- ranef(model_shuffle, condVar = TRUE)
random_effects_l <- random_effects$`Glottocode:Family`
random_effects_l$z_length <- fixef(model_shuffle)[2] +random_effects_l$z_length
random_effects_l <- cbind(language = rownames(random_effects_l), 
                          random_effects_l)
row.names(random_effects_l) <- NULL
random_effects_l <- random_effects_l %>%
  separate(language, into = c("language", "family"), sep = ":")
random_effects_l <- random_effects_l %>% rename(Glottocode = language)
## Get lattitudes and longitudes
languages_info <- read.csv('data/lexibank-lexibank-analysed-a4c0952/cldf/languages.csv') %>%
  distinct(Glottocode, .keep_all = TRUE) %>%
  select(Glottocode, Latitude, Longitude, Name)
## Merge the datasets
languages_info <- languages_info %>%
  inner_join(random_effects_l, by = "Glottocode")
## Number of forms
forms_n <- logistic_data %>%
  distinct(Glottocode, .keep_all = T) %>%
  select(Glottocode, n_forms)
languages_info <- languages_info %>%
  inner_join(forms_n, by = "Glottocode")

## Map data
world <- map_data("world") %>% 
  filter(lat > -62)

## full map 
# Identify the top-10 lowest z_length values
top10_lowest_z_length <- languages_info %>%
  arrange(z_length) %>%
  head(10)
## Create the plot with ggrepel
full_map_results <- ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group), 
               colour = NA, fill = 'grey', alpha = 0.4) +
  geom_point(data = languages_info, aes(x = Longitude, y = Latitude, 
                                        color = z_length, size = n_forms)) +
  scale_color_gradientn(
    colors = c('#c55a01', 'white', '#0097c3'), 
    values = scales::rescale(c(min_z_length, 0, max_z_length)),
    guide = guide_colorbar(title = expression(beta * "-coefficient for length:"))
  ) +
  scale_size_continuous(
    name = "Number of forms:"
  ) +
  geom_text_repel(data = top10_lowest_z_length, aes(x = Longitude, 
                                                    y = Latitude, label = Name),
                  size = 4, # Adjust the size as needed
                  nudge_y = 10, # Adjust the nudging as needed
                  segment.size = 0.2) + # Size of the segment line connecting text to points
  theme_void() +
  theme(
    legend.position = 'bottom',
    legend.direction = "horizontal",
    legend.title = element_text(vjust = 0.8)  # Center align the legend title
  )
ggsave('figures/map_results.png', full_map_results, width = 10, height = 5)
## full basic map 
palette_colors <- wes_palette("Zissou1", 100, type = "continuous")
full_basic_map <- ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group), 
               colour = NA, fill = 'grey', alpha = 0.4) +
  geom_point(data = languages_info, aes(x = Longitude, y = Latitude, color = n_forms)) +
  scale_color_gradientn(
    colours = palette_colors,  # Use the chosen wesanderson palette
    trans = "log10",  # Apply log10 transformation
    name = "Number of forms"  # Title for the colorbar
  ) +
  theme_void() +
  theme(
    legend.position = 'bottom',
    legend.direction = "horizontal",
    legend.title = element_text(vjust = 0.8),  # Center align the legend title
    legend.key.width = unit(2, 'cm')  # Increase the width of the colorbar
  )
ggsave('figures/basic_map.png', full_basic_map, width = 10, height = 5)
