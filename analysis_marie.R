#j'ai appelé les 3 modèles model_null, model_no_rd_slope, model
#et voilà les analyses que j'ai faites même si tu m'as dit qu'avec geom_smooth ça n'allait pas

install.packages("DHARMa")
library(DHARMa)
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