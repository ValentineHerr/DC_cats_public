# ---- analysing density results ---- ###


# clear environement ####

rm(list = ls())# ---- putting catSMR data together ---- ###


# clear environement ####

rm(list = ls())

# load libraries ####
library(raster)
library(coda)
library(progress)
library(rgeos)
library(dplyr)
library(lme4)
library(ggplot2)

## prepare function to write equation of model
model_equation <- function(model, ...) {
  format_args <- list(...)

  model_coeff <- model$coefficients
  format_args$x <- abs(model$coefficients)
  model_coeff_sign <- sign(model_coeff)
  model_coeff_prefix <- case_when(model_coeff_sign == -1 ~ " - ",
                                  model_coeff_sign == 1 ~ " + ",
                                  model_coeff_sign == 0 ~ " + ")
  model_eqn <- paste(strsplit(as.character(model$call$formula), "~")[[2]], # 'y'
                     "=",
                     paste(if_else(model_coeff[1]<0, "- ", ""),
                           do.call(format, format_args)[1],
                           paste(model_coeff_prefix[-1],
                                 do.call(format, format_args)[-1],
                                 " * ",
                                 names(model_coeff[-1]),
                                 sep = "", collapse = ""),
                           sep = ""))
  return(model_eqn)
}

# load data ####
load("C_PostProcessing_of_Densities/density_results.Rdata")


# load("Occupancy_results.Rdata")

# run model ####

fm <- glm(round(Y*10)~pop_density_km2_scaled + I(pop_density_km2_scaled^2) + Median_Federal_Adjusted_Gross_Income_2015_scaled, data = data, family = "poisson", na.action = "na.fail")



summary(fm)

all_models <- MuMIn::dredge(fm, subset = !('I(pop_density_km2_scaled^2)' && !'pop_density_km2_scaled')) # try all models but keepinp first order term if second is in there



View(all_models)

best_model <- eval(attr(all_models, "model.calls")[[1]])
summary(best_model)

newdata <- rbind(
  expand.grid(
    Median_Federal_Adjusted_Gross_Income_2015_scaled = seq(min(data$Median_Federal_Adjusted_Gross_Income_2015_scaled), max(data$Median_Federal_Adjusted_Gross_Income_2015_scaled), length.out = 100),
    distance_to_residential_scaled = round(quantile(data$distance_to_residential_scaled, c(0, 1)), 2),
    pop_density_km2_scaled = median(data$pop_density_km2_scaled),
    variable = "Income",
    Mall = "No"),
  expand.grid(
    Median_Federal_Adjusted_Gross_Income_2015_scaled = median(data$Median_Federal_Adjusted_Gross_Income_2015_scaled),
    distance_to_residential_scaled = round(quantile(data$distance_to_residential_scaled, c(0, 1)), 2),
    pop_density_km2_scaled = seq(min(data$pop_density_km2_scaled), max(data$pop_density_km2_scaled), length.out = 100),
    variable = "Pop_density",
    Mall = "No")
  )


newdata$x <- c(newdata$Median_Federal_Adjusted_Gross_Income_2015_scaled[newdata$variable %in% "Income"] * sd_demographic["Median_Federal_Adjusted_Gross_Income_2015"] + mean_demographic["Median_Federal_Adjusted_Gross_Income_2015"] , newdata$pop_density_km2_scaled[newdata$variable %in% "Pop_density"] * sd_demographic["pop_density_km2"] + mean_demographic["pop_density_km2"])

newdata$Dist_res <- as.factor(ceiling(newdata$distance_to_residential_scaled * sd_distance_to_residential + mean_distance_to_residential ))

newdata <- cbind(newdata, predict(best_model, newdata = newdata,type="link",re.form=NA, se.fit = T))
newdata$exp.fit <- exp(newdata$fit)/10
newdata$LCI <- exp(newdata$fit - 1.96*newdata$se.fit)/10
newdata$UCI <- exp(newdata$fit + 1.96*newdata$se.fit)/10






# lines(exp.fit ~ Median_Federal_Adjusted_Gross_Income_2015_scaled, data = newdata)
p_income <- ggplot(newdata[newdata$variable %in% "Income",], aes(y = exp.fit, x = x)) + geom_line()+
  geom_ribbon(aes(ymin = LCI, ymax = UCI, group =Dist_res), alpha = 0.5) +
  xlab("Median Income ($)") +
  geom_point(data = data.frame(exp.fit =data$Y, #3),
                               x = c(data$Median_Federal_Adjusted_Gross_Income_2015_scaled* sd_demographic["Median_Federal_Adjusted_Gross_Income_2015"] + mean_demographic["Median_Federal_Adjusted_Gross_Income_2015"])), shape=16, alpha = 0.5) +
  ylab(expression(Density~(cat~ha^-1)))+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line =   element_line())

p_pop <- ggplot(newdata[newdata$variable %in% "Pop_density",], aes(y = exp.fit, x = x
                                                              # , group =Dist_res
)) + geom_line() +
  geom_ribbon(aes(ymin = LCI, ymax = UCI, group =Dist_res), alpha = 0.5) +
  xlab("Population density (per km2)") +
  geom_point(data = data.frame(exp.fit =data$Y, #3),
                               x = c(data$pop_density_km2_scaled * sd_demographic["pop_density_km2"] + mean_demographic["pop_density_km2"] )), shape=16, alpha = 0.5) +
  ylab("") +
  theme(panel.grid = element_blank(),
              panel.background = element_blank(),
              axis.line =   element_line())


ggpubr::ggarrange(p_income,p_pop,

                  ncol = 2,

                  nrow = 1,

                  common.legend = T,

                  labels = c("a)", "b)"),

                  font.label = list(size = 11, color = "black", face = "plain"),

                  hjust = -2,

                  vjust = 1)



ggsave("C_PostProcessing_of_Densities/Density_response.png", width = 8, height = 5)


1-(best_model$deviance/best_model$null.deviance) # 0.43 was 0.47 for occupancy >0.5

(abund_AICtab <- all_models)
(abund_best_model <- best_model)
