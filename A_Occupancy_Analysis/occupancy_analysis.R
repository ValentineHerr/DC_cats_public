# --- Run Occupancy models (not Cat ID) --- #

# clear environment ####
rm(list = ls())

# load libraries ####
library(unmarked)

# load data ####
load("Occupancy_Analysis/occupancy_data.Rdata")

# make pretty variable names ####
v_names <- list(
  "p(Int)" = "$p$(int)",
  "psi(Int)" = "$\\psi$(Int)",
  "p(bait)" = "$p$($bait$)",
  "p(placement_type)" = "$p$($PT$)",
  "p(distance_to_obstruction_scaled)"="$p$($dist_{obstr}$)",
  "p(I(distance_to_obstruction_scaled^2))"="$p$($dist_{obstr}^2$)",
  "psi(distance_to_residential_scaled)" = "$\\psi$($dist_{res}$)",
  "psi(pop_density_km2_scaled)" = "$\\psi$($pop_{dens}$)",
  "psi(I(pop_density_km2_scaled^2))" = "$\\psi$($pop_{dens}^2$)",
  "psi(Mall:Median_Federal_Adjusted_Gross_Income_2015_scaled)" = "$\\psi$($mall$:$income$)",
  "psi(distance_to_residential_scaled:Mall:Median_Federal_Adjusted_Gross_Income_2015_scaled)" = "$\\psi$($Mall$:$income$:$dist_{res}$)",
  "df" = "df",
  "logLik" = "logLik",
  "AICc" = "AICc",
  "delta" = "delta",
  "weight" = "weight",
  "(Intercept)" = "(Intercept)",
  "distance_to_residential_scaled" = bquote("dist"[res]),
  "pop_density_km2_scaled" = bquote("pop"[dens]),
  "I(pop_density_km2_scaled^2)" = bquote("pop"[dens]^2),
  "MallNo:Median_Federal_Adjusted_Gross_Income_2015_scaled" = bquote(italic("income")),
  "MallYes:Median_Federal_Adjusted_Gross_Income_2015_scaled" = bquote("Mall"~italic("income")~" (irrelevant)"),
  "distance_to_residential_scaled:MallYes:Median_Federal_Adjusted_Gross_Income_2015_scaled" = bquote("dist"[res]*":mall"~italic("income")~" (irrelevant)"),
  "distance_to_residential_scaled:MallNo:Median_Federal_Adjusted_Gross_Income_2015_scaled" = bquote("dist"[res]*":"*italic("income")),
  # "baitno" = "bait no",
  "bait" = bquote(italic("bait")),
  "baityes" = bquote(italic("bait")~"yes"),
  "distance_to_obstruction_scaled" =bquote("dist"[obst]),
  "I(distance_to_obstruction_scaled^2)" = bquote("dist"[obst]^2),
  "placement_typeprivate" = bquote(italic("PT")~"private"),
  # "placement_typepublic forest" = bquote(italic("PT")~"public forest"),
  # "placement_typepublic open" = bquote(italic("PT")~"public open"),
  "placement_typepublic" = bquote(italic("PT")~"public"),
  "placement_type" = bquote(italic("PT"))
)

## run occupancy models

fm <- occu(~ bait + distance_to_obstruction_scaled + I(distance_to_obstruction_scaled^2) + placement_type  ~ pop_density_km2_scaled + I(pop_density_km2_scaled^2) + distance_to_residential_scaled * Median_Federal_Adjusted_Gross_Income_2015_scaled:Mall , data = umf)  # ○fit a model

all_models <- MuMIn::dredge(fm, subset = !(('p(I(distance_to_obstruction_scaled^2))'  && !'p(distance_to_obstruction_scaled)') | ('psi(I(pop_density_km2_scaled^2))' && !'psi(pop_density_km2_scaled)'))) # try all models but keepinp first order term if second is in there

View(all_models)



# look at best model ####

best_model <- eval(attr(all_models, "model.calls")[[1]])
# best_model <- update(best_model, ~.-1~.-1) # removing intercept
summary(best_model)

## cross validation ####
(cvlist <- crossVal(best_model, method='Kfold'))


## coeficient plot ####
C_occ <- best_model['state']@estimates
C_occ <- data.frame(C_occ, confint(best_model, type = "state"))
names(C_occ) <- c("estimate", "LCI", "UCI")

C_occ <- C_occ[rev(rownames(C_occ)),]

C_occ$parameter <-  factor(rownames(C_occ), levels = rownames(C_occ))


C_p <- best_model['det']@estimates
C_p <- data.frame(C_p, confint(best_model, type = "det"))
names(C_p) <- c("estimate", "LCI", "UCI")
C_p<- C_p[rev(rownames(C_p)),]

C_p$parameter <- factor(rownames(C_p), levels = rownames(C_p))

library(ggplot2)
ggplot(data = C_occ) +
  geom_segment(aes(x = LCI, xend = UCI, y = parameter, yend =parameter)) +
  geom_point(aes(estimate, parameter)) +
  geom_vline(xintercept = 0) +
  theme_classic()

ggplot(data = C_p) +
  geom_segment(aes(x = LCI, xend = UCI, y = parameter, yend =parameter)) +
  geom_point(aes(estimate, parameter)) +
  geom_vline(xintercept = 0) +
  theme_classic()

ggplot(data = C_p) +
  geom_segment(aes(x = LCI, xend = UCI, y = parameter, yend =parameter)) +
  geom_point(aes(estimate, parameter)) +
  geom_vline(xintercept = 0) +
  labs(x = "Estimated value") +
  scale_y_discrete(labels=v_names[as.character(C_p$parameter)]) +
  theme_classic()

