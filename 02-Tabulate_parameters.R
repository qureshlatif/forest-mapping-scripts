library(stringr)
library(dplyr)
library(R.utils)

setwd("C:/Users/Quresh.Latif/files/projects/CEAP_map_tool")
load("Data_compiled.RData")

#__________ Script inputs _____________#
mod.nam <- "mod"
scripts.loc <- "forest-mapping-scripts/"
mod <- loadObject(mod.nam)
out.vals <- c("est", "f")
#______________________________________#

mod$summary %>% write.csv(str_c("Params_all.csv"), row.names = F) # Save full table of parameter estimates

# Compile data to get covariate names #
source(str_c(scripts.loc, "Data_processing_community_occupancy.R"))

params <- c("psi0", "psi.2014", "psi.2016", "psi.2018", str_c("beta.", dimnames(X.beta)[[3]]),
            "theta0", str_c("alpha.", dimnames(X.alpha)[[2]]),
            "pStar", str_c("zeta.", dimnames(X.zeta)[[2]]))
cols <- (expand.grid(out.vals, params, stringsAsFactors = F) %>%
  select(Var2, Var1) %>%
  mutate(Var3 = str_c(Var2, Var1, sep = ".")))$Var3
cols <- cols[-which(str_detect(cols, ".f") & (str_detect(cols, "psi") | str_detect(cols, "theta0") | str_detect(cols, "pStar")))]
out <- matrix(NA, nrow = length(spp.list), ncol = length(cols),
              dimnames = list(spp.list, cols))

parm <- QSLpersonal::expit(mod$sims.list$beta0)
out[, "psi0.est"] <- str_c(
  apply(parm, 2, median) %>% round(digits = 2),
  "(",
  apply(parm, 2, function(x) quantile(x, prob = 0.025, type = 8)) %>% round(digits = 2),
  ",",
  apply(parm, 2, function(x) quantile(x, prob = 0.975, type = 8)) %>% round(digits = 2),
  ")")

parm <- QSLpersonal::expit(mod$sims.list$alpha0)
out[, "theta0.est"] <- str_c(
  apply(parm, 2, median) %>% round(digits = 2),
  "(",
  apply(parm, 2, function(x) quantile(x, prob = 0.025, type = 8)) %>% round(digits = 2),
  ",",
  apply(parm, 2, function(x) quantile(x, prob = 0.975, type = 8)) %>% round(digits = 2),
  ")")

parm <- QSLpersonal::expit(mod$sims.list$zeta0) %>% (function(x) 1 - (1-x)^6)
out[, "pStar.est"] <- str_c(
  apply(parm, 2, median) %>% round(digits = 2),
  "(",
  apply(parm, 2, function(x) quantile(x, prob = 0.025, type = 8)) %>% round(digits = 2),
  ",",
  apply(parm, 2, function(x) quantile(x, prob = 0.975, type = 8)) %>% round(digits = 2),
  ")")

pars <- c("psi.2014", "psi.2016", "psi.2018")
for(i in 1:length(pars)) {
  parm <- QSLpersonal::expit(mod$sims.list$beta0 + mod$sims.list$dev.b0[,,i])
  out[, str_c(pars[i], ".est")] <- str_c(
    apply(parm, 2, median) %>% round(digits = 2),
    "(",
    apply(parm, 2, function(x) quantile(x, prob = 0.025, type = 8)) %>% round(digits = 2),
    ",",
    apply(parm, 2, function(x) quantile(x, prob = 0.975, type = 8)) %>% round(digits = 2),
    ")")
}

pars <- params[which(str_detect(params, "beta."))]
for(i in 1:length(pars)) {
  parm <- mod$sims.list$beta1[,,i]
  out[, str_c(pars[i], ".est")] <- str_c(
    apply(parm, 2, median) %>% round(digits = 2),
    "(",
    apply(parm, 2, function(x) quantile(x, prob = 0.025, type = 8)) %>% round(digits = 2),
    ",",
    apply(parm, 2, function(x) quantile(x, prob = 0.975, type = 8)) %>% round(digits = 2),
    ")")
  out[, str_c(pars[i], ".f")] <-
    apply(parm, 2, function(x) max(c(sum(x > 0), sum(x < 0))) / length(x)) %>%
    round(digits = 2)
}

pars <- params[which(str_detect(params, "alpha."))]
for(i in 1:length(pars)) {
  parm <- mod$sims.list$alpha1[,,i]
  out[, str_c(pars[i], ".est")] <- str_c(
      apply(parm, 2, median) %>% round(digits = 2),
      "(",
      apply(parm, 2, function(x) quantile(x, prob = 0.025, type = 8)) %>% round(digits = 2),
      ",",
      apply(parm, 2, function(x) quantile(x, prob = 0.975, type = 8)) %>% round(digits = 2),
      ")")
  out[, str_c(pars[i], ".f")] <-
      apply(parm, 2, function(x) max(c(sum(x > 0), sum(x < 0))) / length(x)) %>%
      round(digits = 2)
}

pars <- params[which(str_detect(params, "zeta."))]
for(i in 1:length(pars)) {
  parm <- mod$sims.list$zeta1[,,i]
  out[, str_c(pars[i], ".est")] <- str_c(
    apply(parm, 2, median) %>% round(digits = 2),
    "(",
    apply(parm, 2, function(x) quantile(x, prob = 0.025, type = 8)) %>% round(digits = 2),
    ",",
    apply(parm, 2, function(x) quantile(x, prob = 0.975, type = 8)) %>% round(digits = 2),
    ")")
  out[, str_c(pars[i], ".f")] <-
    apply(parm, 2, function(x) max(c(sum(x > 0), sum(x < 0))) / length(x)) %>%
    round(digits = 2)
}

write.csv(out, "Parameter_est.csv", row.names = T)
