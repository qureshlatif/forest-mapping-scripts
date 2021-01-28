# Primary data objects
beta.vars <- c("PACC10", "PACC40", "mnPerArRatio_Opn", "Latitude", "heatload", "TWI")
alpha.vars <- c("CanCov")
alpha.quad <- c(T)
zeta.vars <- c("CanCov", "DayOfYear", "Time_ssr")
zeta.quad <- c(F, T, T)

p.beta <- length(beta.vars) # Number of covariates on grid-level occupancy
p.alpha <- length(alpha.vars) + sum(alpha.quad) # Number of covariates on point-level occupancy
p.zeta <- length(zeta.vars) + sum(zeta.quad) # Number of covariates on detection

# Detection data #
Y <- Y.mat
TPeriod <- TR.mat
gridID <- Cov[, "gridIndex"]
yearID <- Cov[, "YearInd"]
yearID.grid <- landscape_data$YearInd
gridID.grid <- landscape_data$gridIndex
n.grid <- max(gridID)
n.year <- max(yearID)
n.gridXyear <- nrow(landscape_data)
n.point_year <- dim(Y)[1]
n.spp <- dim(Y)[2]

# Grid-level covariates #
X.beta.long <- landscape_data %>%
  select(one_of(beta.vars)) %>%
  data.matrix() %>%
  apply(2, function(x) (x - mean(x)) / sd(x))

X.beta <- array(0, dim = c(n.grid, n.year, p.beta))
for(i in 1:nrow(X.beta.long)) {
  X.beta[gridID.grid[i], yearID.grid[i], ] <- X.beta.long[i,]
}

dimnames(X.beta.long)[[2]] <- dimnames(X.beta)[[3]] <- beta.vars

# Point-level covariates #
X.alpha <- Cov[, alpha.vars] %>%
  as.matrix() %>%
  apply(2, (function(x) (x - mean(x)) / sd(x)))
for(v in 1:length(alpha.vars)) {
  if(alpha.quad[v]) {
    alpha.vars <- c(alpha.vars, str_c(alpha.vars[v], "2"))
    X.alpha <- cbind(X.alpha, X.alpha[,v]^2)
  }
}
dimnames(X.alpha)[[2]] <- alpha.vars

X.zeta <- Cov[, zeta.vars] %>%
  as.matrix() %>%
  apply(2, (function(x) (x - mean(x)) / sd(x)))
for(v in 1:length(zeta.vars)) {
  if(zeta.quad[v]) {
    zeta.vars <- c(zeta.vars, str_c(zeta.vars[v], "2"))
    X.zeta <- cbind(X.zeta, X.zeta[,v]^2)
  }
}
dimnames(X.zeta)[[2]] <- zeta.vars

rm(v)
