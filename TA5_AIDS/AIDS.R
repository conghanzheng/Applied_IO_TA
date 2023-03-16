## ----r_setup, echo = knitr::is_html_output()----------------------------------

## Preliminaries ----
rm(list = ls()) ## Clear the environment

## Install and load packages
packages <- c('rstudioapi','R.utils','dplyr','tidyr','data.table','knitr','kableExtra','Hmisc','systemfit','miscTools') ## Manually list package dependencies
if (!require('pacman')) install.packages('pacman')
library(pacman)
pacman::p_load(packages, character.only = TRUE) ## Check if all packages declared in the previous vector are installed, install if not, and load them

## (Local) Path to the current script
scriptpath <- dirname(rstudioapi::getSourceEditorContext()$path) %>% setwd()

## ----data_clean, echo = knitr::is_html_output(), results = 'markup'-----------

## Online data: Blanciforti et al. (1986) ----
data_wide <- data.table::fread("https://github.com/cran/micEconAids/raw/master/data/Blanciforti86.txt.gz") %>% 
  drop_na() %>%
  rename(Year = V1)

## Number of goods
n <- 4

## Data cleansing ----
## Eventually, we need wide data, but it would be handy if we first reshape it longer
data_long <- tidyr::pivot_longer(data_wide %>% select(c('Year',"xFood",
                                                        'wFood1','wFood2','wFood3','wFood4',
                                                        'pFood1','pFood2','pFood3','pFood4')),
                         cols = !c(Year,xFood), 
                         names_to = c('.value','Food'), 
                         names_pattern = "(\\w+)(\\d{1})") %>%
  group_by(Food) %>%
  rename(w_jt = wFood,
         p_jt = pFood,
         x_t = xFood) %>%
  mutate(w_j = mean(w_jt, na.rm = TRUE),
         p_j = mean(p_jt, na.rm = TRUE),
         lnp_jt = log(p_jt)) %>%
  group_by(Year) %>%
  mutate(lnx_t = log(x_t),
         lnP_stone_t = sum(w_jt*lnp_jt, na.rm = TRUE),
         lnP_stone_t_1 = sum(Hmisc::Lag(w_jt)*lnp_jt, na.rm = TRUE),
         lnx_lnP_stone_t = log(x_t) - lnP_stone_t,
         lnx_lnP_stone_t_1 = log(x_t) - lnP_stone_t_1)

print(head(data_long[,1:7]))

## ----la_functions, echo = knitr::is_html_output()-----------------------------

## LA-AIDS Data ----
la_data <- pivot_wider(data_long %>% select(c('Year','Food','x_t','w_jt','lnp_jt','lnx_lnP_stone_t','lnx_lnP_stone_t_1')),
                       names_from = Food,
                       names_glue = "{.value}_{Food}",
                       values_from = c(w_jt,lnp_jt)
                       )

colnames(la_data) <- gsub("(^\\w+\\_)(j)(t)(\\_)(\\d{1})","\\1\\5\\3",colnames(la_data))

##' Function: return the restrictions on parameters ---- 
##' (our equation system is a constrained system)
##' Symmetry is assumed (notice that symmetry implies homogeneity)
restriction <- function(n) {
  nExoEqn <- n + 2 # number of exogenous variables per equation
  nExo <- (n-1)*nExoEqn # number of exogenous variables 
  
  restriction <- diag(nExo)
  delCols <- NULL
  for (i in 1:(n-1)) { ## homogeneity
    delCol <- (i-1)*(nExoEqn) + 2 + n
    addCol <- ((i-1)*nExoEqn + 3):((i - 1)*nExoEqn + 2 + (n - 1))
    restriction[delCol, addCol] <- -1
    delCols <- c(delCols, delCol)
  }
  for (i in 1:(n - 2)) { ## symmetry
    for (j in (i + 1):(n - 1)) {
     delCol <- (j - 1) * nExoEqn + 2 + i
     addCol <- (i - 1) * nExoEqn + 2 + j
     restriction[,addCol] <- restriction[,addCol] + restriction[,delCol]
     delCols <- c(delCols, delCol)
     }
  }
  restriction <- restriction[,-delCols]
  
  tmpnames <- NULL
   for( i in 1:(n - 1)) {
      tmpnames <- c(tmpnames, paste0("alpha_",i),paste0("beta_",i))
      start <- 1
      stop  <- n
      for( j in 1:n){
         tmpnames <- c(tmpnames, paste0("gamma_",i,j))
      }
   }
  
  rownames(restriction) <- tmpnames
  
  tmpnames <- NULL
   for( i in 1:(n-1)) {
      tmpnames <- c(tmpnames, paste0("alpha_",i),paste0("beta_",i))
      start <- 1
      stop  <- n
      for( j in i:(n-1)){
         tmpnames <- c(tmpnames, paste0("gamma_",i,j))
      }
   }
  
  colnames(restriction) <- tmpnames
  
  return(restriction)
}

##' Function: LA-AIDS estimation ----
##' returning the coefficients
la_coef <- function(est, n, data, stone_lag = FALSE) {
  data <- as.data.frame(data)
  
  nExoEqn <- n + 2 # number of exogenous variables per equation
  nObs <- if(stone_lag == TRUE) c(2:nrow(data)) else c(1:nrow(data))
  
  ## Matrix used to transform the coefficients
  M <- matrix(0, n*nExoEqn, (n - 1)*nExoEqn)
  rownames(M) <- c(paste0("alpha_",c(1:n)),paste0("beta_",c(1:n)),paste0("gamma_",rep(1:n,each=n),rep(1:n,n)))
  tmpnames <- NULL
   for( i in 1:(n - 1)) {
      tmpnames <- c(tmpnames, paste0("alpha_",i),paste0("beta_",i))
      start <- 1
      stop  <- n
      for( j in 1:n){
         tmpnames <- c(tmpnames, paste0("gamma_",i,j))
      }
   }
  colnames(M) <- tmpnames
  alphas <- paste0("alpha_",c(1:n))
  betas <- paste0("beta_",c(1:n))
  gammas <- matrix(paste0("gamma_",rep(1:n,n),rep(1:n,each = n)), nrow = n, ncol = n)
  
  for(i in 1:(n-1) ) {
         M[alphas[i], alphas[i]] <- 1
         M[alphas[n], alphas[i]] <- -1 
         M[betas[i], betas[i]] <- 1 
         M[betas[n], betas[i]] <- -1
         for(j in 1:n) {
            M[gammas[i, j], gammas[i, j]] <- 1 
            M[gammas[n, j], gammas[i, j]] <- -1 
         }
  }
  
  ## Coefficients
  coef <- c(M %*% coef(est))
  names(coef) <- c(paste0("alpha_",c(1:n)),paste0("beta_",c(1:n)),paste0("gamma_",rep(1:n,each=n),rep(1:n,n)))
  coef[alphas[n]]  <- coef[alphas[n]] + 1
  
  ## Variance-covariance matrix
  cov <- M %*% vcov(est) %*% t(M)
  rownames(cov) <- names(coef)
  colnames(cov) <- names(coef)
  
  ## Estimates with se and p values
  stat <- miscTools::coefTable(coef, sqrt(diag(cov)), 1)
  
  ## Fitted values
  wFitted <- as.data.frame(matrix( NA, nrow = nrow(data), ncol = n))
  colnames(wFitted) <- paste0("w_",c(1:n),'t')
  rownames(wFitted) <- rownames(data)
  
  if (stone_lag == FALSE) { ## whether to use lagged shares in Stone's price index
    for( i in 1:nrow(data)) {
         logPrices <- as.numeric(data[i,paste0("lnp_",c(1:n),'t')])
         logTotExp <- as.numeric(log(data[i,'x_t']))
         if(all(!is.na(c(logPrices,logTotExp)))) {
            numerator <- coef[grepl('alpha',names(coef))] + matrix(coef[grepl('gamma',names(coef))], nrow = 4) %*% logPrices + coef[grepl('beta',names(coef))] * logTotExp
            wFitted[i,] <- solve(diag(n) + coef[grepl('beta',names(coef))] %*% t(logPrices), numerator)
         }
      }
  } else {
    for( i in 1:nrow(data)) {
         logPrices <- as.numeric(data[1,paste0("lnp_",c(1:n),'t')])
         logTotExp <- as.numeric(log(data[1,'x_t']))
         if(all(!is.na(c(logPrices,logTotExp)))) {
            numerator <- coef[grepl('alpha',names(coef))] + matrix(coef[grepl('gamma',names(coef))], nrow = 4) %*% logPrices + coef[grepl('beta',names(coef))] * logTotExp
            wFitted[1,] <- solve(diag(n) + coef[grepl('beta',names(coef))] %*% t(logPrices), numerator)
         }
      }
  }
  
  ## Residuals
  wResid <- data.frame(matrix(NA, nrow = nrow(data), ncol = n))
  names(wResid) <- paste0("wResid", as.character(1:n))
  
  for( i in 1:n) {
    wResid[,i] <- data[[paste0("w_",c(1:n),'t')[i]]] - wFitted[,i]
   }
  
  ## R-squared
  r2 <- numeric(n)
  for(i in 1:(n - 1)) {
    r2[i] <- summary(est$eq[[i]])$r.squared
  }
  r2[n] <- miscTools::rSquared(as.numeric(data[nObs, paste0("w_",c(1:n),'t')[n]]), wResid[nObs,n])
  names(r2) <- paste0("w_",c(1:n),'t')
  
  ## Function output
  out <- list()
  out$coef <- coef
  out$alpha <- coef[grepl('alpha',names(coef))]
  out$beta <- coef[grepl('beta',names(coef))]
  out$gamma <- matrix(coef[grepl('gamma',names(coef))], nrow = n)
  out$cov <- cov
  out$stat <- stat
  out$r2 <- r2
  
  return(out)
}

## Function: print the LA-AIDS estimation results ----
la_print <- function(coef) {
  cat( "Demand analysis with the Almost Ideal Demand System (AIDS)\n" )
  cat( "Estimation Method: LA" )
  cat( "Estimated Coefficients:\n" )
  print(coef$stat)
  cat( "\n R-squared:\n" )
  print(coef$r2)
}

## ----la_system, echo = knitr::is_html_output(), restuls = 'markup'------------
## Function: LA equation system ----
la_sur_system <- function(n,stone_lag=FALSE){
  system <- list()
    for(i in 2:n ) {
      if (stone_lag==TRUE) {
        system[[i-1]] <- paste( "w_", as.character(i), 't', " ~ lnx_lnP_stone_t_1", sep = "" )
      } else {
        system[[i-1]] <- paste( "w_", as.character(i), 't', " ~ lnx_lnP_stone_t", sep = "" )
      }
      for( j in 1:n ) {
        system[[i-1]] <- paste(system[[i-1]], " + lnp_", as.character(j), 't', sep = "")
      }
      system[[i-1]] <- as.formula(system[[i-1]])
    }
  return(system)
}

## Print
print(la_sur_system(n=4, stone_lag = FALSE) %>% as.matrix())

## ----la_fit, echo = knitr::is_html_output(), results = 'markup'---------------
## Execute the functions
la_est <- systemfit::systemfit(la_sur_system(n=4, stone_lag = FALSE), 
                               method = "SUR",
                               data = la_data, 
                               restrict.regMat = restriction(n=4))

## Print results
la_print(la_coef(la_est, n=4, data = la_data))

## ----aids_system, echo = knitr::is_html_output(), restuls = 'markup'----------
## Function: AIDS equation system ----
aids_sur_system <- function(n){
  system <- list()
    for(i in 2:n ) {
      system[[i-1]] <- paste( "w_", as.character(i), 't', " ~ lnx_lnP_t", sep = "" )
      for( j in 1:n ) {
        system[[i-1]] <- paste(system[[i-1]], " + lnp_", as.character(j), 't', sep = "")
      }
      system[[i-1]] <- as.formula(system[[i-1]])
    }
  return(system)
}

## Print the system
print(aids_sur_system(n=4) %>% as.matrix())

## ----aids_funtions, echo = knitr::is_html_output(), results = 'markup'--------
## AIDS data ----
aids_data <- pivot_wider(data_long %>% 
                           select(c('Year','Food','x_t','lnx_t','w_jt','lnp_jt')),
                       names_from = Food,
                       names_glue = "{.value}_{Food}",
                       values_from = c(w_jt,lnp_jt)
                       )

colnames(aids_data) <- gsub("(^\\w+\\_)(j)(t)(\\_)(\\d{1})","\\1\\5\\3",colnames(aids_data))

##' Function: to calculate the translog prices lnP ---- 
##' (which depend on parameters)
translog <- function(coef, data) {
  nObs <- nrow(data)
  alpha0 <- 1
  lnP_t <- array(alpha0, c(nObs))
  
  alpha <- coef[grepl('alpha',names(coef))]
  n <- length(alpha)
  beta <- coef[grepl('beta',names(coef))]
  gamma <- matrix(coef[grepl('gamma',names(coef))], nrow = n)

  for( i in 1:n) {
    lnP_t <- lnP_t + alpha[i] * data[[paste0("lnp_",i,'t')]]
    for( j in 1:n) {
      lnP_t <- lnP_t + 0.5 * gamma[i,j] *data[[paste0("lnp_",i,'t')]] * data[[paste0("lnp_",j,'t')]]
      }
  }
  
  names(lnP_t) <- paste0('t=',row.names(data))
  
  return(lnP_t)
}

## Add intial values of translog prices to data (calculated using LA-AIDS estiamtes)
aids_data$lnP_t <- as.numeric(translog(coef = la_coef(la_est, n=4, data = la_data)$coef, data = la_data))

aids_data <- aids_data %>%
  mutate(lnx_lnP_t = lnx_t - lnP_t)

## Function: AIDS estimation ----
aids_est <- function(maxIter, tol = 1e-5, b_start = la_coef(la_est, n=4, data = la_data)$coef, data = aids_data, n) {

## starting values
b <- b_start
b_diff <- b
data_updated <- data
count <- 1

for (nIter in 1:maxIter) {
  b_last <- b ## save the coef from last iteration
  
  lnP_update <- translog(b, data = aids_data)
  data_updated$lnP_t <- as.numeric(lnP_update)
  data_updated <- data_updated %>% mutate(lnx_lnP_t = lnx_t - lnP_t)
  
  ## SUR
  est <- systemfit::systemfit(aids_sur_system(n), method = "SUR",data = data_updated,restrict.regMat = restriction(n))
  
  nExoEqn <- n + 2 # number of exogenous variables per equation
  nObs <- c(1:nrow(data_updated))
  
  ## Matrix used to transform the coefficients
  M <- matrix(0, n*nExoEqn, (n - 1)*nExoEqn)
  rownames(M) <- c(paste0("alpha_",c(1:n)),paste0("beta_",c(1:n)),paste0("gamma_",rep(1:n,each=n),rep(1:n,n)))
  tmpnames <- NULL
   for( i in 1:(n - 1)) {
      tmpnames <- c(tmpnames, paste0("alpha_",i),paste0("beta_",i))
      start <- 1
      stop  <- n
      for( j in 1:n){
         tmpnames <- c(tmpnames, paste0("gamma_",i,j))
      }
   }
  colnames(M) <- tmpnames
  alphas <- paste0("alpha_",c(1:n))
  betas <- paste0("beta_",c(1:n))
  gammas <- matrix(paste0("gamma_",rep(1:n,n),rep(1:n,each = n)), nrow = n, ncol = n)
  
  for(i in 1:(n-1) ) {
         M[alphas[i], alphas[i]] <- 1
         M[alphas[n], alphas[i]] <- -1 
         M[betas[i], betas[i]] <- 1 
         M[betas[n], betas[i]] <- -1
         for(j in 1:n) {
            M[gammas[i, j], gammas[i, j]] <- 1 
            M[gammas[n, j], gammas[i, j]] <- -1 
         }
  }
  
  ## Coefficients
  b <- c(M %*% coef(est))
  names(b) <- c(paste0("alpha_",c(1:n)),paste0("beta_",c(1:n)),paste0("gamma_",rep(1:n,each=n),rep(1:n,n)))
  b[alphas[n]]  <- b[alphas[n]] + 1
  
  alpha <- b[grepl('alpha',names(b))]
  beta <- b[grepl('beta',names(b))]
  gamma <- matrix(b[grepl('gamma',names(b))], nrow = n)
  
  b_diff <- b - b_last
  
  count <- count + 1
  
  if (((t(b_diff) %*% b_diff) / (t(b) %*% b))^0.5 <= tol) break()
}
out <- list()
out$b <- b
out$iter <- nIter
return(out)
}

## Execute the estimations, adjust the upper limit of number of iterations
out <- cbind(la_coef(la_est, n=4, data = la_data)$coef,
             aids_est(maxIter = 1, tol = 1e-5, b_start = la_coef(la_est, n=4, data = la_data)$coef, data = aids_data, n=4)$b,
             aids_est(maxIter = 5, tol = 1e-5, b_start = la_coef(la_est, n=4, data = la_data)$coef, data = aids_data, n=4)$b,
             aids_est(maxIter = 10, tol = 1e-5, b_start = la_coef(la_est, n=4, data = la_data)$coef, data = aids_data, n=4)$b) %>% as.data.frame()

## Get the actual number of runs (which is likely to be lower than our preset upper limit)
colnames(out) <- c("LA-AIDS",paste0("No. of iterations = ",aids_est(maxIter = 1, tol = 1e-5, b_start = la_coef(la_est, n=4, data = la_data)$coef, data = aids_data, n=4)$iter), paste0("No. of iterations = ", aids_est(maxIter = 5, tol = 1e-5, b_start = la_coef(la_est, n=4, data = la_data)$coef, data = aids_data, n=4)$iter),paste0("No. of iterations = ", aids_est(maxIter = 10, tol = 1e-5, b_start = la_coef(la_est, n=4, data = la_data)$coef, data = aids_data, n=4)$iter,' (convergence)'))

## Print, compare with LA results
knitr::kable(out, 
             align = rep('c',ncol(out)),
             booktabs = T) %>%
  column_spec(2:5, width = "8em") %>%
  add_header_above(c(" " = 2, "AIDS" = 3))
  

## ----aids_elas, echo = knitr::is_html_output(), restuls = 'markup'------------
## Function: to calcualte expenditure and price elasticities of demand (Marshallian and Hicksian) ----
aids_elas <- function(coef,data,shares,prices) {
  alpha <- coef[grepl('alpha',names(coef))]
  n <- length(alpha)
  beta <- coef[grepl('beta',names(coef))]
  gamma <- matrix(coef[grepl('gamma',names(coef))], nrow = n)
  
  ones <- rep(1,n)
  
  expenditure <- ones + beta/shares
  
  marshall <- -diag(1,n,n) + gamma/(shares %*% t(ones)) - beta %*% t(ones) * (ones %*% t(alpha) + ones %*% t(gamma %*% log(prices))) / (shares %*% t(ones))
  
  hicks <- marshall + (expenditure %*% t(ones)) * (ones %*% t(shares))
  
  names(expenditure) <- paste0("q_",c(1:n))
  rownames(hicks) <- paste0("q_",c(1:n))
  colnames(hicks) <- paste0("p_",c(1:n))
  rownames(marshall) <- paste0("q_",c(1:n))
  colnames(marshall) <- paste0("p_",c(1:n))
  
  ## Dimensions of the output
  out <- list()
  out$expenditure <- expenditure
  out$hicks <- hicks
  out$marshall <- marshall
  
  ## Print
  cat("Expenditure Elasticities \n")
  print(out$expenditure)
  cat("\n Marshallian Price Elasticities \n")
  print(out$marshall)
  cat("\n Hicksian (compensated) Price Elasticities \n")
  print(out$hicks)
}

## Excute the functions
aids_elas(coef = aids_est(maxIter = 10, tol = 1e-5, b_start = la_coef(la_est, n=4, data = la_data)$coef, data = aids_data, n=4)$b,
          data = aids_data,
          shares = unique(data_long %>% ungroup %>% select(c('Food','w_j')))$w_j,
          prices = unique(data_long %>% ungroup %>% select(c('Food','p_j')))$p_j)

