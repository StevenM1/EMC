ll_mv <- function(pars, data){
  y <- data[,ncol(data)-1]
  X <- as.matrix(data[,1:(ncol(data)-2)])
  betas <- pars[1:(length(pars) -2)]
  sigma <- exp(pars[length(pars) -1])
  rho <- pars[length(pars)]
  order <- toeplitz(0:(length(y) -1))
  V <- rho ** order * (sigma**2)
  
  return(max(mvtnorm::dmvnorm(y, mean=X %*% betas, sigma = V, log = T, checkSymmetry = F), -10000))
}

ll_uv <- function(pars, data){
  y <- data[,ncol(data)-1]
  X <- as.matrix(data[,1:(ncol(data)-2)])
  betas <- pars[1:(length(pars) -1)]
  sigma <- exp(pars[length(pars)])
  return(max(sum(dnorm(y, mean=X %*% betas, sd = sigma, log = T)), -10000))
}

behav_data <- read.table(file = 'bids/behavior.tsv', sep = '\t', header = TRUE)
behav_data <- behav_data %>% rename(subject = pp)
behav_data <- behav_data[behav_data$subject %in% as.numeric(subjects),]
behav_data$response <- ifelse(behav_data$response == "left", 1, 2)
behav_data$stimulus <- ifelse(behav_data$stimulus == "stimleft", 1, 2)
