library(mclust)
# log-sum-exp trick
LogSumExp <- function(vec) {
  max_vec <- max(vec)
  trick <- max_vec + log(sum(exp(vec - max_vec)))
  return(trick)
}

# K should be >= 2
check_K <- function(K) {
  for(x in K) {
    if(x != as.integer(x)) {
      stop("K should be a vector of integer")
    }
  }
  if(min(K) < 2) {
    stop("each element in K should be greater or equal than 2")
  }
}


# initialize Beta
initialize_Beta <- function(K, nG) {
  nOmics <- length(K)
  Beta <- vector(mode = "list", length = nOmics)
  for(i in 1:nOmics) {
    # row represent cluster, column represents variable
    Beta[[i]] <- matrix(runif((nG + 1) * (K[i] - 1), min = -1, max = 1),
                        nrow = K[i] - 1)
  }
  return(Beta)
}

# initialize Mu
initialize_Mu <- function(K, nZ) {
  nOmics <- length(K)
  Mu <- vector(mode = "list", length = nOmics)
  for(i in 1:nOmics) {
    # column represents cluster, row represents variable
    Mu[[i]] <- matrix(runif(K[i] * nZ[i], min = -1, max = 1),
                      nrow = nZ[i])
  }
  return(Mu)
}


# initialize Sigma
initialize_Sigma <- function(K, nZ) {
  nOmics <- length(K)
  Sigma <- vector(mode = "list", length = nOmics)
  for(i in 1:nOmics) {
    # each element is an array of nZ[i] x nZ[i] x K[i]
    Sigma[[i]] <- array(0, dim = c(nZ[i], nZ[i], K[i]))
    for(j in 1:K[i]) {
      Sigma[[i]][, , j] <- diag(nZ[i])
    }
  }
  return(Sigma)
}



# initialize Mu and Sigma
initialize_Mu_Sigma <- function(K, Z, modelNames) {
  nOmics <- length(K)
  Mu <- vector(mode = "list", length = nOmics)
  Sigma <- vector(mode = "list", length = nOmics)
  z <- vector(mode = "list", length = nOmics)
  for(i in 1:nOmics) {
    temp_fit <- Mclust(data = Z[[i]],
                       G = K[i],
                       modelNames = modelNames[i])
    Mu[[i]] <- temp_fit$parameters$mean
    Sigma[[i]] <- temp_fit$parameters$variance$sigma
    z[[i]] <- temp_fit$z
  }
  return(list(Mu = Mu,
              Sigma = Sigma,
              z = z))
}

# initialize Delta
# for normal outcome, Delta is a list, with beta + sigma
# for binary outcome, Delta is a vector
initialize_Delta <- function(K, nCoY = 0, family = c("gaussian", "binomial"),
                             z, Y) {
  family <- match.arg(family)
  if(family == "gaussian") {
    
    # if 2 omics layers
    if(length(K) == 2) {
      r_matrix <- cbind(z[[1]], z[[2]])
      r_fit <- r_matrix[, -c(1, K[1] + 1)]
      fit <- lm(Y ~ r_fit)
      mu <- as.numeric(coef(fit))
      sd <- sd(resid(fit))
      x <- list(mu = mu,
                sd = sd,
                K = K)
    }
    
    # if 3 omics layers
    if(length(K) == 3) {
      r_matrix <- cbind(z[[1]], z[[2]], z[[3]])
      r_fit <- r_matrix[, -c(1, K[1] + 1, K[1] + K[2] + 1)]
      fit <- lm(Y ~ r_fit)
      mu <- as.numeric(coef(fit))
      sd <- sd(resid(fit))
      x <- list(mu = mu,
                sd = sd,
                K = K)
    }
    
    # if 4 omics layers
    if(length(K) == 4) {
      r_matrix <- cbind(z[[1]], z[[2]], z[[3]], z[[4]])
      r_fit <- r_matrix[, -c(1, K[1] + 1, K[1] + K[2] + 1, K[1] + K[2] + K[3] + 1)]
      fit <- lm(Y ~ r_fit)
      mu <- as.numeric(coef(fit))
      sd <- sd(resid(fit))
      x <- list(mu = mu,
                sd = sd,
                K = K)
    }
    
    # if 5 omics layers
    if(length(K) == 5) {
      r_matrix <- cbind(z[[1]], z[[2]], z[[3]], z[[4]], z[[5]])
      r_fit <- r_matrix[, -c(1, K[1] + 1, K[1] + K[2] + 1, K[1] + K[2] + K[3] + 1, K[1] + K[2] + K[3] + K[4] + 1)]
      fit <- lm(Y ~ r_fit)
      mu <- as.numeric(coef(fit))
      sd <- sd(resid(fit))
      x <- list(mu = mu,
                sd = sd,
                K = K)
    }
    
  }
  
  
  if(family == "binomial") {
    
    # if 2 omics layers
    if(length(K) == 2) {
      r_matrix <- cbind(z[[1]], z[[2]])
      r_fit <- r_matrix[, -c(1, K[1] + 1)]
      fit <- glm(Y ~ r_fit, family = "binomial")
      b <- as.numeric(coef(fit))
      b_array <- vec_to_array(K = K, mu = b)
      p <- 1 / (1 + exp(-b_array))
      x <- list(mu = p,
                K = K)
    }
    
    
    # if 3 omics layers
    if(length(K) == 3) {
      r_matrix <- cbind(z[[1]], z[[2]], z[[3]])
      r_fit <- r_matrix[, -c(1, K[1] + 1, K[1] + K[2] + 1)]
      fit <- glm(Y ~ r_fit, family = "binomial")
      b <- as.numeric(coef(fit))
      b_array <- vec_to_array(K = K, mu = b)
      p <- 1 / (1 + exp(-b_array))
      x <- list(mu = p,
                K = K)
    }
    
    # if 4 omics layers
    if(length(K) == 4) {
      r_matrix <- cbind(z[[1]], z[[2]], z[[3]], z[[4]])
      r_fit <- r_matrix[, -c(1, K[1] + 1, K[1] + K[2] + 1, K[1] + K[2] + K[3] + 1)]
      fit <- glm(Y ~ r_fit, family = "binomial")
      b <- as.numeric(coef(fit))
      b_array <- vec_to_array(K = K, mu = b)
      p <- 1 / (1 + exp(-b_array))
      x <- list(mu = p,
                K = K)
    }
    
    # if 5 omics layers
    if(length(K) == 4) {
      r_matrix <- cbind(z[[1]], z[[2]], z[[3]], z[[4]], z[[4]])
      r_fit <- r_matrix[, -c(1, K[1] + 1, K[1] + K[2] + 1, K[1] + K[2] + K[3] + 1, K[1] + K[2] + K[3] + K[4] + 1)]
      fit <- glm(Y ~ r_fit, family = "binomial")
      b <- as.numeric(coef(fit))
      b_array <- vec_to_array(K = K, mu = b)
      p <- 1 / (1 + exp(-b_array))
      x <- list(mu = p,
                K = K)
    }
    
    
  }
  return(x)
}


# indicator function
indicator <- function(x) {
  m <- 0
  if(x > 1) {
    m <- 1
  }
  return(m)
}

# transform mu to an array, each element corresponds to a mean for a combination
# of clusters
vec_to_array <- function(K, mu) {
  res <- array(data = rep(0, prod(K)),
               dim = K)
  # if nK = 2, transform mu to a matrix
  if(length(K) == 2) {
    for(i in 1:K[1]) {
      for(j in 1:K[2]) {
        res[i, j] <- mu[1] + indicator(i) * mu[i] + indicator(j) * mu[K[1] + j - 1]
      }
    }
  }
  
  # if nK = 3, transform mu to a 3d array
  if(length(K) == 3) {
    for(i in 1:K[1]) {
      for(j in 1:K[2]) {
        for (k in 1:K[3]) {
          res[i, j, k] <- mu[1] + indicator(i) * mu[i] + indicator(j) * mu[K[1] + j - 1] + indicator(k) * mu[K[1] + K[2] + k - 2]
        }
      }
    }
  }
  
  # if nK = 4, transform mu to a 4d array
  if(length(K) == 4) {
    for(i in 1:K[1]) {
      for(j in 1:K[2]) {
        for (k in 1:K[3]) {
          for(l in 1:K[4]) {
            res[i, j, k, l] <- mu[1] + indicator(i) * mu[i] + indicator(j) * mu[K[1] + j - 1] + indicator(k) * mu[K[1] + K[2] + k - 2] + indicator(l) * mu[K[1] + K[2] + K[3] + l - 3]
          }
        }
      }
    }
  }
  
  # if nK = 5, transform mu to a 5d array
  if(length(K) == 5) {
    for(i in 1:K[1]) {
      for(j in 1:K[2]) {
        for (k in 1:K[3]) {
          for(l in 1:K[4]) {
            for(m in 1:K[5]) {
              res[i, j, k, l] <- mu[1] + indicator(i) * mu[i] + indicator(j) * mu[K[1] + j - 1] + indicator(k) * mu[K[1] + K[2] + k - 2] + indicator(l) * mu[K[1] + K[2] + K[3] + l - 3] + indicator(m) * mu[K[1] + K[2] + K[3] + K[4] + m - 4]
            }
          }
        }
      }
    }
  }
  
  return(res)
}

# rearrange cluster order
#
# for continuous outcome - use the cluster combination corresponding to smallest
# mean as the reference cluster
get_ref_cluster <- function(Delta) {
  K <- Delta$K
  mu <- Delta$mu
  mu_matrix <- vec_to_array(K = K, mu = mu)
  ref_index <- which(mu_matrix == min(mu_matrix))
  ref <- arrayInd(ref_index, .dim = K)
  return(ref)
}


# re-arrange parameters for Delta
reorder_Delta <- function(ref, Delta) {
  K <- Delta$K
  mu_matrix <- vec_to_array(K = K, mu = Delta$mu)
  mu <- mu_matrix[ref]
  
  # if 2 omics layers
  if(length(K) == 2) {
    # reorder K1
    mu_k1 <- rep(0, K[1])
    for(i in 1:K[1]) {
      mu_k1[i] <- mu_matrix[i, 1] - mu_matrix[ref[1], 1]
    }
    mu_k1_sort <- sort(mu_k1)
    mu <- c(mu, mu_k1_sort[mu_k1_sort != 0])
    # order of re-arranged cluster for omics 1
    k1 <- order(mu_k1)
    k1_order <- c(ref[1], k1[k1 != ref[1]])
    
    # reorder K2
    mu_k2 <- rep(0, K[2])
    for(i in 1:K[2]) {
      mu_k2[i] <- mu_matrix[1, i] - mu_matrix[1, ref[2]]
    }
    mu_k2_sort <- sort(mu_k2)
    mu <- c(mu, mu_k2_sort[mu_k2_sort != 0])
    # order of re-arranged cluster for omics 2
    k2 <- order(mu_k2)
    k2_order <- c(ref[2], k2[k2 != ref[2]])
    
    K_order <- list(K1 = k1_order,
                    K2 = k2_order)
  }
  
  
  # if 3 omics layers
  if(length(K) == 3) {
    # reorder K1
    mu_k1 <- rep(0, K[1])
    for(i in 1:K[1]) {
      mu_k1[i] <- mu_matrix[i, 1, 1] - mu_matrix[ref[1], 1, 1]
    }
    mu_k1_sort <- sort(mu_k1)
    mu <- c(mu, mu_k1_sort[mu_k1_sort != 0])
    # order of re-arranged cluster for omics 1
    k1 <- order(mu_k1)
    k1_order <- c(ref[1], k1[k1 != ref[1]])
    
    # reorder K2
    mu_k2 <- rep(0, K[2])
    for(i in 1:K[2]) {
      mu_k2[i] <- mu_matrix[1, i, 1] - mu_matrix[1, ref[2], 1]
    }
    mu_k2_sort <- sort(mu_k2)
    mu <- c(mu, mu_k2_sort[mu_k2_sort != 0])
    # order of re-arranged cluster for omics 2
    k2 <- order(mu_k2)
    k2_order <- c(ref[2], k2[k2 != ref[2]])
    
    # reorder K3
    mu_k3 <- rep(0, K[3])
    for(i in 1:K[3]) {
      mu_k3[i] <- mu_matrix[1, 1, i] - mu_matrix[1, 1, ref[3]]
    }
    mu_k3_sort <- sort(mu_k3)
    mu <- c(mu, mu_k3_sort[mu_k3_sort != 0])
    # order of re-arranged cluster for omics 3
    k3 <- order(mu_k3)
    k3_order <- c(ref[3], k3[k3 != ref[3]])
    
    
    K_order <- list(K1 = k1_order,
                    K2 = k2_order,
                    K3 = k3_order)
  }
  
  
  # if 4 omics layers
  
  # if 5 omics layers
  if(length(K) == 5) {
    # reorder K1
    mu_k1 <- rep(0, K[1])
    for(i in 1:K[1]) {
      mu_k1[i] <- mu_matrix[i, 1, 1, 1, 1] - mu_matrix[ref[1], 1, 1, 1, 1]
    }
    mu_k1_sort <- sort(mu_k1)
    mu <- c(mu, mu_k1_sort[mu_k1_sort != 0])
    # order of re-arranged cluster for omics 1
    k1 <- order(mu_k1)
    k1_order <- c(ref[1], k1[k1 != ref[1]])
    
    # reorder K2
    mu_k2 <- rep(0, K[2])
    for(i in 1:K[2]) {
      mu_k2[i] <- mu_matrix[1, i, 1, 1, 1] - mu_matrix[1, ref[2], 1, 1, 1]
    }
    mu_k2_sort <- sort(mu_k2)
    mu <- c(mu, mu_k2_sort[mu_k2_sort != 0])
    # order of re-arranged cluster for omics 2
    k2 <- order(mu_k2)
    k2_order <- c(ref[2], k2[k2 != ref[2]])
    
    # reorder K3
    mu_k3 <- rep(0, K[3])
    for(i in 1:K[3]) {
      mu_k3[i] <- mu_matrix[1, 1, i, 1, 1] - mu_matrix[1, 1, ref[3], 1, 1]
    }
    mu_k3_sort <- sort(mu_k3)
    mu <- c(mu, mu_k3_sort[mu_k3_sort != 0])
    # order of re-arranged cluster for omics 3
    k3 <- order(mu_k3)
    k3_order <- c(ref[3], k3[k3 != ref[3]])
    
    # reorder K4
    mu_k4 <- rep(0, K[4])
    for(i in 1:K[4]) {
      mu_k4[i] <- mu_matrix[1, 1, 1, i, 1] - mu_matrix[1, 1, 1, ref[4], 1]
    }
    mu_k4_sort <- sort(mu_k4)
    mu <- c(mu, mu_k4_sort[mu_k4_sort != 0])
    # order of re-arranged cluster for omics 4
    k4 <- order(mu_k4)
    k4_order <- c(ref[4], k4[k4 != ref[4]])
    
    # reorder K5
    mu_k5 <- rep(0, K[5])
    for(i in 1:K[5]) {
      mu_k5[i] <- mu_matrix[1, 1, 1, 1, i] - mu_matrix[1, 1, 1, 1, ref[5]]
    }
    mu_k5_sort <- sort(mu_k5)
    mu <- c(mu, mu_k5_sort[mu_k5_sort != 0])
    # order of re-arranged cluster for omics 5
    k5 <- order(mu_k5)
    k5_order <- c(ref[5], k5[k5 != ref[5]])
    
    
    K_order <- list(K1 = k1_order,
                    K2 = k2_order,
                    K3 = k3_order,
                    K4 = k4_order,
                    K5 = k5_order)
  }
  
  
  Delta$mu <- mu
  return(list(Delta = Delta,
              K_order = K_order))
}



reorder_Mu_Sigma <- function(Mu_Sigma, K_order) {
  for(i in 1:length(K_order)) {
    temp_Mu <- Mu_Sigma$Mu[[i]]
    temp_Sigma <- Mu_Sigma$Sigma[[i]]
    # reorder Mu
    Mu_Sigma$Mu[[i]] <- temp_Mu[, K_order[[i]]]
    Mu_Sigma$Sigma[[i]] <- temp_Sigma[, , K_order[[i]]]
  }
  return(Mu_Sigma)
}



reorder_Beta <- function(Beta, K_order) {
  for(i in 1:length(K_order)) {
    temp_Beta <- Beta[[i]]
    temp_Beta <- rbind(rep(0, ncol(temp_Beta)),
                       temp_Beta)
    temp_Beta_reorder <- temp_Beta[K_order[[i]], ]
    ref <- temp_Beta_reorder[1, ]
    for(j in 1:nrow(temp_Beta_reorder)) {
      temp_Beta_reorder[j, ] <- temp_Beta_reorder[j, ] - ref
    }
    Beta[[i]] <- temp_Beta_reorder[-1, ]
  }
  
  return(Beta)
}



reorder_z <- function(z, K_order) {
  if(length(K_order) == 2) {
    z <- z[K_order[[1]], K_order[[2]], ]
  }
  return(z)
}


#' function to reorder all model parameters
#'
#' @param model A model returned by EM_lucid
#'
#' @return A LUCID model reordered by effect size of outcome
#' @export
#'
reorder_lucid <- function(model) {
  ref <- get_ref_cluster(Delta = model$res_Delta$Delta)
  r_Delta <- reorder_Delta(ref = ref,
                           Delta = model$res_Delta$Delta)
  r_Mu_Sigma <- reorder_Mu_Sigma(model$res_Mu_Sigma,
                                 K_order = r_Delta$K_order)
  r_Beta <- reorder_Beta(Beta = model$res_Beta$Beta,
                         K_order = r_Delta$K_order)
  model$res_Delta$Delta <- r_Delta
  model$res_Mu_Sigma$Mu <- r_Mu_Sigma$Mu
  model$res_Mu_Sigma$Sigma <- r_Mu_Sigma$Sigma
  model$res_Beta$Beta <- r_Beta
  model$z <- reorder_z(model$z, K_order = r_Delta$K_order)
  return(model)
}



# Reorder lucid M 

# ----------- reorder LUCIDusM model by specifying reference cluster ------------
# note: only works for K = 2 in each omic layer
# reference = c(1,1,2)
# lucidus_fit <- fit_reordered
reorder_lucidM <- function(lucidus_fit,
                           reference = NULL) {
  if(is.null(reference)) {
    warning("no reference specified, return the original model")
    return(lucidus_fit)
  }
  
  n_omic <- length(reference)
  
  # reorder beta
  GtoX <- lucidus_fit$res_Beta$Beta
  lucidus_fit$res_Beta$Beta <- lapply(1:n_omic, function(i) {
    (-1)^(reference[i] - 1) * GtoX[[i]] # if reference = 1, no changes; 
    # if reference = 2, flip the reference and negate the estimates
  })
  # reorder mu
  XtoZ <- lucidus_fit$res_Mu_Sigma$Mu
  lucidus_fit$res_Mu_Sigma$Mu <- lapply(1:n_omic, function(i) {
    x <- c(1, 2) # order of clusters
    if(reference[i] == 2) {
      x <- c(2, 1)
      XtoZ[[i]][, x]
    } else{
      XtoZ[[i]][, x]
    }
  }) 
  # reorder gamma
  # XtoY <- lucidus_fit$res_Delta$Delta$mu
  # XtoY[1] <- XtoY[1] + sum(XtoY[-1] * (reference - 1)) # reference level using the new reference
  # XtoY[-1] <- (-1)^(reference - 1) * XtoY[-1] # if reference = 2, flip the estimates
  # lucidus_fit$res_Delta$Delta$mu <- XtoY
  XtoY_rfit <- lucidus_fit$res_Delta$fit$coefficients
  XtoY_rfit[1] <- XtoY_rfit[1] + sum(XtoY_rfit[-1] * (reference - 1)) # reference level using the new reference
  XtoY_rfit[-1] <- (-1)^(reference - 1) * XtoY_rfit[-1] # if reference = 2, flip the estimates
  lucidus_fit$res_Delta$fit$coefficients <- XtoY_rfit
  
  # return the object using the new reference
  return(lucidus_fit)
}




# function to calculate BIC
cal_bic <- function(model) {
  nOmics <- length(model$K)
  Beta <- model$res_Beta$Beta
  Mu <- model$res_Mu_Sigma$Mu
  Sigma <- model$res_Mu_Sigma$Sigma
  Delta <- model$res_Delta$Delta
  
  # calculate number of parameters
  n_Beta <- sum(sapply(1:nOmics, function(i) {
    length(Beta[[i]])
  }))
  n_Mu <- sum(sapply(1:nOmics, function(i) {
    length(Mu[[i]])
  }))
  n_Sigma <- sum(sapply(1:nOmics, function(i) {
    length(Sigma[[i]])
  }))
  n_Delta <- length(Delta$mu) + 1
  n_par <- n_Beta + n_Mu + n_Sigma + n_Delta
  
  # calculate BIC
  bic <- -2 * model$loglik + n_par * log(model$N)
  return(bic)
}


predict_lucid <- function(model,
                          G = NULL,
                          Z = NULL,
                          Y = NULL) {
  K <- model$K
  nOmics <- length(K)
  N <- model$N
  # initialize container for predicted value
  pred_X <- vector(mode = "list", length = nOmics)
  pred_z <- vector(mode = "list", length = nOmics)
  # prediction of X and Y based on fitted data
  if(all(is.null(G), is.null(Z), is.null(Y))) {
    r <- model$z
    # 1 - prediction for X
    for (i in 1:nOmics) {
      r_margin <- t(sapply(1:N, function(j) {
        marginSums(r[, , j], margin = i)
      }))
      pred_X[[i]] <- map(r_margin)
      pred_z[[i]] <- r_margin
    }
    # 2 - prediction for Y
    pred_Y <- predict(model$res_Delta$fit)
  } else { # prediction of X and Y based on new data
    # predict X and Y based on G
    
    # predict X and Y based on Z
    
    # predict X and Y based on X and G
    
    # predict X and Y based
    
  }
  
  return(list(pred_X = pred_X,
              pred_z = pred_z,
              pred_Y = pred_Y))
}


Mstep_GtoX <- function(G, r, K, N) {
  nOmics <- length(K)
  # store multinomial logistic regression model with corresponding coefficients
  fit <- vector(mode = "list", length = nOmics)
  Beta <- vector(mode = "list", length = nOmics)
  
  # if 2 omics layers
  if(nOmics == 2) {
    for(i in 1:nOmics) {
      r_margin <- t(sapply(1:N, function(j) {
        marginSums(r[, , j], margin = i)
      }))
      invisible(capture.output(temp_fit <- nnet::multinom(r_margin ~ G)))
      fit[[i]] <- temp_fit
      Beta[[i]] <- coef(temp_fit)
    }
  }
  
  # if 3 omics layers
  if(nOmics == 3) {
    for(i in 1:nOmics) {
      r_margin <- t(sapply(1:N, function(j) {
        marginSums(r[, , , j], margin = i)
      }))
      invisible(capture.output(temp_fit <- nnet::multinom(r_margin ~ G)))
      fit[[i]] <- temp_fit
      Beta[[i]] <- coef(temp_fit)
    }
  }
  
  
  # if 4 omics layers
  if(nOmics == 4) {
    for(i in 1:nOmics) {
      r_margin <- t(sapply(1:N, function(j) {
        marginSums(r[, , , , j], margin = i)
      }))
      invisible(capture.output(temp_fit <- nnet::multinom(r_margin ~ G)))
      fit[[i]] <- temp_fit
      Beta[[i]] <- coef(temp_fit)
    }
  }
  
  
  # if 5 omics layers
  if(nOmics == 5) {
    for(i in 1:nOmics) {
      r_margin <- t(sapply(1:N, function(j) {
        marginSums(r[, , , , , j], margin = i)
      }))
      invisible(capture.output(temp_fit <- nnet::multinom(r_margin ~ G)))
      fit[[i]] <- temp_fit
      Beta[[i]] <- coef(temp_fit)
    }
  }
  
  return(list(fit = fit,
              Beta = Beta))
}

Mstep_XtoZ <- function(Z, r, K, modelNames, N) {
  nOmics <- length(K)
  # store GMM model with corresponding model
  fit <- vector(mode = "list", length = nOmics)
  Mu <- vector(mode = "list", length = nOmics)
  Sigma <- vector(mode = "list", length = nOmics)
  
  # if 2 omics data
  if(nOmics == 2) {
    for(i in 1:nOmics) {
      r_margin <- t(sapply(1:N, function(j) {
        marginSums(r[, , j], margin = i)
      }))
      r_margin <- round(r_margin, digits = 8)
      temp_fit <- mstep(data = Z[[i]],
                        G = K[i],
                        z = r_margin,
                        modelName = modelNames[i])
      fit[[i]] <- temp_fit
      Mu[[i]] <- temp_fit$parameters$mean
      Sigma[[i]] <- temp_fit$parameters$variance$sigma
    }
  }
  
  
  # if 3 omics data
  if(nOmics == 3) {
    for(i in 1:nOmics) {
      r_margin <- t(sapply(1:N, function(j) {
        marginSums(r[, , , j], margin = i)
      }))
      r_margin <- round(r_margin, digits = 8)
      temp_fit <- mstep(data = Z[[i]],
                        G = K[i],
                        z = r_margin,
                        modelName = modelNames[i])
      fit[[i]] <- temp_fit
      Mu[[i]] <- temp_fit$parameters$mean
      Sigma[[i]] <- temp_fit$parameters$variance$sigma
    }
  }
  
  # if 4 omics data
  if(nOmics == 4) {
    for(i in 1:nOmics) {
      r_margin <- t(sapply(1:N, function(j) {
        marginSums(r[, , , , j], margin = i)
      }))
      r_margin <- round(r_margin, digits = 8)
      temp_fit <- mstep(data = Z[[i]],
                        G = K[i],
                        z = r_margin,
                        modelName = modelNames[i])
      fit[[i]] <- temp_fit
      Mu[[i]] <- temp_fit$parameters$mean
      Sigma[[i]] <- temp_fit$parameters$variance$sigma
    }
  }
  
  # if 5 omics data
  if(nOmics == 5) {
    for(i in 1:nOmics) {
      r_margin <- t(sapply(1:N, function(j) {
        marginSums(r[, , , , , j], margin = i)
      }))
      r_margin <- round(r_margin, digits = 8)
      temp_fit <- mstep(data = Z[[i]],
                        G = K[i],
                        z = r_margin,
                        modelName = modelNames[i])
      fit[[i]] <- temp_fit
      Mu[[i]] <- temp_fit$parameters$mean
      Sigma[[i]] <- temp_fit$parameters$variance$sigma
    }
  }
  
  
  return(list(fit = fit,
              Mu = Mu,
              Sigma = Sigma))
}





Mstep_XtoY <- function(Y, r, K, N, family) {
  
  
  # if 2 omics layers
  if(length(K) == 2) {
    r_matrix <- t(sapply(1:N, function(i) {
      c(rowSums(r[, , i]), colSums(r[, , i]))
    }))
    r_fit <- r_matrix[, -c(1, K[1] + 1)]
    
    if(family == "gaussian") {
      fit <- lm(Y ~ r_fit)
      mu <- as.numeric(coef(fit))
      sd <- sd(resid(fit))
    }
    
    if(family == "binomial") {
      fit <- glm(Y ~ r_fit, family = "binomial")
      mu <- as.numeric(coef(fit))
      mu_array <- vec_to_array(K = K, mu = mu)
      p <- exp(mu_array) / sum(exp(mu_array))
      sd <- NULL
      # fit <- NULL
      # p <- matrix(rep(0, prod(K)),
      #             nrow = K[1],
      #             ncol = K[2])
      # for(i in 1:K[1]) {
      #   for(j in 1:K[2]) {
      #     p[i, j] <- sum(r[i, j, ] * Y) / sum(r[i, j, ])
      #   }
      # }
      fit <- fit
      mu <- p
      sd <- NULL
    }
    
    if(any(is.na(mu))) {
      na_index <- which(is.na(mu))
      if(na_index <= K[1]) {
        stop("no cluster strucutre is defined for Z1")
      } else{
        stop("no cluster structure is defined for Z2")
      }
    }
  }
  
  
  # if 3 omics layers
  if(length(K) == 3) {
    r_matrix <- t(sapply(1:N, function(i) {
      c(marginSums(r[, , , i], margin = 1),
        marginSums(r[, , , i], margin = 2),
        marginSums(r[, , , i], margin = 3))
    }))
    r_fit <- r_matrix[, -c(1, K[1] + 1, K[1] + K[2] + 1)]
    
    if(family == "gaussian") {
      fit <- lm(Y ~ r_fit)
      mu <- as.numeric(coef(fit))
      sd <- sd(resid(fit))
    }
    
    if(family == "binomial") {
      fit <- glm(Y ~ r_fit, family = "binomial")
      mu <- as.numeric(coef(fit))
      mu_array <- vec_to_array(K = K, mu = mu)
      p <- exp(mu_array) / sum(exp(mu_array))
      sd <- NULL
      fit <- fit
      mu <- p
      sd <- NULL
    }
    
    if(any(is.na(mu))) {
      na_index <- which(is.na(mu))
      if(na_index <= K[1]) {
        stop("no cluster strucutre is defined for Z1")
      } else{
        stop("no cluster structure is defined for Z2")
      }
    }
  }
  
  
  # if 4 omics layers
  if(length(K) == 4) {
    r_matrix <- t(sapply(1:N, function(i) {
      c(marginSums(r[, , , , i], margin = 1),
        marginSums(r[, , , , i], margin = 2),
        marginSums(r[, , , , i], margin = 3),
        marginSums(r[, , , , i], margin = 4))
    }))
    r_fit <- r_matrix[, -c(1, K[1] + 1, K[1] + K[2] + 1, K[1] + K[2] + K[3] + 1)]
    
    if(family == "gaussian") {
      fit <- lm(Y ~ r_fit)
      mu <- as.numeric(coef(fit))
      sd <- sd(resid(fit))
    }
    
    if(family == "binomial") {
      fit <- glm(Y ~ r_fit, family = "binomial")
      mu <- as.numeric(coef(fit))
      mu_array <- vec_to_array(K = K, mu = mu)
      p <- exp(mu_array) / sum(exp(mu_array))
      sd <- NULL
      fit <- fit
      mu <- p
      sd <- NULL
    }
    
    if(any(is.na(mu))) {
      na_index <- which(is.na(mu))
      if(na_index <= K[1]) {
        stop("no cluster strucutre is defined for Z1")
      } else{
        stop("no cluster structure is defined for Z2")
      }
    }
  }
  
  
  # if 5 omics layers
  if(length(K) == 5) {
    r_matrix <- t(sapply(1:N, function(i) {
      c(marginSums(r[, , , , , i], margin = 1),
        marginSums(r[, , , , , i], margin = 2),
        marginSums(r[, , , , , i], margin = 3),
        marginSums(r[, , , , , i], margin = 4),
        marginSums(r[, , , , , i], margin = 5))
    }))
    r_fit <- r_matrix[, -c(1,
                           K[1] + 1,
                           K[1] + K[2] + 1,
                           K[1] + K[2] + K[3] + 1,
                           K[1] + K[2] + K[3] + K[4] + 1)]
    
    if(family == "gaussian") {
      fit <- lm(Y ~ r_fit)
      mu <- as.numeric(coef(fit))
      sd <- sd(resid(fit))
    }
    
    if(family == "binomial") {
      fit <- glm(Y ~ r_fit, family = "binomial")
      mu <- as.numeric(coef(fit))
      mu_array <- vec_to_array(K = K, mu = mu)
      p <- exp(mu_array) / sum(exp(mu_array))
      sd <- NULL
      fit <- fit
      mu <- p
      sd <- NULL
    }
    
    if(any(is.na(mu))) {
      na_index <- which(is.na(mu))
      if(na_index <= K[1]) {
        stop("no cluster strucutre is defined for Z1")
      } else{
        stop("no cluster structure is defined for Z2")
      }
    }
  }
  
  
  Delta <- list(mu = mu,
                sd = sd,
                K = K)
  return(list(fit = fit,
              Delta = Delta))
}
# use bic to conduct model selection
tune_lucid <- function(K_list, G, Z, Y, ...) {
  # change K_list into a matrix
  K_matrix <- as.matrix(expand.grid(K_list))
  if(min(K_matrix) < 2) {
    stop("minimum K should be 2")
  }
  bic <- rep(0, nrow(K_matrix))
  model_list <- vector(mode = "list", length = nrow(K_matrix))
  for(i in 1:nrow(K_matrix)) {
    model_list[[i]] <- EM_lucid(G = G, Z = Z, Y = Y,
                                K = K_matrix[i, ],
                                ...)
    bic[i] <- cal_bic(model_list[[i]])
  }
  model_opt_index <- which(bic == min(bic))
  K_matrix <- cbind(K_matrix, bic)
  return(list(tune_K = K_matrix,
              model_list = model_list,
              model_opt = model_list[[model_opt_index]]))
}

f_GtoX <- function(G, Beta_matrix) {
  N <- nrow(G)
  Beta_matrix <- rbind(rep(0, ncol(Beta_matrix)),
                       Beta_matrix)
  xb <- cbind(rep(1, N), G) %*% t(Beta_matrix)
  xb_LSE <- apply(xb, 1, LogSumExp)
  return(xb - xb_LSE)
}

f_XtoZ <- function(Z, Mu_matrix, Sigma_matrix) {
  N <- nrow(Z)
  K <- ncol(Mu_matrix)
  XtoZ <- matrix(rep(0, N * K), nrow = N)
  for (i in 1:K) {
    XtoZ[, i] <- mclust::dmvnorm(data = Z,
                                 mean = Mu_matrix[, i],
                                 sigma = Sigma_matrix[, , i],
                                 log = TRUE)
  }
  return(XtoZ)
}


#' @title log f(Y|X)
#'
#' @description Calculate the log likelihood of outcome Y given all latent variables
#' X
#'
#' @param Y an N by 1 matrix
#' @param Delta a list with parameters related to outcome
#' @param family a string, either gaussian or binomial
#'
#' @return for 2 omics layers, an K1 by K2 by N matrix; for 3 omics layers, an
#' K1 by K2 by K3 by N matrix
#'
#' @examples
#' Y <- matrix(rnorm(100), nrow = 100)
#' delta <- list(mu = rnorm(4), sd = 1, K = c(2, 3))
#' f_XtoY(Y = Y, Delta = delta, family = "gaussian")
f_XtoY <- function(Y, Delta, family) {
  
  if(!is.matrix(Y)) {
    stop("Y should be a matrix")
  }
  N <- nrow(Y)
  K <- Delta$K
  XtoY <- array(rep(0, prod(K) * N), dim = c(K, N))
  
  # if the outcome is continuous
  if(family == "gaussian") {
    
    # if 2 omics layers
    if(length(K) == 2) {
      mu <- vec_to_array(K = K, mu = Delta$mu)
      for(i in 1:K[1]) {
        for(j in 1:K[2]) {
          XtoY[i, j, ] <- dnorm(Y, mean = mu[i, j], sd = Delta$sd, log = TRUE)
        }
      }
    }
    
    # if 3 omics layers
    if(length(K) == 3) {
      mu <- vec_to_array(K = K, mu = Delta$mu)
      for (i in 1:K[1]) {
        for(j in 1:K[2]) {
          for(k in 1:K[3]) {
            XtoY[i, j, k, ] <- dnorm(Y, mean = mu[i, j, k], sd = Delta$sd, log = TRUE)
          }
        }
      }
    }
    
    # if 4 omics layers
    if(length(K) == 4) {
      mu <- vec_to_array(K = K, mu = Delta$mu)
      for (i in 1:K[1]) {
        for(j in 1:K[2]) {
          for(k in 1:K[3]) {
            for(l in 1:K[4]) {
              XtoY[i, j, k, l, ] <- dnorm(Y, mean = mu[i, j, k, l], sd = Delta$sd, log = TRUE)
            }
          }
        }
      }
    }
    
    # if 5 omics layers
    if(length(K) == 5) {
      mu <- vec_to_array(K = K, mu = Delta$mu)
      for (i in 1:K[1]) {
        for(j in 1:K[2]) {
          for(k in 1:K[3]) {
            for(l in 1:K[4]) {
              for(m in 1:K[5]) {
                XtoY[i, j, k, l, m, ] <- dnorm(Y, mean = mu[i, j, k, l, m], sd = Delta$sd, log = TRUE)
              }
            }
          }
        }
      }
    }
    
  }
  
  # if the outcome is binary
  if(family == "binomial") {
    p <- Delta$mu
    
    # if 2 omics layers
    if(length(K) == 2) {
      for(i in 1:K[1]) {
        for(j in 1:K[2]) {
          XtoY[i, j, ] <- dbinom(Y, size = 1, p = p[i, j], log = TRUE)
        }
      }
    }
    
    # if 3 omics layers
    if(length(K) == 3) {
      for(i in 1:K[1]) {
        for(j in 1:K[2]) {
          for(k in 1:K[3]) {
            XtoY[i, j, k, ] <- dbinom(Y, size = 1, p = p[i, j, k], log = TRUE)
          }
        }
      }
    }
    
    # if 4 omics layers
    if(length(K) == 4) {
      for(i in 1:K[1]) {
        for(j in 1:K[2]) {
          for(k in 1:K[3]) {
            for(l in 1:K[4]) {
              XtoY[i, j, k, l, ] <- dbinom(Y, size = 1, p = p[i, j, k, l], log = TRUE)
            }
          }
        }
      }
    }
    
    # if 5 omics layers
    if(length(K) == 5) {
      for(i in 1:K[1]) {
        for(j in 1:K[2]) {
          for(k in 1:K[3]) {
            for(l in 1:K[4]) {
              for(m in 1:K[5]) {
                XtoY[i, j, k, l, m, ] <- dbinom(Y, size = 1, p = p[i, j, k, l, m], log = TRUE)
              }
            }
          }
        }
      }
    }
    
  }
  return(XtoY)
}


# Estep
Estep <- function(G, Z, Y, Beta, Mu, Sigma, Delta, family, useY) {
  N <- nrow(Y)
  K <- Delta$K
  res <- array(rep(0, prod(K) * N),
               dim = c(K, N))
  
  # E step for 2 omics data
  if(length(K) == 2) {
    f1 <- lapply(1:2, function(i) {
      f_GtoX(G = G, Beta_matrix = Beta[[i]])
    })
    f2 <- lapply(1:2, function(i) {
      f_XtoZ(Z = Z[[i]], Mu_matrix = Mu[[i]], Sigma_matrix = Sigma[[i]])
    })
    if(useY) {
      f3 <- f_XtoY(Y = Y, Delta = Delta, family = family)
    }
    
    for(i in 1:K[1]) {
      for(j in 1:K[2]) {
        res[i, j, ] <- f1[[1]][, i] + f1[[2]][, j] + f2[[1]][, i] + f2[[2]][, j]
        if(useY) {
          res[i, j, ] <- res[i, j, ] + f3[i, j, ]
        }
      }
    }
  }
  
  # E step for 3 omics data
  if(length(K) == 3) {
    f1 <- lapply(1:3, function(i) {
      f_GtoX(G = G, Beta_matrix = Beta[[i]])
    })
    f2 <- lapply(1:3, function(i) {
      f_XtoZ(Z = Z[[i]], Mu_matrix = Mu[[i]], Sigma_matrix = Sigma[[i]])
    })
    if(useY) {
      f3 <- f_XtoY(Y = Y, Delta = Delta, family = family)
    }
    
    for(i in 1:K[1]) {
      for(j in 1:K[2]) {
        for(k in 1:K[3]) {
          res[i, j, k, ] <- f1[[1]][, i] + f1[[2]][, j] + f1[[3]][, k] + f2[[1]][, i] + f2[[2]][, j] + f2[[3]][, k]
          if(useY) {
            res[i, j, k, ] <- res[i, j, k, ] + f3[i, j, k, ]
          }
        }
      }
    }
  }
  
  
  # E step for 4 omics layers
  if(length(K) == 4) {
    f1 <- lapply(1:4, function(i) {
      f_GtoX(G = G, Beta_matrix = Beta[[i]])
    })
    f2 <- lapply(1:4, function(i) {
      f_XtoZ(Z = Z[[i]], Mu_matrix = Mu[[i]], Sigma_matrix = Sigma[[i]])
    })
    if(useY) {
      f3 <- f_XtoY(Y = Y, Delta = Delta, family = family)
    }
    
    for(i in 1:K[1]) {
      for(j in 1:K[2]) {
        for(k in 1:K[3]) {
          for(l in 1:K[4]) {
            res[i, j, k, l, ] <- f1[[1]][, i] + f1[[2]][, j] + f1[[3]][, k] + f1[[4]][, l] + f2[[1]][, i] + f2[[2]][, j] + f2[[3]][, k] + f2[[4]][, l]
            if(useY) {
              res[i, j, k, l, ] <- res[i, j, k, l, ] + f3[i, j, k, l, ]
            }
          }
        }
      }
    }
  }
  
  
  
  # E step for 5 omics layers
  if(length(K) == 5) {
    f1 <- lapply(1:5, function(i) {
      f_GtoX(G = G, Beta_matrix = Beta[[i]])
    })
    f2 <- lapply(1:5, function(i) {
      f_XtoZ(Z = Z[[i]], Mu_matrix = Mu[[i]], Sigma_matrix = Sigma[[i]])
    })
    if(useY) {
      f3 <- f_XtoY(Y = Y, Delta = Delta, family = family)
    }
    
    for(i in 1:K[1]) {
      for(j in 1:K[2]) {
        for(k in 1:K[3]) {
          for(l in 1:K[4]) {
            for(m in 1:K[5]) {
              res[i, j, k, l, ] <- f1[[1]][, i] + f1[[2]][, j] + f1[[3]][, k] + f1[[4]][, l] + f1[[5]][, m] + f2[[1]][, i] + f2[[2]][, j] + f2[[3]][, k] + f2[[4]][, l] + f2[[5]][, m]
              if(useY) {
                res[i, j, k, l, ] <- res[i, j, k, l, m, ] + f3[i, j, k, l, m, ]
              }
            }
          }
        }
      }
    }
  }
  
  return(res)
}



Estep_to_r <- function(Estep_array, K, N) {
  
  # E step for 2 omics layers
  if(length(K) == 2) {
    for(i in 1:N) {
      a <- as.vector(Estep_array[, , i])
      r <- exp(a - LogSumExp(a))
      Estep_array[, , i] <- array(r, dim = K)
    }
  }
  
  # E step for 3 omics layers
  if(length(K) == 3) {
    for(i in 1:N) {
      a <- as.vector(Estep_array[, , , i])
      r <- exp(a - LogSumExp(a))
      Estep_array[, , , i] <- array(r, dim = K)
    }
  }
  
  
  # E step for 4 omics layers
  if(length(K) == 4) {
    for(i in 1:N) {
      a <- as.vector(Estep_array[, , , , i])
      r <- exp(a - LogSumExp(a))
      Estep_array[, , , , i] <- array(r, dim = K)
    }
  }
  
  # E step for 5 omics layers
  if(length(K) == 5) {
    for(i in 1:N) {
      a <- as.vector(Estep_array[, , , , , i])
      r <- exp(a - LogSumExp(a))
      Estep_array[, , , , , i] <- array(r, dim = K)
    }
  }
  
  return(Estep_array)
}



cal_loglik <- function(Estep_array, Estep_r) {
  return(sum(Estep_array * Estep_r))
}


est_lucid_parallel <- function(G, Z, Y, K,
                     modelNames = rep("VVV", length(K)),
                     useY = TRUE,
                     init_par = NULL,
                     tol = 1e-3,
                     max_itr = 1e3,
                     family = c("gaussian", "binomial"),
                     seed = 123) {
  # check data format
  if(!is.matrix(G)) {
    G <- as.matrix(G)
  }
  for(i in 1:length(Z)) {
    if(!is.matrix(Z[[i]])) {
      Z[[i]] <- as.matrix(Z[[i]])
    }
  }
  if(!is.matrix(Y)) {
    Y <- as.matrix(Y)
  }


  # basic setup
  N <- nrow(G)
  nOmics <- length(Z)
  nG <- ncol(G)
  nZ <- as.integer(sapply(Z, ncol))
  family <- match.arg(family)

  # initialize model parameters
  set.seed(seed)
  Mu_Sigma <- initialize_Mu_Sigma(K = K, Z = Z, modelNames = modelNames)
  Mu <- Mu_Sigma$Mu
  Sigma <- Mu_Sigma$Sigma
  Beta <- vector(mode = "list", length = nOmics)
  for(i in 1:nOmics) {
    invisible(capture.output(temp_fit <- nnet::multinom(Mu_Sigma$z[[i]] ~ G)))
    Beta[[i]] <- coef(temp_fit)
  }
  # Beta <- initialize_Beta(K = K, nG = nG)
  Delta <- initialize_Delta(K = K, nCoY = 0, family = family,
                            z = Mu_Sigma$z, Y = Y)
  loglik <- -Inf

  # start EM algorithm
  flag_converge <- FALSE
  itr <- 0
  while(!flag_converge & itr < max_itr) {
    # E-step
    Estep_array <- Estep(G = G, Z = Z, Y = Y,
                         Beta = Beta, Mu = Mu, Sigma = Sigma, Delta = Delta,
                         family = family, useY = useY)
    Estep_r <- Estep_to_r(Estep_array = Estep_array,
                          K = K,
                          N = N)

    # M-step
    res_Beta <- Mstep_GtoX(G = G, r = Estep_r, K = K, N = N)
    res_Mu_Sigma <- Mstep_XtoZ(Z = Z, r = Estep_r, K = K,
                               modelNames = modelNames, N = N)
    if(useY) {
      res_Delta <- Mstep_XtoY(Y = Y, r = Estep_r, K = K, N = N,
                              family = family)
    }

    # update parameters
    Beta <- res_Beta$Beta
    Mu <- res_Mu_Sigma$Mu
    Sigma <- res_Mu_Sigma$Sigma
    if(useY) {
      Delta <- res_Delta$Delta
    }



    # check convergence
    loglik_update <- cal_loglik(Estep_array = Estep_array,
                                Estep_r = Estep_r)
    if(abs(loglik - loglik_update) < tol) {
      flag_converge <- TRUE
      cat("converge!\n")
    } else {
      itr <- itr + 1
      loglik <- loglik_update
      cat(paste0("iteration ", itr, ": log-likelihood = ", loglik_update, "\n"))
    }
  }


  if(!useY) {
    res_Delta <- Mstep_XtoY(Y = Y, r = Estep_r, K = K, N = N,
                            family = family)
    Delta <- res_Delta$Delta
  }

  return(list(res_Beta = res_Beta,
              res_Mu_Sigma = res_Mu_Sigma,
              res_Delta = res_Delta,
              loglik = loglik_update,
              z = Estep_r,
              K = K,
              N = N))
}
