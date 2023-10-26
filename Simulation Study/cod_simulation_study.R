source("./StatComput/pwpowerlaw.R")
source("./StatComput/Simulation Study/scenarios.R")
library(data.table)
library(rlist)
library(xtable)
library(numDeriv)
library(maxLik)

mean.sqrt <- function(x) sqrt(mean(x^2))
count_dj <- function (x, z) n_each_interval(x$time[x$status == 1], z)
diff_c <- function(x) c(x[1] - x[2], x[2] - x[3], x[3] - x[4])

Table_compl <- function(n.changepoints) {
  n.grid <- c(50, 100, 500, 1000)
  
  if (n.changepoints == 2){
    p <- p_2
    Theta <- Theta_2
  } else {
    p <- p_3
    Theta <- Theta_3
  }
  
  block_1 <- matrix(c(rep(NA, n.changepoints + 1), paste0("$\\alpha_", 1:(n.changepoints + 1), "$"), rep(NA, (n.changepoints + 1)^3)), ncol = n.changepoints + 1, byrow = T)
  block_1[seq(3, 5* n.changepoints + 1, by = n.changepoints + 1),] <- sprintf(Theta, fmt = '%.1f')
  
  block_2 <- matrix(c(NA, "Estimators", rep(paste0("$\\hat\\alpha_", 1:(n.changepoints + 1), "$"), n.changepoints + 1)))
  
  block_3 <- matrix(NA, ncol = 3 * 4, nrow = nrow(block_1))
  block_3[1, seq(2, 12, by = 3)] <- paste0("$n=", n.grid, "$")
  block_3[2, ] <- rep(c("Bias", "RMSE", "CP"), 4)
  
  for (i in 1:4) { #size
    for (j in 1:3) { #case
      load(paste0("./StatComput/Simulation Study/case", j, "_compl_", n.grid[i], "_", n.changepoints, ".RData"))
      R <- 10000
      Data <- Data[, 1:R]
      aux_bc <- apply(Data, 2, function(x) mle_cens_pwpowerlaw(x, p, T))
      aux_or <- apply(Data, 2, function(x) mle_cens_pwpowerlaw(x, p, F))
      
      #bias
      temp1 <- sprintf(apply(Theta[j, ] - aux_bc, 1, mean), fmt = '%.3f')
      temp2 <- sprintf(apply(Theta[j, ] - aux_or, 1, mean), fmt = '%.3f')
      block_3[(3*j):(3*(j + 1) - 1), 3*(i - 1) + 1] <- paste(temp1, "(", temp2, ")")
      
      #rmse
      temp1 <- sprintf(apply(Theta[j, ] - aux_bc, 1, mean.sqrt), fmt = '%.3f')
      temp2 <- sprintf(apply(Theta[j, ] - aux_or, 1, mean.sqrt), fmt = '%.3f')
      block_3[(3*j):(3*(j + 1) - 1), 3*(i - 1) + 2] <- paste(temp1, "(", temp2, ")")
      
      #cp
      aux1 <- 
        (n.grid[i] * apply(aux_bc, 2, function(x) diff_c(auxiliar(p, x))) / (aux_bc - 1)^2)^(-1/2)
      aux2 <- 
        (n.grid[i] * apply(aux_or, 2, function(x) diff_c(auxiliar(p, x))) / (aux_or - 1)^2)^(-1/2)
      
      lim_inf_1 <- aux_bc - 1.96 * aux1
      lim_sup_1 <- aux_bc + 1.96 * aux1
      lim_inf_2 <- aux_or - 1.96 * aux2
      lim_sup_2 <- aux_or + 1.96 * aux2
      
      xxx <- matrix(Theta[j, ], ncol = R, nrow = 3, byrow = F)
      xxx1 <- data.table::between(xxx, lim_inf_1, lim_sup_1)
      xxx1 <- matrix(xxx1, nrow = 3)
      xxx2 <- data.table::between(xxx, lim_inf_2, lim_sup_2)
      xxx2 <- matrix(xxx2, nrow = 3)
      temp1 <- sprintf(apply(xxx1, 1, mean), fmt = '%.3f')
      temp2 <- sprintf(apply(xxx2, 1, mean), fmt = '%.3f')
      block_3[(3*j):(3*(j + 1) - 1), 3*(i - 1) + 3] <- paste(temp1, "(", temp2, ")")
    }
  }
  out1 <- cbind(block_1, block_2, block_3[,1:6])
  out2 <- cbind(block_1, block_2, block_3[,7:12])
  out <- rbind(out1, out2)
  tbl <- xtable(out, align = rep("c", ncol(out) + 1))
  print.xtable(tbl, sanitize.text.function = function(x) x,
               include.rownames = F, include.colnames = F)
}

Table_cens <- function(n.changepoints) {
  n.grid <- c(50, 100, 500, 1000)
  
  if (n.changepoints == 2){
    p <- p_2
    Theta <- Theta_2
  } else {
    p <- p_3
    Theta <- Theta_3
  }
  
  block_1 <- matrix(c(rep(NA, n.changepoints + 1), paste0("$\\alpha_", 1:(n.changepoints + 1), "$"), rep(NA, (n.changepoints + 1)^3)), ncol = n.changepoints + 1, byrow = T)
  block_1[seq(3, 5* n.changepoints + 1, by = n.changepoints + 1),] <- sprintf(Theta, fmt = '%.1f')
  
  block_2 <- matrix(c(NA, "Estimators", rep(paste0("$\\hat\\alpha_", 1:(n.changepoints + 1), "$"), n.changepoints + 1)))
  
  block_3 <- matrix(NA, ncol = 3 * 4, nrow = nrow(block_1))
  block_3[1, seq(2, 12, by = 3)] <- paste0("$n=", n.grid, "$")
  block_3[2, ] <- rep(c("Bias", "RMSE", "CP"), 4)
  
  for (i in 1:4) { #size
    for (j in 1:3) { #case
      load(paste0("./StatComput/Simulation Study/case", j, "_cens_", n.grid[i], "_", n.changepoints, ".RData"))
      R <- 10000
      Data <- Data[, 1:R]
      aux_or <- apply(Data, 2, function(x) mle_cens_pwpowerlaw(x, p, F))
      
      #bias
      temp2 <- sprintf(apply(Theta[j, ] - aux_or, 1, mean), fmt = '%.3f')
      block_3[(3*j):(3*(j + 1) - 1), 3*(i - 1) + 1] <- temp2
      
      #rmse
      temp2 <- sprintf(apply(Theta[j, ] - aux_or, 1, mean.sqrt), fmt = '%.3f')
      block_3[(3*j):(3*(j + 1) - 1), 3*(i - 1) + 2] <- temp2
      
      #cp
      fun_aux <- function(x, y) {
        log_lik <- function(ttt) loglik_cens_pwpowerlaw(ttt, x, p)
        sqrt(diag(solve(-numDeriv::hessian(log_lik, y))))
      }
      
      aux2 <- (aux_or - 1)/sqrt(apply(Data, 2, function(x) n_each_interval(x$time[x$status == 1], p)))

      lim_inf_2 <- aux_or - 1.96 * aux2
      lim_sup_2 <- aux_or + 1.96 * aux2
      
      xxx <- matrix(Theta[j, ], ncol = R, nrow = 3, byrow = F)
      xxx2 <- data.table::between(xxx, lim_inf_2, lim_sup_2)
      xxx2 <- matrix(xxx2, nrow = 3)
      temp2 <- sprintf(apply(xxx2, 1, mean), fmt = '%.3f')
      block_3[(3*j):(3*(j + 1) - 1), 3*(i - 1) + 3] <- temp2
    }
  }
  out1 <- cbind(block_1, block_2, block_3[,1:6])
  out2 <- cbind(block_1, block_2, block_3[,7:12])
  out <- rbind(out1, out2)
  tbl <- xtable(out, align = rep("c", ncol(out) + 1))
  print.xtable(tbl, sanitize.text.function = function(x) x,
               include.rownames = F, include.colnames = F)
}

Table_compl(n.changepoints = 2)

Table_cens(n.changepoints = 2)

