source("./StatComput/pwpowerlaw.R")
source("./StatComput/Simulation Study/scenarios.R")
set.seed(2022)

# Generator ---------------------------------------------------------------

Data_matrix_compl <- function(R, n, r.theta, p)
{
  replicate(n = R, data.frame(time = rpwpowerlaw(
    n = n, p = p, alpha = r.theta
  ),
  status = 1))
}

Data_matrix_cens <- function(R, n, r.theta, p, PI = 0.3)
{
  replicate(n = R, try(censrand_pwpowerlaw(
    n = n,
    PI = PI,
    p = p,
    alpha = r.theta
  )))
}

# Run (complete data) -----------------------------------------------------

Data <- Data_matrix_compl(R, n = 50, r.theta = Theta_2[1,], p = p_2)
save(Data, file = "./StatComput/Simulation Study/case1_compl_50_2.RData")
Data <- Data_matrix_compl(R, n = 50, r.theta = Theta_2[2,], p = p_2)
save(Data, file = "./StatComput/Simulation Study/case2_compl_50_2.RData")
Data <- Data_matrix_compl(R, n = 50, r.theta = Theta_2[3,], p = p_2)
save(Data, file = "./StatComput/Simulation Study/case3_compl_50_2.RData")

Data <- Data_matrix_compl(R, n = 100, r.theta = Theta_2[1,], p = p_2)
save(Data, file = "./StatComput/Simulation Study/case1_compl_100_2.RData")
Data <- Data_matrix_compl(R, n = 100, r.theta = Theta_2[2,], p = p_2)
save(Data, file = "./StatComput/Simulation Study/case2_compl_100_2.RData")
Data <- Data_matrix_compl(R, n = 100, r.theta = Theta_2[3,], p = p_2)
save(Data, file = "./StatComput/Simulation Study/case3_compl_100_2.RData")

Data <- Data_matrix_compl(R, n = 500, r.theta = Theta_2[1,], p = p_2)
save(Data, file = "./StatComput/Simulation Study/case1_compl_500_2.RData")
Data <- Data_matrix_compl(R, n = 500, r.theta = Theta_2[2,], p = p_2)
save(Data, file = "./StatComput/Simulation Study/case2_compl_500_2.RData")
Data <- Data_matrix_compl(R, n = 500, r.theta = Theta_2[3,], p = p_2)
save(Data, file = "./StatComput/Simulation Study/case3_compl_500_2.RData")

Data <- Data_matrix_compl(R, n = 1000, r.theta = Theta_2[1,], p = p_2)
save(Data, file = "./StatComput/Simulation Study/case1_compl_1000_2.RData")
Data <- Data_matrix_compl(R, n = 1000, r.theta = Theta_2[2,], p = p_2)
save(Data, file = "./StatComput/Simulation Study/case2_compl_1000_2.RData")
Data <- Data_matrix_compl(R, n = 1000, r.theta = Theta_2[3,], p = p_2)
save(Data, file = "./StatComput/Simulation Study/case3_compl_1000_2.RData")

# Run (censored data) -----------------------------------------------------

Data <- Data_matrix_cens(R, n = 50, r.theta = Theta_2[1,], p = p_2)
save(Data, file = "./StatComput/Simulation Study/case1_cens_50_2.RData")
Data <- Data_matrix_cens(R, n = 50, r.theta = Theta_2[2,], p = p_2)
save(Data, file = "./StatComput/Simulation Study/case2_cens_50_2.RData")
Data <- Data_matrix_cens(R, n = 50, r.theta = Theta_2[3,], p = p_2)
save(Data, file = "./StatComput/Simulation Study/case3_cens_50_2.RData")

Data <- Data_matrix_cens(R, n = 100, r.theta = Theta_2[1,], p = p_2)
save(Data, file = "./StatComput/Simulation Study/case1_cens_100_2.RData")
Data <- Data_matrix_cens(R, n = 100, r.theta = Theta_2[2,], p = p_2)
save(Data, file = "./StatComput/Simulation Study/case2_cens_100_2.RData")
Data <- Data_matrix_cens(R, n = 100, r.theta = Theta_2[3,], p = p_2)
save(Data, file = "./StatComput/Simulation Study/case3_cens_100_2.RData")

Data <- Data_matrix_cens(R, n = 500, r.theta = Theta_2[1,], p = p_2)
save(Data, file = "./StatComput/Simulation Study/case1_cens_500_2.RData")
Data <- Data_matrix_cens(R, n = 500, r.theta = Theta_2[2,], p = p_2)
save(Data, file = "./StatComput/Simulation Study/case2_cens_500_2.RData")
Data <- Data_matrix_cens(R, n = 500, r.theta = Theta_2[3,], p = p_2)
save(Data, file = "./StatComput/Simulation Study/case3_cens_500_2.RData")

Data <- Data_matrix_cens(R, n = 1000, r.theta = Theta_2[1,], p = p_2)
save(Data, file = "./StatComput/Simulation Study/case1_cens_1000_2.RData")
Data <- Data_matrix_cens(R, n = 1000, r.theta = Theta_2[2,], p = p_2)
save(Data, file = "./StatComput/Simulation Study/case2_cens_1000_2.RData")
Data <- Data_matrix_cens(R, n = 1000, r.theta = Theta_2[3,], p = p_2)
save(Data, file = "./StatComput/Simulation Study/case3_cens_1000_2.RData")

