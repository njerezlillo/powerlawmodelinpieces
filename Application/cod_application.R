source("./StatComput/pwpowerlaw.R")
library(survival)
library(dplyr)
library(maxLik)
library(xtable)
library(ggplot2)
library(poweRlaw)

dat  <- rio::import("./StatComput/Application/dataset.xlsx") %>% 
  mutate(DurationOfPower = StartGovernment1 - EndGovernment1 + 1)
x <- na.omit(dat$DurationOfPower)
n <- length(x)

# Change points estimation ------------------------------------------------

target <- function(y) profile_loglik_pwpowerlaw(p = c(1, y), x)

# k = 1

A_1 <- matrix(1, ncol = 1, byrow = T)
d_1 <- -1
fit_p_1 <- maxSANN(target, start = 10,
                   constraints = list(ineqA = A_1, ineqB = d_1))
fit_p_1 <- c(1, fit_p_1$estimate)
fit_alpha_1 <- mle_pwpowerlaw(x, fit_p_1)

# k = 2

A_2 <- matrix(c(1, 0, -1, 1), ncol = 2, byrow = T)
d_2 <- c(-1, 0)
fit_p_2 <- maxSANN(target, start = c(5, 15),
                   constraints = list(ineqA = A_2, ineqB = d_2))
fit_p_2 <- c(1, fit_p_2$estimate)
fit_alpha_2 <- mle_pwpowerlaw(x, fit_p_2)

# k = 3

A_3 <- matrix(c(1, 0, 0, -1, 1, 0, 0, -1, 1), ncol = 3, byrow = T)
d_3 <- c(-1, 0, 0)
fit_p_3 <- maxSANN(target, start = c(6, 10, 30),
                   constraints = list(ineqA = A_3, ineqB = d_3))
fit_p_3 <- c(1, fit_p_3$estimate)
fit_alpha_3 <- mle_pwpowerlaw(x, fit_p_3)

# k = 4

A_4 <- matrix(c(1, 0, 0, 0, -1, 1, 0, 0, 0, -1, 1, 0, 0, 0, -1, 1), ncol = 4, byrow = T)
d_4 <- c(-1, 0, 0, 0)
fit_p_4 <- maxSANN(target, start = c(3, 10, 25, 30),
                   constraints = list(ineqA = A_4, ineqB = d_4))
fit_p_4 <- c(1, fit_p_4$estimate)
fit_alpha_4 <- mle_pwpowerlaw(x, fit_p_4)

# AIC and BIC -------------------------------------------------------------

AIC_1 <- 2 * 3 - 2 * loglik_pwpowerlaw(fit_alpha_1, x, fit_p_1)
AIC_2 <- 2 * 5 - 2 * loglik_pwpowerlaw(fit_alpha_2, x, fit_p_2)
AIC_3 <- 2 * 7 - 2 * loglik_pwpowerlaw(fit_alpha_3, x, fit_p_3)

BIC_1 <- log(n) * 3 - 2 * loglik_pwpowerlaw(fit_alpha_1, x, fit_p_1)
BIC_2 <- log(n) * 5 - 2 * loglik_pwpowerlaw(fit_alpha_2, x, fit_p_2)
BIC_3 <- log(n) * 7 - 2 * loglik_pwpowerlaw(fit_alpha_3, x, fit_p_3)

gof <- data.frame(
  "k" = 1:3,
  AIC = c(AIC_1, AIC_2, AIC_3),
  BIC = c(BIC_1, BIC_2, BIC_3)
)

print(xtable(gof), include.rownames = F)

# Coefficients ------------------------------------------------------------

nj <- n_each_interval(x, fit_p_2)
se <- (fit_alpha_2 - 1)/sqrt(sum(nj)*diff_c(auxiliar(fit_p_2, fit_alpha_2)))
ic_l <- round(fit_alpha_2 - 1.96 * se, 2)
ic_u <- round(fit_alpha_2 + 1.96 * se, 2)

fit <- data.frame(j = 1:3,
           "estimation" = fit_alpha_2,
           "CI" = paste0("(", ic_l, ", ", ic_u, ")"))

print(xtable(fit), include.rownames = F)

# Sequential selection ----------------------------------------------------

w <- (fit_alpha_2[1] - fit_alpha_2[2])^2/sum(se[1:2]^2)
w
pchisq(w, 1, lower.tail = F)

significance <- 0.01
var <- (fit_alpha_1 - 1)^2/(sum(nj)*diff_c(auxiliar(fit_p_1, fit_alpha_1))[1:2])
w <- (fit_alpha_1[1] - fit_alpha_1[2])^2/sum(var[1:2])
pchisq(w, 1, lower.tail = F) < significance

var <- (fit_alpha_2 - 1)^2/(sum(nj)*diff_c(auxiliar(fit_p_2, fit_alpha_2)))
w1 <- (fit_alpha_2[1] - fit_alpha_2[2])^2/sum(var[1:2])
w2 <- (fit_alpha_2[2] - fit_alpha_2[3])^2/sum(var[2:3])
pchisq(w1, 1, lower.tail = F) < significance/2
pchisq(w2, 1, lower.tail = F) < significance/2

var <- (fit_alpha_3 - 1)^2/(sum(nj)*c(diff_c(auxiliar(fit_p_3, fit_alpha_3)), 0.1463))
w1 <- (fit_alpha_3[1] - fit_alpha_3[2])^2/sum(var[1:2])
w2 <- (fit_alpha_3[2] - fit_alpha_3[3])^2/sum(var[2:3])
w3 <- (fit_alpha_3[3] - fit_alpha_3[4])^2/sum(var[3:4])
pchisq(w1, 1, lower.tail = F) < significance/4
pchisq(w2, 1, lower.tail = F) < significance/4
pchisq(w3, 1, lower.tail = F) < significance/4

# Histogram ---------------------------------------------------------------

ggplot(dat, aes(x = DurationOfPower)) + 
  geom_histogram(color = "black", fill = "gray90") +
  labs(x = "Duration of Government (Years)", y = "Frequency") +
  theme_classic()
#ggsave("./Application/Fig3.pdf", height = 6, width = 10)
#ggsave("./Application/Fig3.eps", height = 6, width = 10)

# Fig 4 -------------------------------------------------------------------

#cleopatra 22, tutancamon 10, Ramessés II 68

aux <- displ$new(x)
xmin <- estimate_xmin(aux)$xmin
mle <- 1 + (sum(x > xmin) - 1) * (sum(log(x[x > xmin]/xmin)))^-1
C <- mean(x > xmin)
S <- function(z) poweRlaw::ppldis(z, xmin, mle, lower.tail = F) * C
S <- Vectorize(S)
S_fit_2 <- function(t) spwpowerlaw(t, fit_p_2, fit_alpha_2)
S_fit_2 <- Vectorize(S_fit_2)

KS_fit <- survfit(Surv(x, rep(1, n)) ~ 1) 
dat %>% filter(Ruler %in% c("Cleopatra VII", "Tutancãmon", "Ramessés II"))

#postscript("./Application/Fig4.eps", height = 6, width = 10, horizontal = F,onefile = FALSE, paper = "special")
#pdf("./Application/Fig4.pdf", height = 6, width = 10)
par(mar = c(5, 4.5, 2, 2))
plot(KS_fit, xlab = "Time (Years)", ylab = "Survival",
   conf.int = F, bty = "n", xlim = c(0, 100))
curve(S, xmin, 100, add = T, col = "red", lwd = 2, lty = 3)
segments(xmin, 0, xmin, S(xmin), col = "red", lty = 2)
curve(S_fit_2, 0.2, 100, add = T, col = "darkblue", lwd = 2, lty = 2)
abline(v = fit_p_2, lty = 2, col = "gray")
legend(60, 0.89, c("power-law", "piecewise (k = 2)"), bty = "n",
     lwd = c(2, 2), col = c("red", "darkblue"), lty = c(3, 2))
points(c(10, 22, 68), c(0.43, 0.21, 0), pch = 16)
text(18, 0.5, "Tutankhamun")
arrows(15, 0.47, 11, 0.44, length = 0.08)
text(68, 0.08, "Ramessés II")
arrows(68, 0.05, 68, 0.02, length = 0.08)
text(30, 0.25, "Cleopatra VII")
arrows(28, 0.22, 25, 0.215, length = 0.08)
#dev.off()

