source("./StatComput/pwpowerlaw.R")

tovec <- function (z) {
  out <- paste(sprintf("%.01f", z), collapse = ', ')
  out <- c("(", out, ")")
  print(paste(out, collapse = ''))
}

# p = 3 -------------------------------------------------------------------

p <- c(1, 3, 5)
alpha1 <- c(1.5, 4, 6); alpha2 <- c(3, 2, 4); alpha3 <- c(2, 6, 2)

s1 <- function(x) spwpowerlaw(x, p, alpha1)
s2 <- function(x) spwpowerlaw(x, p, alpha2)
s3 <- function(x) spwpowerlaw(x, p, alpha3)
s1 <- Vectorize(s1, "x")
s2 <- Vectorize(s2, "x")
s3 <- Vectorize(s3, "x")

d1 <- function(x) dpwpowerlaw(x, p, alpha1)
d2 <- function(x) dpwpowerlaw(x, p, alpha2)
d3 <- function(x) dpwpowerlaw(x, p, alpha3)
d1 <- Vectorize(d1, "x")
d2 <- Vectorize(d2, "x")
d3 <- Vectorize(d3, "x")

grid <- seq(1.01, 10, by = 0.01)

#postscript("Fig1.eps", width = 12, height = 6,horizontal = F,onefile = FALSE, paper = "special")
#pdf("Fig1.pdf", width = 12, height = 6)
par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))

plot(grid, s1(grid), xlim = c(0, 8), ylim = c(0, 1), main = "(a)",
     type = "l", col = "red", ylab = "S(t)", xlab = "t", bty = "n")
lines(grid, s2(grid), col = "darkblue", lty = 2)
lines(grid, s3(grid), col = "darkgreen", lty = 3)
segments(p, rep(-0.2, 3), p, rep(0, 3), col = "gray")
legend("topright", c(tovec(alpha1), tovec(alpha2), tovec(alpha3)),
       col = c("red","darkblue", "darkgreen"), lty = c(1, 2, 3), bty = "n")

plot(grid, d1(grid), xlim = c(0, 8), main = "(b)", ylim = c(0, 0.6),
     type = "l", col = "red", ylab = "f(t)", xlab = "t", bty = "n")
lines(grid, d2(grid), col = "darkblue", lty = 2)
lines(grid, d3(grid), col = "darkgreen", lty = 3)
segments(p, rep(-0.2, 3), p, rep(0, 3), col = "gray")
legend("topright", c(tovec(alpha1), tovec(alpha2), tovec(alpha3)),
       col = c("red","darkblue", "darkgreen"), lty = c(1, 2, 3), bty = "n")

par(mfrow = c(1, 1))
#dev.off()

# p = 4 -------------------------------------------------------------------

p <- c(1, 3, 5, 8)
alpha1 <- c(1.2, 3, 6, 3); alpha2 <- c(2, 2, 4, 9); alpha3 <- c(1.5, 2.5, 9, 2)

s1 <- function(x) spwpowerlaw(x, p, alpha1)
s2 <- function(x) spwpowerlaw(x, p, alpha2)
s3 <- function(x) spwpowerlaw(x, p, alpha3)
s1 <- Vectorize(s1, "x")
s2 <- Vectorize(s2, "x")
s3 <- Vectorize(s3, "x")

d1 <- function(x) dpwpowerlaw(x, p, alpha1)
d2 <- function(x) dpwpowerlaw(x, p, alpha2)
d3 <- function(x) dpwpowerlaw(x, p, alpha3)
d1 <- Vectorize(d1, "x")
d2 <- Vectorize(d2, "x")
d3 <- Vectorize(d3, "x")

grid <- seq(1.01, 15, by = 0.01)

#postscript("Fig2.eps", width = 12, height = 6, horizontal = F,onefile = FALSE, paper = "special")
#pdf("Fig2.pdf", width = 12, height = 6)
par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))

plot(grid, s1(grid), xlim = c(0, 15), ylim = c(0, 1), main = "(a)",
     type = "l", col = "red", ylab = "S(t)", xlab = "t", bty = "n")
lines(grid, s2(grid), col = "darkblue", lty = 2)
lines(grid, s3(grid), col = "darkgreen", lty = 3)
segments(p, rep(-0.2, 3), p, rep(0, 3), col = "gray")
legend("topright", c(tovec(alpha1), tovec(alpha2), tovec(alpha3)),
       col = c("red","darkblue", "darkgreen"), lty = c(1, 2, 3), bty = "n")

plot(grid, d1(grid), xlim = c(0, 15), main = "(b)", ylim = c(0, 0.6),
     type = "l", col = "red", ylab = "f(t)", xlab = "t", bty = "n")
lines(grid, d2(grid), col = "darkblue", lty = 2)
lines(grid, d3(grid), col = "darkgreen", lty = 3)
segments(p, rep(-0.2, 3), p, rep(0, 3), col = "gray")
legend("topright", c(tovec(alpha1), tovec(alpha2), tovec(alpha3)),
       col = c("red","darkblue", "darkgreen"), lty = c(1, 2, 3), bty = "n")

par(mfrow = c(1, 1))
#dev.off()
