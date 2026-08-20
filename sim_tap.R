# Simulation driver for TAP estimator
# - This script runs Monte Carlo simulations using the TAP core functions
# - Sources: tap_functions.R (core functions)
# - Keeps original simulation logic; comments and filenames cleaned

rm(list = ls())

library(reshape2)
library(ggplot2)
library(gridExtra)
library(rootSolve)
library(corrplot)
library(doParallel)
library(survey)
library(dplyr)
library(latex2exp)
library(plotly)
library(ContourFunctions)
library(RColorBrewer)
library(gginference)

# Source cleaned core functions (renamed from FUNs.R)
source('tap_functions.R')

# ------------------------
# Toy example: TAP MSE surface
LambdaCgammaMSE.plot <- function(lambda, c.gamma,
                                 V.sigma,
                                 eta)
{
  LambdaCgammaMSE(theta = c(lambda, c.gamma),
                   p = 1,
                   V.sigma = V.sigma,
                   eta = eta)
}
LambdaCgammaMSE.plot_v <- Vectorize(LambdaCgammaMSE.plot)

# initialization
eta.vect <- c(0, .5, 1.5)
V.a <- 2; V.dr <- 1; Gamma.A.dr <- .5
V.sigma <- matrix(c(V.a, Gamma.A.dr,
                    Gamma.A.dr, V.dr), ncol = 2)
lambda <- seq(0.5, 40, length.out = 100)
c.gamma <- seq(0.1, 40, length.out = 100)

png('plot_mse_lambda_c_gamma.png', width = 1258, height = 826)
par(mfrow = c(1, 3), cex.main = 2, cex.lab = 1.8,
    mar = c(2, 4, 2, 4))
main.labels <- c('A', 'B', 'C')
adj.vects <- c(.15, .5, .85)
for (i in 1:3) {
  mse <- outer(lambda, c.gamma,
               FUN = LambdaCgammaMSE.plot_v,
               eta = eta.vect[i],
               V.sigma = list(V.sigma))
  library(GA)
  plot3D::persp3D(lambda, c.gamma, mse,
                  zlim = c(0.5, 3.5),
                  theta = 140, phi = 20, expand = 1, d = 30,
                  col.palette = colorRampPalette(rev(brewer.pal(11, 'RdYlBu'))),
                  col = colorRampPalette(rev(brewer.pal(11, 'RdYlBu')))(100),
                  border = 'NA',
                  nlevels = 100,
                  xaxs = "i",
                  xlab = '', ylab = '', zlab = 'MSE',
                  colkey = list(length = 0.3, cex.axis = 1.8)
  )
  mtext(paste0('(', main.labels[i], ')'), side = 3, line = -25, outer = TRUE, adj = adj.vects[i],
        cex = 2)
  plot3D::text3D(15, -7, -1,
                 labels = expression(c[gamma]), add = TRUE,
                 cex = 1.8)
  plot3D::text3D(-5, 12, -1,
                 labels = expression(Lambda), add = TRUE,
                 cex = 1.8)
}
dev.off()

# --------------------------------------
# Figure for adaptive confidence interval
test.T <- function(c.gamma, muTmu) {
  sapply(muTmu, function(x) 1 - pchisq(c.gamma, df = 1, ncp = x / 2))
}
x <- 1:100
h.line <- 30
eps <- 1e-8

png('HT_plot.png', width = 500, height = 400)
par(mar = c(4.1, 2.1, 2.1, 2.1), xpd = TRUE,
    cex = 1.1,
    cex.lab = 1)
curve(test.T(c.gamma = qchisq(1 - .05, df = 1), x), 0, 60, ylim = c(0, 1),
      col = 1,
      lwd = 2,
      xlab = '',
      ylab = '', bty = "n", xaxt = 'n', yaxt = 'n')
axis(1, tck = 0, at = c(-5, 60), labels = c('', expression(paste(mu[2]^T, mu[2]))),
     lty = 1)
axis(2, tck = 0, at = c(-.1, 1), labels = c('', 'P'),
     lty = 1, las = 2)
u <- par("usr")
arrows(u[1], u[3], u[2], u[3], code = 2, xpd = TRUE)
arrows(u[1], u[3], u[1], u[4], code = 2, xpd = TRUE)
rect(0, 0, h.line, 1, col = 'grey90', border = NA)
rect(10, .90, 20, .95, col = 'white', border = NA)
text(55, .92, expression(1 - epsilon), cex = 1.2)
text(55, .06, expression(tilde(alpha)), cex = 1.2)
text(15, .93, paste('B'))
abline(h = test.T(c.gamma = qchisq(1 - .05, df = 1), h.line), lty = 2)
abline(h = test.T(c.gamma = qchisq(1 - eps, df = 1), h.line), lty = 2)
curve(test.T(c.gamma = qchisq(1 - .05, df = 1), x), 0, 60, ylim = c(0, 1),
      col = 1,
      lwd = 2,
      xlab = expression(paste(mu[2]^T, mu[2])),
      ylab = '', bty = "n", add = TRUE)
curve(test.T(c.gamma = qchisq(1 - eps, df = 1), x), 0, 60, add = TRUE,
      col = 1,
      lwd = 2, lty = 2)

legend('bottom',
       legend = c(expression(paste('P(', W[2]^T, W[2], '>', c[gamma], '|', mu[2]^T, mu[2], ')')),
                  expression(paste('P(', T >= v[n], '|', mu[2]^T, mu[2], ')'))),
       col = 1,
       lwd = 2,
       lty = 1:2,
       bty = "n",
       inset = c(0, -0.1), horiz = TRUE,
       text.width = 15)
dev.off()

# ----------------------
# Begin TAP estimation simulations
# Set up parallel computation
sims <- makeCluster(detectCores() - 2)
clusterExport(sims, ls())
registerDoParallel(sims)

set.seed(27695)
# Monte Carlo sample
niter <- 2000
b <- c(0, 200, 2000)
# the core functions are in tap_functions.R
# run simulations (uses GenerateSimuDta and TAP.Est from tap_functions.R)
sim.results.TAP <- sapply(b, function(x) {
  foreach(i = 1:niter, .combine = cbind, .packages = 'dplyr') %dopar% {
    Data.list <- GenerateSimuDta(N = 1e5,
                                 n.A = 1e3, n.B = 1e4,
                                 b = x,
                                 m_model = 'continuous',
                                 seed = x * (i + 1))
    TAP.Est(X.A.all = Data.list$X.A.all,
            weight.A = Data.list$weight.A,
            X.B.all = Data.list$X.B.all,
            Y.name = 'Y',
            m_model = 'continuous',
            B = 1000, K = 100, tau = c(0.125, 0.25, 0.5, 1, 2, 4),
            theta.mean = NULL, theta.usingB.mean = NULL, theta.KH.mean = NULL)
  }
})

stopCluster(sims)
