### Core TAP functions (cleaned and comment-fixed)
# NOTE: Function names kept the same as in the original code for compatibility
# of existing calls (e.g. GenerateSimuDta). Only comments and minor formatting
# were fixed.

# --- simulated data --------------------------------------------------------
GenerateSimuDta <- function(N=1e5,
                            n.A=100, n.B=5e3,
                            b=.1,
                            m_model='continuous',
                            seed)
{
  set.seed(27695)
  # Build a population covariate matrix X (intercept + two covariates)
  X <- matrix(c(rep(1,N),
                rnorm(N),
                rnorm(N, mean=1)),
              ncol = 3)
  # Unmeasured confounder
  u <- rnorm(N, mean=0)
  X.all <- cbind(X,u)

  link.all <- X%*%c(1,1,1)+.5*u*b*n.B^(-1/2)+
    rnorm(N)

  if(m_model=='continuous')
  {
    Y.all <- link.all
  }
  if(m_model=='binary')
  {
    Y.all <- exp(link.all)/(1+exp(link.all))
  }

  # Sample A (non-simple sampling)
  select.prob.A <- exp(X%*%c(.5, .2,.1))/
    (1+exp(X%*%c(.5, .2,.1)))
  select.prob.A <- select.prob.A/sum(select.prob.A)

  # Use provided seed for sampling reproducibility
  set.seed(seed)
  A.chosen <- sample(1:N,
                     size=n.A,
                     prob = select.prob.A)
  weight.A <- 1/(n.A*select.prob.A[A.chosen])
  X.Ac <- X[A.chosen,]
  Y.Ac <- Y.all[A.chosen]

  # Sample B (non-probability sample with additional confounder u)
  select.prob.B <- exp(cbind(X, u)%*%c(.5,.1,.2, .1))/
    (1+exp(cbind(X, u)%*%c(.5,.1,.2, .1)))
  B.chosen <- sample(1:N,
                     size = n.B,
                     prob = select.prob.B)
  X.Bc <- X[B.chosen,]
  Y.Bc <- Y.all[B.chosen]

  return(list(
    X.A.all=data.frame(X.A=X.Ac,
                       Y=Y.Ac),
    X.B.all = data.frame(X.B=X.Bc,
                         Y=Y.Bc),
    weight.A=weight.A
  ))
}

# --- Truncated normal first moment ----------------------------------------
NormalTruncatedFirstMom <- function(mu,
                                    p=1,
                                    a=0,
                                    b=Inf)
{
  if(a==b) return(0)
  else {
    ncp.mu <- mu^2
    dem <- pchisq(b, df=p, ncp = ncp.mu)-
      pchisq(a, df=p, ncp = ncp.mu)
    num <- mu*(
      pchisq(b, df=p+2, ncp = ncp.mu)-
        pchisq(a, df=p+2, ncp = ncp.mu)
    )
    return(num/dem)
  }
}

# --- Truncated normal second moment ---------------------------------------
NormalTruncatedSecondMom <- function(mu,
                                     p=1,
                                     a=0,
                                     b=Inf)
{
  if(a==b) return(0)
  else {
    ncp.mu <- mu^2
    dem <- pchisq(b, df=p, ncp = ncp.mu)-
      pchisq(a, df=p, ncp = ncp.mu)
    num1 <- pchisq(b, df=p+2, ncp = ncp.mu)-
      pchisq(a, df=p+2, ncp = ncp.mu)
    num2 <- ncp.mu*(
      pchisq(b, df=p+4, ncp = ncp.mu)-
        pchisq(a, df=p+4, ncp = ncp.mu))
    (num1+num2)/dem
  }
}

# --- Asymptotic bias for TAP estimator ------------------------------------
LambdaCgammaBias <- function(theta,
                             p=1,
                             V.sigma,
                             eta)
{
  lambda <- theta[1]
  c.gamma <- theta[2]
  V.a <- V.sigma[1,1]
  Gamma.A.dr <- V.sigma[1,2]
  V.dr <- V.sigma[2,2]

  mu.1 <- eta*(Gamma.A.dr-V.a)/
    ((V.a+V.dr-2*Gamma.A.dr)*(V.a*V.dr-Gamma.A.dr^2))^(1/2)
  mu.2 <- -eta/(V.a+V.dr-2*Gamma.A.dr)^(1/2)

  p1 <- -(V.a*V.dr-Gamma.A.dr^2)^(1/2)/
    (V.a+V.dr-2*Gamma.A.dr)^(1/2)*
    NormalTruncatedFirstMom(mu.1)
  p2 <- (lambda*(Gamma.A.dr-V.dr)-(Gamma.A.dr-V.a))/
    (1+lambda)/(V.a+V.dr-2*Gamma.A.dr)^(1/2)*
    NormalTruncatedSecondMom(mu.2, b=c.gamma)

  p3 <- p1
  p4 <- -(Gamma.A.dr-V.a)/
    (V.a+V.dr-2*Gamma.A.dr)^(1/2)*
    NormalTruncatedFirstMom(mu.2, a=c.gamma)

  xi <- pchisq(c.gamma, df=p, ncp = mu.2^2)

  (p1+p2)*xi+(p3+p4)*(1-xi)
}

# --- Asymptotic variance for TAP estimator -------------------------------
LambdaCgammaVar <- function(theta,
                            p=1,
                            V.sigma,
                            eta)
{
  lambda <- theta[1]
  c.gamma <- theta[2]
  V.a <- V.sigma[1,1]
  Gamma.A.dr <- V.sigma[1,2]
  V.dr <- V.sigma[2,2]

  mu.1 <- eta*(Gamma.A.dr-V.a)/
    ((V.a+V.dr-2*Gamma.A.dr)*(V.a*V.dr-Gamma.A.dr^2))^(1/2)
  mu.2 <- -eta/(V.a+V.dr-2*Gamma.A.dr)^(1/2)

  p1 <- (V.a*V.dr-Gamma.A.dr^2)/
    (V.a+V.dr-2*Gamma.A.dr)

  p2 <- ((lambda*(Gamma.A.dr-V.dr)-(Gamma.A.dr-V.a))/
           (1+lambda)/(V.a+V.dr-2*Gamma.A.dr)^(1/2))^2*
    (NormalTruncatedSecondMom(mu.2, b=c.gamma)-
       NormalTruncatedFirstMom(mu.2, b=c.gamma)^2)

  p3 <- p1

  p4 <- (-(Gamma.A.dr-V.a)/
           (V.a+V.dr-2*Gamma.A.dr)^(1/2))^2*
    (NormalTruncatedSecondMom(mu.2, a=c.gamma)-
       NormalTruncatedFirstMom(mu.2, a=c.gamma)^2)

  xi <- pchisq(c.gamma, df=p, ncp = mu.2^2)

  (p1+p2)*xi+(p3+p4)*(1-xi)
}

# --- Asymptotic MSE (analytic form) --------------------------------------
LambdaCgammaMSE <- function(theta,
                            p=1,
                            V.sigma,
                            eta, rho=1)
{
  lambda <- theta[1]
  c.gamma <- theta[2]
  V.a <- V.sigma[1,1]
  Gamma.A.dr <- V.sigma[1,2]
  V.dr <- V.sigma[2,2]

  sigma.T <- rho*(V.a+V.dr-2*Gamma.A.dr)
  sigma.S <- rho*(V.a+V.dr-2*Gamma.A.dr)*
    (V.a*V.dr-Gamma.A.dr^2)

  mu.1 <- sigma.S^(-1/2)*(Gamma.A.dr-V.a)*eta
  mu.2 <- -sigma.T^(-1/2)*eta

  p1.bias <- -(V.a*V.dr-Gamma.A.dr^2)^(1/2)/
    (V.a+V.dr-2*Gamma.A.dr)^(1/2)*
    NormalTruncatedFirstMom(mu.1)

  p2.bias <- (abs((Gamma.A.dr-V.a))- 
                lambda*(abs(Gamma.A.dr-V.dr)))/
    (1+lambda)/(V.a+V.dr-2*Gamma.A.dr)^(1/2)*
    NormalTruncatedFirstMom(mu.2, b=c.gamma)

  p1.var <- (V.a*V.dr-Gamma.A.dr^2)/
    (V.a+V.dr-2*Gamma.A.dr)

  p2.var <- ((abs((Gamma.A.dr-V.a))- 
                lambda*(abs(Gamma.A.dr-V.dr)))/
               (1+lambda)/(V.a+V.dr-2*Gamma.A.dr)^(1/2))^2*
    (NormalTruncatedSecondMom(mu.2, b=c.gamma)-
       NormalTruncatedFirstMom(mu.2, b=c.gamma)^2)

  p3.bias <- p1.bias
  p4.bias <- abs(Gamma.A.dr-V.a)/
    (V.a+V.dr-2*Gamma.A.dr)^(1/2)*
    NormalTruncatedFirstMom(mu.2, a=c.gamma)

  p3.var <- p1.var
  p4.var <- (abs(Gamma.A.dr-V.a)/
               (V.a+V.dr-2*Gamma.A.dr)^(1/2))^2*
    (NormalTruncatedSecondMom(mu.2, a=c.gamma)-
       NormalTruncatedFirstMom(mu.2, a=c.gamma)^2)

  xi <- pchisq(c.gamma, df=p, ncp = mu.2^2)

  ((p1.bias+p2.bias)^2+p1.var+p2.var)*xi+
    ((p3.bias+p4.bias)^2+p3.var+p4.var)*(1-xi)
}

# --- Monte Carlo approximation of MSE -------------------------------------
LambdaCgammaMSE.MC <- function(lambda.star,
                               c.gamma,
                               W1.series,
                               W2.series,
                               V.sigma,
                               eta)
{
  V.a <- V.sigma[1,1]
  Gamma.A.dr <- V.sigma[1,2]
  V.dr <- V.sigma[2,2]
  V.eff <- (V.a+V.dr-2*Gamma.A.dr)^(-1)*
    (V.a*V.dr-Gamma.A.dr^2)
  V.a_eff <- V.a-V.eff
  V.dr_eff <- V.dr-V.eff

  mu.tap.series <- -V.eff^(1/2)*W1.series+
    (1+lambda.star)^(-1)*(V.a_eff^(1/2)-lambda.star*V.dr_eff^(1/2))*W2.series+
    (W2.series^2>=c.gamma)*
    (1+lambda.star)^(-1)*lambda.star*(V.dr_eff^(1/2)+V.a_eff^(1/2))*W2.series
  return(sd(mu.tap.series))
}

OptimLambdaCgamma.MC <- function(V.sigma, eta,
                                 W1.series,
                                 W2.series,
                                 lamba.range,
                                 c.gamma.range)
{
  lambda.gamma.grid <- expand.grid(lambda.star=lamba.range,
                                   c.gamma=c.gamma.range)
  out <- do.call(mapply, 
                 c(list(FUN=LambdaCgammaMSE.MC,
                        MoreArgs = list(V.sigma=V.sigma, eta=eta,
                                        W1.series=W1.series,
                                        W2.series=W2.series)),
                   lambda.gamma.grid))
  theta.opt <- lambda.gamma.grid[which.min(out),]%>%as.numeric()
  return(theta.opt)
}

# --- Main TAP estimator ---------------------------------------------------
TAP.Est <- function(X.A.all, weight.A,
                    X.B.all,
                    Y.name='age',
                    m_model='continuous',
                    B=100,
                    K=1, # number of double bootstrap datasets
                    tau,  # candidates for double bootstrap
                    theta.mean, theta.usingB.mean, theta.KH.mean
)
{
  library(BB)
  n.A <- nrow(X.A.all)
  n.B <- nrow(X.B.all)
  n <- n.A+n.B
  rho <- n.B/n

  X.A.all.r <- (cbind(X.A.all[, !colnames(X.A.all)%in%Y.name])%>%
                  as.matrix())
  Y.A <- X.A.all[, colnames(X.A.all)%in%Y.name]
  X.B.all.r <- (cbind(X.B.all[, !colnames(X.B.all)%in%Y.name])%>%
                  as.matrix())
  Y.B <- X.B.all[, colnames(X.B.all)%in%Y.name]
  p <- ncol(X.A.all)
  p.r <- ncol(X.A.all.r)

  # Solve propensity score for non-probability sample (NPR)
  Prop_NPR <- function(theta)
  {
    sum(as.matrix(X.B.all.r)%*%theta)-
      sum(weight.A*
            log(1+exp(as.matrix(X.A.all.r)%*%
                        theta)))
  }

  Prop_NPR_gr <- function(theta)
  {
    pi.A.temp <- exp(as.matrix(X.A.all.r)%*%theta)/
      (1+exp(as.matrix(X.A.all.r)%*%theta))
    gr <- apply(X.B.all.r, 2, sum)-
      apply(as.numeric(weight.A*pi.A.temp)*
              as.matrix(X.A.all.r),2,sum)
    gr
  }

  theta.B <- dfsane(par = rep(-.5,p.r), # initial points
                    fn=Prop_NPR_gr,
                    control = list(trace=T,
                                   NM=T,
                                   BFGS=T,
                                   tol=1.e-2))$par

  pi.B.est <- exp(as.matrix(X.B.all.r)%*%as.matrix(theta.B))/
    (1+exp(as.matrix(X.B.all.r)%*%as.matrix(theta.B)))

  N.A <- sum(weight.A)
  pi.A <- 1/weight.A

  X.n <- rbind(X.A.all.r,
               X.B.all.r)%>%as.matrix()
  Y.n <- c(Y.A,Y.B)%>%as.matrix()

  if(m_model=='continuous')
  {
    beta.v <- solve(t(X.n)%*%X.n)%*%
      t(X.n)%*%Y.n
    beta.usingB <- solve(t(X.B.all.r)%*%
                           X.B.all.r)%*%
      t(X.B.all.r)%*%Y.B

    KH.estimator <- function(theta)
    {
      r <- rep(NA, 2*p.r)
      pi.B.temp <- exp(as.matrix(X.B.all.r)%*%theta[1:p.r])/
        (1+exp(as.matrix(X.B.all.r)%*%theta[1:p.r]))

      r[1:p.r] <- apply(X.B.all.r*(1/1/pi.B.temp)%>%as.numeric(),2,sum)-
        apply(X.A.all.r*(1/1/pi.A)%>%as.numeric(),2,sum)

      r[(p.r+1):(2*p.r)] <- apply(X.B.all.r*
                                    ((1/pi.B.temp-1)*
                                       (Y.B-X.B.all.r%*%
                                          theta[(p.r+1):(2*p.r)]))%>%
                                    as.numeric(),2,sum)
      r
    }

    theta.KH <- dfsane(par = c(theta.B, beta.usingB),
                       fn=KH.estimator,
                       control = list(trace=T,
                                      NM=T,
                                      BFGS=T,
                                      tol=1.e-2))$par
    theta.B.KH <- theta.KH[1:p.r]
    beta.KH <- theta.KH[(p.r+1):(2*p.r)]

    Y.A.est <- X.A.all.r%*%beta.v
    Y.B.est <- X.B.all.r%*%beta.v
    Y.A.est.usingB <- X.A.all.r%*%beta.usingB
    Y.B.est.usingB <- X.B.all.r%*%beta.usingB
    Y.A.est.KH <- X.A.all.r%*%beta.KH
    Y.B.est.KH <- X.B.all.r%*%beta.KH

    pi.B.est.KH <- exp(as.matrix(X.B.all.r)%*%as.matrix(theta.B.KH))/
      (1+exp(as.matrix(X.B.all.r)%*%as.matrix(theta.B.KH)))
  }

  if(m_model=='binary')
  {
    beta.model.v <- glm(Y.n~.+0,
                        data=data.frame(Y.n,
                                        X.n),
                        family = binomial(link = 'logit'))
    beta.model.usingB <- glm(Y.B~.+0,
                             data=data.frame(Y.B,
                                             X.B.all.r),
                             family = binomial(link = 'logit'))

    KH.estimator.logit <- function(theta)
    {
      r <- rep(NA, 2*p.r)
      pi.B.temp <- exp(as.matrix(X.B.all.r)%*%theta[1:p.r])/
        (1+exp(as.matrix(X.B.all.r)%*%theta[1:p.r]))
      m.b <-exp(X.B.all.r%*%theta[(p.r+1):(2*p.r)])/
        (1+exp(X.B.all.r%*%theta[(p.r+1):(2*p.r)]))
      m.A.gr <-  X.A.all.r*as.numeric(exp(X.A.all.r%*%theta[(p.r+1):(2*p.r)])/
                                        (1+exp(X.A.all.r%*%theta[(p.r+1):(2*p.r)]))^2)
      m.B.gr <-  X.B.all.r*as.numeric(exp(X.B.all.r%*%theta[(p.r+1):(2*p.r)])/
                                        (1+exp(X.B.all.r%*%theta[(p.r+1):(2*p.r)]))^2)

      r[1:p.r] <- apply(m.B.gr*(1/pi.B.temp)%>%as.numeric(),2,sum)-
        apply(m.A.gr*(1/pi.A)%>%as.numeric(),2,sum)
      r[(p.r+1):(2*p.r)] <- apply(m.B.gr*((1/pi.B.temp-1)*(Y.B-m.b))%>%
                                    as.numeric(),2,sum)
      r
    }

    theta.KH <- dfsane(par = c(theta.B, coef(beta.model.v)%>%as.numeric()),
                       fn=KH.estimator.logit,
                       control = list(trace=T,
                                      NM=T,
                                      BFGS=T,
                                      tol=1.e-2))$par
    theta.B.KH <- theta.KH[1:p.r]
    beta.KH <- theta.KH[(p.r+1):(2*p.r)]

    Y.A.est <- predict(beta.model.v, type='response',
                       newdata=data.frame(X.A.all.r))%>%
      as.numeric()
    Y.B.est <- predict(beta.model.v, type='response',
                       newdata=data.frame(X.B.all.r))%>%
      as.numeric()
    Y.A.est.usingB <- predict(beta.model.usingB, type='response',
                              newdata=data.frame(X.A.all.r))%>%
      as.numeric()
    Y.B.est.usingB <- predict(beta.model.usingB, type='response',
                              newdata=data.frame(X.B.all.r))%>%
      as.numeric()
    Y.A.est.KH <- exp(X.A.all.r%*%beta.KH)/(1+exp(X.A.all.r%*%beta.KH))
    Y.B.est.KH <- exp(X.B.all.r%*%beta.KH)/(1+exp(X.B.all.r%*%beta.KH))

    pi.B.est.KH <- exp(as.matrix(X.B.all.r)%*%as.matrix(theta.B.KH))/
      (1+exp(as.matrix(X.B.all.r)%*%as.matrix(theta.B.KH)))
  }

  # --- Point estimates ----------------------------------------------------
  mu.A <-  sum(Y.A/pi.A)/N.A
  mu.dr <- sum((Y.B-Y.B.est)/pi.B.est)/sum(1/pi.B.est)+
    sum(Y.A.est/pi.A)/N.A
  mu.dr.usingB <- sum((Y.B-Y.B.est.usingB)/pi.B.est)/sum(1/pi.B.est)+
    sum(Y.A.est.usingB/pi.A)/N.A
  mu.dr.KH <- sum((Y.B-Y.B.est.KH)/pi.B.est.KH)/sum(1/pi.B.est.KH)+
    sum(Y.A.est.KH/pi.A)/N.A

  # --- Bootstrap series --------------------------------------------------
  mu.A.boot <- list()
  mu.dr.boot <- list()
  mu.dr.usingB.boot <- list()
  mu.dr.KH.boot <- list()

  for (k in 1:K) {
    mu.A.series <- numeric(B)
    mu.dr.series <- numeric(B)
    mu.dr.usingB.series <- numeric(B)
    mu.dr.KH.series <- numeric(B)
    for (i in 1:B) 
    {
      A.chosen.boot <- sample(1:n.A, replace = T, size = n.A)
      B.chosen.boot <- sample(1:n.B, replace = T, size = n.B)
      X.A.boot <- X.A.all.r[A.chosen.boot,]
      X.B.boot <- X.B.all.r[B.chosen.boot,]
      Y.A.boot <- Y.A[A.chosen.boot]
      Y.B.boot <- Y.B[B.chosen.boot]

      X.n.boot <- rbind(X.A.boot, X.B.boot)
      Y.n.boot <- c(Y.A.boot, Y.B.boot)

      mu.A.temp <- sum(Y.A.boot/pi.A[A.chosen.boot])/sum(weight.A[A.chosen.boot])

      mu.dr.temp <- sum((Y.B.boot-Y.B.est[B.chosen.boot])/
                          pi.B.est[B.chosen.boot,])/sum(1/pi.B.est[B.chosen.boot,])+ 
        sum(Y.A.est[A.chosen.boot]/pi.A[A.chosen.boot])/sum(weight.A[A.chosen.boot])

      mu.dr.usingB.temp <- sum((Y.B.boot-Y.B.est.usingB[B.chosen.boot])/
                                 pi.B.est[B.chosen.boot,])/sum(1/pi.B.est[B.chosen.boot,])+ 
        sum(Y.A.est.usingB[A.chosen.boot]/pi.A[A.chosen.boot])/sum(weight.A[A.chosen.boot])

      mu.dr.KH.temp <- sum((Y.B.boot-Y.B.est.KH[B.chosen.boot])/
                             pi.B.est.KH[B.chosen.boot,])/sum(1/pi.B.est.KH[B.chosen.boot,])+ 
        sum(Y.A.est.KH[A.chosen.boot]/pi.A[A.chosen.boot])/sum(weight.A[A.chosen.boot])

      mu.dr.series[i] <- mu.dr.temp
      mu.A.series[i] <- mu.A.temp
      mu.dr.usingB.series[i] <- mu.dr.usingB.temp
      mu.dr.KH.series[i] <- mu.dr.KH.temp
    }
    mu.A.boot[[k]] <- mu.A.series
    mu.dr.boot[[k]] <- mu.dr.series
    mu.dr.usingB.boot[[k]] <- mu.dr.usingB.series
    mu.dr.KH.boot[[k]] <- mu.dr.KH.series
  }

  TAP_variance <- function(V.mat,
                           MU.A, MU.B,
                           MU.A.series,
                           MU.B.series,
                           rho)
  {
    V.a <- V.mat[1,1]
    Gamma.A.dr <- V.mat[1,2]
    V.dr <- V.mat[2,2]
    sigma.T <- rho*(V.a+V.dr-2*Gamma.A.dr)
    sigma.S <- rho*(V.a+V.dr-2*Gamma.A.dr)*
      (V.a*V.dr-Gamma.A.dr^2)

    V.eff <- (V.a+V.dr-2*Gamma.A.dr)^(-1)*
      (V.a*V.dr-Gamma.A.dr^2)
    V.a_eff <- V.a-V.eff
    V.dr_eff <- V.dr-V.eff

    W2.hat <- rho^(1/2)*sigma.T^(-1/2)*n^(1/2)*(mu.A-mu.dr)

    W1.series <- rho^(1/2)*sigma.S^(-1/2)*n^(1/2)*
      ((Gamma.A.dr-V.dr)*(MU.A.series-mean(MU.A.series))+ 
         (Gamma.A.dr-V.a)*(MU.B.series-mean(MU.A.series)))
    W2.series <- rho^(1/2)*sigma.T^(-1/2)*n^(1/2)*(MU.A.series-MU.B.series)

    eta.est <- -sigma.T^(1/2)*mean(W2.series)

    mu.1 <- sigma.S^(-1/2)*(Gamma.A.dr-V.a)*eta.est
    mu.2 <- -sigma.T^(-1/2)*eta.est

    T.statistic.series <- W2.series^2
    T.statistic <- rho*sigma.T^(-1)*n*(MU.A-MU.B)^2

    theta.hat <- constrOptim(c(1,qchisq(1-.05, df=1)),
                             grad = NULL,
                             method='Nelder-Mead',
                             f=LambdaCgammaMSE,
                             V.sigma=V.mat,
                             eta=eta.est,
                             rho=rho,
                             ui=rbind(c(1,0),
                                      c(-1,0),
                                      c(0,1),
                                      c(0,-1)),
                             ci=c(0,-10000,.01,-50))$par

    lambda.star <- ifelse(theta.hat[1]<0, 0, theta.hat[1])
    c.gamma <- theta.hat[2]

    mu.pool.series <- ((V.dr-Gamma.A.dr)*MU.A.series+ 
                         (V.a-Gamma.A.dr)*MU.B.series)/
      (V.a+V.dr-2*Gamma.A.dr)
    mu.pool <- ((V.dr-Gamma.A.dr)*MU.A+ 
                  (V.a-Gamma.A.dr)*MU.B)/
      (V.a+V.dr-2*Gamma.A.dr)

    mu.tap.series <- (T.statistic.series>c.gamma)*MU.A.series+
      (T.statistic.series<=c.gamma)*(MU.A.series+
                                       lambda.star*MU.B.series)/
      (1+lambda.star)
    mu.tap <- ifelse(T.statistic>c.gamma,
                     MU.A,(MU.A+lambda.star*MU.B)/
                       (1+lambda.star))

    return(list(bootstrap_samples=cbind(W1.series, W2.series,
                                        mu.pool.series, mu.tap.series,
                                        T.statistic.series),
                theta=theta.hat,
                output=c(mu.pool=mu.pool,
                         mu.tap=mu.tap,
                         'T'=T.statistic,
                         V.eff=V.eff,
                         V.a_eff=V.a_eff,
                         V.dr_eff=V.dr_eff,
                         W2.hat=W2.hat,
                         eta.est=eta.est)))
  }

  # --- Methods: pooled, using B, KH --------------------------------------
  Sigma.boot <- mapply(function(X,Y) {(n.A+n.B)*var(cbind(X,Y))},
                       X=mu.A.boot, Y=mu.dr.boot,SIMPLIFY = F)
  Sigma.boot.mean <- apply(Sigma.boot%>%simplify2array, 1:2, mean)

  results_normal <- mapply(TAP_variance,
                           V.mat=Sigma.boot,
                           MU.A=mu.A,
                           MU.B=mu.dr,
                           MU.A.series=mu.A.boot,
                           MU.B.series=mu.dr.boot,
                           rho=rho,SIMPLIFY = F)

  if(is.null(theta.mean))
  {
    theta.mean <- lapply(results_normal, function(x)x$theta)%>%
      simplify2array()%>%apply(1, mean)
  }
  T.statistic.mean <- lapply(results_normal, function(x)x$output['T'])%>%
    unlist()%>%mean

  Sigma.boot.usingB <- mapply(function(X,Y) {(n.A+n.B)*var(cbind(X,Y))},
                              X=mu.A.boot, Y=mu.dr.usingB.boot,SIMPLIFY = F)
  Sigma.boot.usingB.mean <- apply(Sigma.boot.usingB%>%simplify2array, 1:2, mean)
  results_usingB <- mapply(TAP_variance,
                           V.mat=Sigma.boot.usingB,
                           MU.A=mu.A,
                           MU.B=mu.dr.usingB,
                           MU.A.series=mu.A.boot,
                           MU.B.series=mu.dr.usingB.boot,
                           rho=rho,SIMPLIFY = F)
  if(is.null(theta.usingB.mean))
  {
    theta.usingB.mean <- lapply(results_usingB, function(x)x$theta)%>%
      simplify2array()%>%apply(1, mean)
  }
  T.statistic.usingB.mean <- lapply(results_usingB, function(x)x$output['T'])%>%
    unlist()%>%mean

  Sigma.boot.KH <- mapply(function(X,Y) {(n.A+n.B)*var(cbind(X,Y))},
                          X=mu.A.boot, Y=mu.dr.KH.boot,SIMPLIFY = F)
  Sigma.boot.KH.mean <- apply(Sigma.boot.KH%>%simplify2array, 1:2, mean)
  results_KH <- mapply(TAP_variance,
                       V.mat=Sigma.boot.KH,
                       MU.A=mu.A,
                       MU.B=mu.dr.KH,
                       MU.A.series=mu.A.boot,
                       MU.B.series=mu.dr.KH.boot,
                       rho=rho,SIMPLIFY = F)
  if(is.null(theta.KH.mean))
  {
    theta.KH.mean <- lapply(results_KH, function(x)x$theta)%>%
      simplify2array()%>%apply(1, mean)
  }
  T.statistic.KH.mean <- lapply(results_KH, function(x)x$output['T'])%>%
    unlist()%>%mean

  mu.pool.boot <- lapply(results_normal,
                         function(x)x$bootstrap_samples[,'mu.pool.series'])
  mu.pool.usingB.boot <- lapply(results_usingB,
                                function(x)x$bootstrap_samples[,'mu.pool.series'])
  mu.pool.KH.boot <- lapply(results_KH,
                            function(x)x$bootstrap_samples[,'mu.pool.series'])

  mu.tap.boot <- lapply(results_normal,
                        function(x)x$bootstrap_samples[,'mu.tap.series'])
  mu.tap.usingB.boot <- lapply(results_usingB,
                               function(x)x$bootstrap_samples[,'mu.tap.series'])
  mu.tap.KH.boot <- lapply(results_KH,
                           function(x)x$bootstrap_samples[,'mu.tap.series'])

  mu.pool <- ((Sigma.boot.mean[2,2]-Sigma.boot.mean[2,1])*mu.A+
                (Sigma.boot.mean[1,1]-Sigma.boot.mean[1,2])*mu.dr)/
    (Sigma.boot.mean[2,2]-Sigma.boot.mean[2,1]+
       Sigma.boot.mean[1,1]-Sigma.boot.mean[1,2])
  mu.tap <- ifelse(T.statistic.mean>theta.mean[2],
                   mu.A,(mu.A+theta.mean[1]*mu.dr)/
                     (1+theta.mean[1]))

  mu.pool.usingB <- ((Sigma.boot.usingB.mean[2,2]-Sigma.boot.usingB.mean[2,1])*mu.A+
                       (Sigma.boot.usingB.mean[1,1]-Sigma.boot.usingB.mean[1,2])*mu.dr)/
    (Sigma.boot.usingB.mean[2,2]-Sigma.boot.usingB.mean[2,1]+
       Sigma.boot.usingB.mean[1,1]-Sigma.boot.usingB.mean[1,2])

  mu.tap.usingB <- ifelse(T.statistic.usingB.mean>theta.usingB.mean[2],
                          mu.A,(mu.A+theta.usingB.mean[1]*mu.dr)/
                            (1+theta.usingB.mean[1]))

  mu.pool.KH <- ((Sigma.boot.KH.mean[2,2]-Sigma.boot.KH.mean[2,1])*mu.A+
                   (Sigma.boot.KH.mean[1,1]-Sigma.boot.KH.mean[1,2])*mu.dr)/
    (Sigma.boot.KH.mean[2,2]-Sigma.boot.KH.mean[2,1]+
       Sigma.boot.KH.mean[1,1]-Sigma.boot.KH.mean[1,2])

  mu.tap.KH <- ifelse(T.statistic.KH.mean>theta.KH.mean[2],
                      mu.A,(mu.A+theta.KH.mean[1]*mu.dr)/
                        (1+theta.KH.mean[1]))

  Sigma.boot.all <- mapply(function(mu.A,mu.dr,mu.pool,mu.tap) {(n.A+n.B)*
      var(cbind(mu.A,mu.dr,mu.pool,mu.tap))},
      mu.A=mu.A.boot, mu.dr=mu.dr.boot,
      mu.pool=mu.pool.boot, mu.tap=mu.tap.boot,
      SIMPLIFY = F)%>%simplify2array()%>%apply(1:2, mean)

  Sigma.boot.usingB.all <- mapply(function(mu.A,mu.dr,mu.pool,mu.tap) {(n.A+n.B)*
      var(cbind(mu.A,mu.dr,mu.pool,mu.tap))},
      mu.A=mu.A.boot, mu.dr=mu.dr.usingB.boot,
      mu.pool=mu.pool.usingB.boot, mu.tap=mu.tap.usingB.boot,
      SIMPLIFY = F)%>%simplify2array()%>%apply(1:2, mean)

  Sigma.boot.KH.all <- mapply(function(mu.A,mu.dr,mu.pool,mu.tap) {(n.A+n.B)*
      var(cbind(mu.A,mu.dr,mu.pool,mu.tap))},
      mu.A=mu.A.boot, mu.dr=mu.dr.KH.boot,
      mu.pool=mu.pool.KH.boot, mu.tap=mu.tap.KH.boot,
      SIMPLIFY = F)%>%simplify2array()%>%apply(1:2, mean)

  # --- Adaptive confidence interval (ACI) --------------------------------
  ACI <- function(results_variance, alpha.ci,
                  tau = 10^-seq(100,50, length.out = 10),
                  mu.tap=mu.tap)
  {
    alpha.tilde.ci <- alpha.ci/2
    K <- length(results_variance)

    W2.hat.mean <- lapply(results_variance,
                          function(x)x$output['W2.hat'])%>%
      unlist()%>%mean()

    mu.1 <- lapply(results_variance, function(x)mean(x$bootstrap_samples[,'W1.series']))%>%
      unlist()%>%mean()
    mu.2 <- lapply(results_variance, function(x)mean(x$bootstrap_samples[,'W2.series']))%>%
      unlist()%>%mean()

    V.eff <- lapply(results_variance, function(x)x$output['V.eff'])%>%unlist()%>%
      mean()
    V.a_eff <- lapply(results_variance, function(x)x$output['V.a_eff'])%>%unlist()%>%
      mean()
    V.dr_eff <- lapply(results_variance, function(x)x$output['V.dr_eff'])%>%unlist()%>%
      mean()

    lambda <- lapply(results_variance, function(x)x$theta[1])%>%unlist()%>%
      mean()
    c.gamma <- lapply(results_variance, function(x)x$theta[2])%>%unlist()%>%
      mean()

    W2.hat <- lapply(results_variance, function(x)x$output['W2.hat'])%>%unlist()%>%
      mean()

    W2.inf <- W2.hat+ qnorm(alpha.tilde.ci/2)
    W2.sup <- W2.hat+ qnorm(1-alpha.tilde.ci/2)

    mu.tap.projection.W2.method <- function(mu,
                                            alpha.tilde)
    {
      W1.series.new <- rnorm(500, mean = mu.1)
      W2.series.new <- rnorm(500, mean = mu.2)
      W2.series.project <- rnorm(500, mean = mu)
      mu.tap.series <- -V.eff^(1/2)*W1.series.new+
        (1+lambda)^(-1)*(V.a_eff^(1/2)-lambda*V.dr_eff^(1/2))*W2.series.new+
        (W2.series.new^2>=c.gamma&W2.series.new^2<log(log(n)))*
        (1+lambda)^(-1)*lambda*(V.dr_eff^(1/2)+V.a_eff^(1/2))*W2.series.new+
        (W2.series.project^2>log(log(n)))*
        (1+lambda)^(-1)*lambda*(V.dr_eff^(1/2)+V.a_eff^(1/2))*W2.series.project
      mu.tap-quantile(mu.tap.series, alpha.tilde)/sqrt(n)
    }

    mu.tap.inf1 <- optimize(
      function(x)mu.tap.projection.W2.method(x,
                                             1-alpha.tilde.ci/2),
      interval=c(W2.inf, W2.sup),
      maximum=FALSE)$objective

    mu.tap.sup1 <- optimize(
      function(x)mu.tap.projection.W2.method(x,
                                             alpha.tilde.ci/2),
      interval=c(W2.inf, W2.sup),
      maximum=TRUE)$objective

    ACI.bound_based <- function(tau,
                                selected.index=1:K)
    {
      mu.tap.inf2.series <- c()
      mu.tap.sup2.series <- c()
      for (kk in selected.index) {
        W1.series <- results_variance[[kk]]$bootstrap_samples[,'W1.series']
        W2.series <- results_variance[[kk]]$bootstrap_samples[,'W2.series']
        T.statistics <- results_variance[[kk]]$output['T']
        mu.1 <- mean(W1.series)
        mu.2 <- mean(W2.series)
        V.eff <- results_variance[[kk]]$output['V.eff']
        V.a_eff <- results_variance[[kk]]$output['V.a_eff']
        V.dr_eff <- results_variance[[kk]]$output['V.dr_eff']
        mu.pool <- results_variance[[kk]]$output['mu.pool']
        lambda <- results_variance[[kk]]$theta[1]
        c.gamma <- results_variance[[kk]]$theta[2]

        mu_tap_upper_bound <- -V.eff^(1/2)*W1.series+
          sapply(W2.series, function(x)
          {
            optimize(
              function(mu)mu.tap.bound.W2.method(x,mu,
                                                 lambda, c.gamma,
                                                 tau*log(log(n))),
              interval=c(-100, 100),
              maximum=TRUE)$objective
          })

        mu_tap_inf_bound <- -V.eff^(1/2)*W1.series+
          sapply(W2.series, function(x)
          {
            optimize(
              function(mu)mu.tap.bound.W2.method(x,mu,
                                                 lambda, c.gamma,
                                                 tau*log(log(n))),
              interval=c(-100, 100),
              maximum=FALSE)$objective
          })
        mu.tap.inf2 <- mu.tap-quantile(mu_tap_upper_bound, 1-alpha.tilde.ci)/sqrt(n)
        mu.tap.sup2 <- mu.tap-quantile(mu_tap_inf_bound, alpha.tilde.ci)/sqrt(n)

        mu.tap.inf2.series <- c(mu.tap.inf2.series, mu.tap.inf2)
        mu.tap.sup2.series <- c(mu.tap.sup2.series, mu.tap.sup2)
      }
      cbind(mu.tap.inf2=mu.tap.inf2.series,
            mu.tap.sup2=mu.tap.sup2.series)
    }

    groud_truth <- lapply(mu.A.boot, mean)%>%unlist()
    tau.double_boot <-  which.max(sapply(tau, 
                                         function(x)
                                         {
                                           selected.index <- sample(1:K,ceiling(K/3))
                                           ACI_bounds <- ACI.bound_based(tau = x,
                                                                         selected.index)
                                           mean(ACI_bounds[,1]<groud_truth[selected.index]&
                                                  groud_truth[selected.index]<ACI_bounds[,2])
                                         })>(1-alpha.ci)) 
    mu.tap.bounds.fixed <- ACI.bound_based(tau = 1)%>%apply(2, mean)
    mu.tap.bounds.double <- ACI.bound_based(tau = tau[tau.double_boot])%>%apply(2, mean)
    mu.tap.inf2 <- mu.tap.bounds.fixed[1]
    mu.tap.sup2 <- mu.tap.bounds.fixed[2]

    mu.tap.inf3 <- mu.tap.bounds.double[1]
    mu.tap.sup3 <- mu.tap.bounds.double[2]
    return(c(lower1=mu.tap.inf1,
             upper1=mu.tap.sup1,
             lower2=mu.tap.inf2,
             upper2=mu.tap.sup2,
             lower3=mu.tap.inf3,
             upper3=mu.tap.sup3))
  }

  ACI_normal <- ACI(results_normal, alpha.ci = .05,
                    tau = tau, mu.tap = mu.tap)
  names(ACI_normal) <- c('mu.tap_inf1','mu.tap_sup1',
                         'mu.tap_inf2','mu.tap_sup2',
                         'mu.tap_inf3','mu.tap_sup3')
  ACI_usingB <- ACI(results_usingB, alpha.ci = .05,
                    tau = tau, mu.tap = mu.tap.usingB)
  names(ACI_usingB) <- c('mu.tap.usingB_inf1','mu.tap.usingB_sup1',
                         'mu.tap.usingB_inf2','mu.tap.usingB_sup2',
                         'mu.tap.usingB_inf3','mu.tap.usingB_sup3')
  ACI_KH <-ACI(results_KH, alpha.ci = .05,
               tau = tau, mu.tap = mu.tap.KH)
  names(ACI_KH) <- c('mu.tap.KH_inf1','mu.tap.KH_sup1',
                     'mu.tap.KH_inf2','mu.tap.KH_sup2',
                     'mu.tap.KH_inf3','mu.tap.KH_sup3')

  regular_inf <- function(mu.boot, alpha.ci = .05)
  {
    mean(lapply(mu.boot, function(x)mean(x)-qnorm(1-alpha.ci/2)*sd(x))%>%
           unlist())}
  regular_sup <- function(mu.boot, alpha.ci = .05)
  {
    mean(lapply(mu.boot, function(x)mean(x)+qnorm(1-alpha.ci/2)*sd(x))%>%
           unlist())}

  regular_infs <- lapply(list(mu.A_inf=mu.A.boot,
                              mu.dr_inf=mu.dr.boot,
                              mu.dr.usingB_inf=mu.dr.usingB.boot,
                              mu.dr.KH_inf=mu.dr.KH.boot,
                              mu.pool_inf=mu.pool.boot,
                              mu.pool.usingB_inf=mu.pool.usingB.boot,
                              mu.pool.KH_inf=mu.pool.KH.boot),
                         function(x)regular_inf(x, alpha.ci = .05))%>%data.frame()
  regular_sups <- lapply(list(mu.A_sup=mu.A.boot,
                              mu.dr_sup=mu.dr.boot,
                              mu.dr.usingB_sup=mu.dr.usingB.boot,
                              mu.dr.KH_sup=mu.dr.KH.boot,
                              mu.pool_sup=mu.pool.boot,
                              mu.pool.usingB_sup=mu.pool.usingB.boot,
                              mu.pool.KH_sup=mu.pool.KH.boot),
                         function(x)regular_sup(x, alpha.ci = .05))%>%data.frame()

  return(
    list(
      point.est=data.frame(mu.A=mu.A,
                           mu.dr=mu.dr,
                           mu.eff=mu.pool,
                           mu.eff.usingB=mu.pool.usingB,
                           mu.eff.KH=mu.pool.KH,
                           mu.tap=mu.tap,
                           mu.tap.usingB=mu.tap.usingB,
                           mu.tap.KH=mu.tap.KH),
      Test.stat=data.frame(T.statistic = T.statistic.mean,
                           T.statistic.usingB = T.statistic.usingB.mean,
                           T.statistic.KH = T.statistic.KH.mean),
      CIs=c(regular_infs,
            regular_sups,
            ACI_normal,
            ACI_usingB,
            ACI_KH),
      theta.hat=data.frame(theta.mean,
                           theta.usingB.mean,
                           theta.KH.mean),
      Sigma.boot=Sigma.boot.all,
      Sigma.boot.usingB = Sigma.boot.usingB.all,
      Sigma.boot.KH = Sigma.boot.KH.all
    )
  )
}
