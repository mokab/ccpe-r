library(minpack.lm)
library(MASS)


curveFit <- function(x_centre, y_centre, a, x, y, z) {
  v <- 1
  t <- y + z
  #f <- function(v, x, a, x_centre, y_centre) {
  #  a*cos(x/v) + x_centre + a*sin(x/v) + y_centre
  #}
  #ff <- a*cos(x/v) + x_centre + a*sin(x/v) + y_centre
  result <- nlsLM(t ~ a*cos(x/v) + x_centre + a*sin(x/v) + y_centre, start = c(v=v), trace = FALSE)
  v <- coef(result)[1]
  th <- x / v
  y_new <- a*cos(th) + x_centre
  z_new <- a*sin(th) + y_centre
  Z_p <- rbind(x, y_new, z_new)
  resnorm <- sum(result$residuals^2)
  return(list(v = v, th = th, resnorm = resnorm, Z_p = Z_p))
}

HelixFit <- function(Z, vaxis) {
  x <- Z[,1]
  y <- Z[,2]
  z <- Z[,3]
  circleFit_output <- circleFit(cbind(y, z))
  x_centre <- circleFit_output$x_centre
  y_centre <- circleFit_output$y_centre
  a <- circleFit_output$r
  #cat('vaxis is x default\n')
  curveFit_output <- curveFit(x_centre, y_centre, a, x, y, z)
  v <- curveFit_output$par[1]
  th <- curveFit_output$par[2]
  resnorm <- curveFit_output$resnorm
  Z_p <- curveFit_output$Z_p
  return(list(x_centre = x_centre, y_centre = y_centre, Z_p = Z_p, a = a, v = v, th = th))
}


circleFit <- function(P) {
  per_error <- 0.1/100
  
  # initial estimates
  X <- colMeans(P)
  r <- sqrt(mean(rowSums((matrix(X, nrow = nrow(P), ncol = ncol(P), byrow = TRUE) - P)^2)))
  
  v_cen2points <- matrix(0, nrow = nrow(P), ncol = ncol(P))
  niter <- 0
  
  # looping until convergence
  while (niter < 1 || per_diff > per_error) {
    
    # vector from centre to each point
    v_cen2points[, 1] <- P[, 1] - X[1]
    v_cen2points[, 2] <- P[, 2] - X[2]
    
    # distacnes from centre to each point
    centre2points <- sqrt(rowSums(v_cen2points^2))
    
    # distances from edge of circle to each point
    d <- centre2points - r
    
    # computing 3x3 jacobean matrix J, and solvign matrix eqn.
    R <- v_cen2points / matrix(c(centre2points, centre2points), ncol = 2, byrow = TRUE)
    J <- cbind(-rep(1, nrow(P)), -R)
    D_rXY <- -1 * ginv(t(J) %*% J) %*% t(J) %*% d
    
    # updating centre and radius
    r_old <- r; X_old <- X
    r <- r + D_rXY[1]
    X <- X + D_rXY[2:3]
    
    # calculating maximum percentage change in values
    per_diff <- max(abs(c((r_old - r) / r, (X_old - X) / X))) * 100
    
    # prevent endless looping
    niter <- niter + 1
    if (niter > 1000) {
      stop("Convergence not met in 1000 iterations!")
    }
  }
  
  x_centre <- X[1]
  y_centre <- X[2]
  sq_error <- sum(d^2)
  
  return(list(x_centre = x_centre, y_centre = y_centre, r = r, sq_error = sq_error))
}

get_R <- function(Z, Y, sig) {
  N <- nrow(Z)
  K <- nrow(Y)
  R <- matrix(nrow=N, ncol=K)
  
  for (i in 1:N) {
    sum <- 0
    for (k in 1:K) {
      sum <- sum + exp(-sum((Z[i,] - Y[k,])^2) / sig)
    }
    for (j in 1:K) {
      R[i, j] <- exp(-sum((Z[i,] - Y[j,])^2) / sig) / sum
    }
  }
  
  return(R)
}

  
CCPE <- function(X,lambda,gam,sig) {
  pca <- prcomp(X, scale = TRUE, center = TRUE)
  W <- pca$rotation
  Z <- pca$x
  N <- nrow(Z)
  k <- N
  vaxis <- "x"
  Y <- Z
  
  iter <- 0
  for(i in 1:200) {
    cat(paste0('\r', i))
    HelixFitOut <- HelixFit(Z,vaxis)
    x_centre <- HelixFitOut$x_centre
    y_centre <- HelixFitOut$y_centre
    Z_p <- HelixFitOut$Z_p
    a <- HelixFitOut$a
    v <- HelixFitOut$v
    th <- HelixFitOut$th
    R <- get_R(Z,Y,sig)
    Tmat <- diag(diag(matrix(1, k, N) * R))
    Q <- ginv((1 + lambda + gam) * diag(N) - gam * R %*% ginv(Tmat) %*% t(R))
    A <- Z_p %*% Q %*% X
    svdout <- svd(A)
    U <- svdout$u
    S <- diag(svdout$d)
    V <- svdout$v
    I <- diag(3)
    W <- V %*% I %*% t(U)
    Z <- Q %*% (X %*% W + lambda * t(Z_p)) 
    Y <- ginv(Tmat) %*% R %*% Z
    error1 <- sum((X - Z %*% t(W))^2)
    error2 <- lambda * sum((Z - t(Z_p))^2)
    error3 <- gam * (sum(Z^2) - 2 * sum(R * Z %*% t(Y)) + sum(t(Y) %*% Tmat %*% Y))
    MSE_error <- error1 + error2 + error3
  }
  pseudotime <- Z[,1]
  return(pseudotime)
}


CCPE(X=data.matrix(iris[,-5]),lambda=1,gam=1,sig=1)