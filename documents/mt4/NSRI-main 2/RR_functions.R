
RR1 <- function(df){
  
  y <- df$Freq
  
  n <- sum(y)
  
  pobs <- y / n

  P <- matrix(c(5, 1, 1, 5) / 6, 2, 2)
  
  phat <- solve(P) %*% pobs
  
  varp <- prod(pobs)/(n*(P[2,2]-P[2,1])^2)
  
  data.frame(df, phat = phat, min95 = phat - 1.96*sqrt(varp), max95 = phat + 1.96*sqrt(varp))

}

RR2 <- function(df){
  
  df <- na.omit(df)
  
  y  <- df$Freq
  
  n  <- sum(y)
  
  pobs <- y / n
  
  P <- matrix(c(5, 1, 1, 5) / 6, 2, 2)
  
  Q <- P %x% P
  
  phat <- solve(Q) %*% pobs
  
  varp <- (n - 1)^-1 * solve(Q) %*% (diag(pobs) -  pobs %*% t(pobs)) %*% t(solve(Q))
  
  data.frame(df, phat = phat, min95 = phat - 1.96*sqrt(diag(varp)), max95 = phat + 1.96*sqrt(diag(varp)))
  
}

RRlogreg <- function(df){
  
  df <- na.omit(df)
  
  y <- c(df[[1]])
  
  x <- model.matrix(~ ., df)[, -2]
  

  P <- matrix(c(5, 1, 1, 5) / 6, 2, 2)
  
  logl <- function(b){
    
    prev  <- exp(x %*% b)/(1 + exp(x %*% b))
    
    p     <- t(cbind(1 - prev, prev))
    
    logl  <- log((P[1, ] %*% p)^(1 - y)) + log((P[2, ] %*% p)^y)
    
    return(-sum(logl) + 0.01*sum(b[1]^2))
  }
  
  fit  <- optim(c(-6, rep(0, ncol(x) - 1)), logl, hessian = T, method = "BFGS")
  
  bhat <- fit$par
  
  bvar <- diag(solve(fit$hessian))

  data.frame(pars = colnames(x), bhat, bvar)
  
}

RRpropodds <- function(df){
  
  
  df <- na.omit(df)
  
  y  <- df[, 1:2]
  
  sy <- rowSums(y)
  
  n  <- nrow(y)
  
  M  <- ncol(y)
  
  x <- model.matrix(~ ., df)[, -(1:3), drop = F]
    
  P <- matrix(c(5, 1, 1, 5) / 6, 2, 2)
  
  Q <- P %x% P
  
  Q[2, ] <- colSums(Q[2:3, ])
  
  Q <- Q[-3, -3]
  
  
  logl <- function(theta){
    
    alpha     <- matrix(theta[1:M], n, M, byrow = T)
    
    b         <- theta[(M + 1):length(theta)]
    
    if(ncol(df) > M){
      
    xb        <- x %*% b
    
    cum.probs <- exp(alpha - c(xb)) / (1 + exp(alpha - c(xb)))
      
    } else {
      
      cum.probs <- exp(alpha) / (1 + exp(alpha))
    }
    
    probs     <- cbind(cum.probs, 1) - cbind(0, cum.probs)
    
    logl  <- 0
    
    for(j in 0:M){
      
      logl <- logl - sum(log(Q[j + 1, ] %*% t(probs[sy == j, ])))
      
    }
    
    return(logl + sum(b^2))
  }
  
  if(ncol(df) > M){
    
    par0 <- c(2, 5, rep(0, ncol(x))) 
    
  } else {
    
    par0 <- c(2, 5)
    
  }
  
  fit  <- optim(par0, logl, hessian = T, method = "BFGS")
  
  bhat  <- fit$par
  
  bvar  <- diag(solve(fit$hessian))
  
  data.frame(pars = c(paste0("alpha", 1:M), if(ncol(df) > M)colnames(x)), bhat, bvar)
  
}


































