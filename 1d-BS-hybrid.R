# install.packages("QuantLib")
# install.packages("RQuantLib")
# library("RQuantLib")
library(MASS)

C <- function(coeff, x, ctrs, x_star=Inf,K){
  if (x<x_star) {
    sum(coeff* LHS_xlarge(x,ctrs))
  }
  else{
    K*exp(x)-K
  }
}

# When x > x*
LHS_xlarge <- function (x, ctrs){
  # input x > x^* 
  # output vector of length n
  exp(-(x-ctrs)^2)
}
RHS_xlarge <- function (x, K){
  K*exp(x)-K
}

# When x < x*
LHS_xsmall<- function (x, ctrs, sig, r , D0 , dt){
  d1 <- -2 * (x-ctrs)
  d2 <- -2 + 4 * (x-ctrs)^2
  An <- 1 - (  sig^2 /2 * d2 + (r-D0-sig^2 / 2) *d1 - r )*dt
  LHS_xlarge(x,ctrs) * An
}



LHS_xstar <- function(x, ctrs){
  d1 <- -2 * (x-ctrs)
  LHS_xlarge(x,ctrs) * d1
}

RHS_xstar <- function(x, K){
  K*exp(x)
}



# Construct centers for each x* in the region of x
# Achive equivalent comparibability


# Given a guess of x*, find a set of coeff
linsys_a2 <- function(x_star, x_mat, ctrs_mat, dt, dx, coeff_pre , x_star_pre, K,D0){
  # x_star_pre = 23
  # x_star = 30
  # coeff_pre = NULL
  x_all = x_mat[x_star,]
  ctrs_post = ctrs_mat[x_star,]
  ctrs_pre = ctrs_mat[x_star_pre,]
  if (is.null(coeff_pre)){
    C_pre = sapply(x_all, function(y){max((K*exp(y)-K), 0)})
  }else{
    C_pre = sapply(x_all, function(y){C(coeff_pre,y,ctrs_pre, x[x_star_pre], K)})
  }
  N = length(x_all)
  LHS <- c()
  RHS <- c()
  for (i in 2 : (N-1)){
    LHS = rbind(LHS, LHS_xsmall(x= x_all[i], ctrs_post, sig, r, D0,dt))
    RHS = append(RHS,C_pre[i]) # RHS_xsmall(i , C_pre)
  }
  # Boundary Condition (x -> -inf , c -> 0)
  # 0-order
  LHS = rbind(LHS, LHS_xlarge(x = x_all[1] , ctrs_post)) # x_all[1]
  RHS = append(RHS,0) 
  # 1-order
  LHS = rbind(LHS, LHS_xstar(x = x_all[1], ctrs_post))
  RHS = append(RHS,0)
  
  # At the Exercising boundary: x*; 
  # 0-order smoothness
  LHS = rbind(LHS, LHS_xlarge(x = x_all[N], ctrs_post))
  RHS = append(RHS,RHS_xlarge(x = x_all[N], K))
  
  # # 1-order smoothness
  # LHS = rbind(LHS, LHS_xstar(x = x_all[N], ctrs_post))
  # RHS = append(RHS,RHS_xstar(x = x_all[N], K))
  
  coeff_post = ginv(LHS)%*%RHS
  err1 <- abs(LHS_xlarge(x = x_all[N], ctrs_post)%*% coeff_post - RHS_xlarge(x = x_all[N], K))
  err2 <- abs(LHS_xstar(x = x_all[N], ctrs_post)%*% coeff_post - RHS_xstar(x = x_all[N], K))
  err  <- abs(err1) +abs(err2)
  
  sapply(x_all, function(y){C(coeff_post,y,ctrs_post, x[x_star], K)})
  # err  <- (LHS %*% coeff - RHS)
  # te <-  sum((err)^2)
  # mse  <- te/(nrow(LHS))
  return(list(lhs = LHS, rhs = RHS, coeff = coeff_post, err =err, err_order0=err1, err_order1 = err2))
}



# Plot the difference between theoretical and numerical results 
plot_compare_a2 <- function(xs, coeff_pre, x_star_pre, tau, K){
  library("NMOF")
  ctrs_pre = ctrs_mat[x_star_pre,]
  collocation = sapply(xs, function(y){C(coeff_pre,y,ctrs_pre,x_star=x[x_star_pre],K)} )
  theoretical = sapply(xs, function(x){vanillaOptionAmerican(S= K*exp(x), X=K, tau=tau, r, q =D0, v=sig^2, type = "call", greeks = FALSE, M = 400)})
  theoretical2 = sapply(xs, function(x){vanillaOptionEuropean(S= K*exp(x), X=K, tau=tau, r, q =D0, v=sig^2, type = "call", greeks = FALSE)})
  plot(xs,collocation, type="l", pch=21, lty=1,lwd= 2.5,col="black" ,axes = FALSE,ann=FALSE)
  lines(xs,theoretical, type="l", pch=22, lty=2, lwd =2.5,col="lightblue",axes=FALSE, ann=FALSE)
  lines(xs,theoretical2, type="l", pch=23, lty=2, lwd =2.5,cex =1, col="pink",axes=FALSE, ann=FALSE)
  lines(xs, sapply(c(K*exp(xs)-K), function(x){max(x,0)}),type="l", pch=24, lty=1, lwd =2,col="gray",axes=FALSE, ann=FALSE)
  # Make x axis using Mon-Fri labels
  axis(1, las=1, at =seq(min(xs),max(xs), by=0.1))
  axis(2, las=1, at=seq(min(floor(collocation)),max(collocation), by=5))
  box()
  title(xlab= "x  = log(S/K)", col.lab=rgb(0,0,0))
  title(ylab= "American Call Price", col.lab=rgb(0,0,0))
  title(main = "American Call: B-S Model", family = "serif")
  legend(-0.5,range(collocation)[2], c(paste0("r: ",r, ", sig: ",sig, ", T: ", T, ", D0: ", D0 ),"Collocation American","Theoretical American","Theoretical European","Call Payoff"), cex=0.8, 
         col=c("black","black","lightblue","pink","gray"), pch=c(15, 15,15,15,15), lty=c(0, 1,2,2,1),lwd=c(1,2.5,2.5,2.5,2))
  
  return ( list (e = sum((collocation-theoretical) ^2 )/ K) ) 
}


# Parameter Setting
r   <- 0.03
sig <- 0.4
D0  <- 0.07
K   <- 100
T 	<- 3

rho <- 1
a   <- -1
x_p <- 1
M   <- 10
dt  <- T / M
dx  <- sqrt(T)/M
N   <- ceiling((x_p-a)/dx)
x   <- append(seq(from = a, to = 0, length = N/2), seq(from = dx, to = x_p, length = N/2))
N   <- length(x)
# ctrs <- c(seq(from = a*2, x_p*2, length = N))


ctrs_mat <- c()
for (i in 1:N) {
  ctrs_mat <- rbind(ctrs_mat,c(seq(from = a*1.1, x[i]*1.1, length = N)) )
}

x_mat <- c()
for (i in 1:N) {
  x_mat <- rbind(x_mat,c(seq(from = a, x[i], length = N)) )
}

xs = seq(-0.5,0.5 ,length=11)
# xs
# xs = log(seq(from =40, to = 120, by =10)/100)
x_star_pre = match(0,x)
coeff_pre = NULL
ptm <- proc.time()
# i=31
# x[i]
for (j in 1:M){
  # Find x* for each tau
  se     = 10^10
  for (i in 1:N){
    if (x[i] >= 0 ){
      try = linsys_a2( i , x_mat, ctrs_mat, dt, dx, coeff_pre, x_star_pre, K, D0)
      if (try$err_order1< se ){
        se = try$err_order1
        output = try
        x_star_pre1 = i
        coeff_pre1 = output$coeff
        C_pre1 = sapply(xs,function(y){C(coeff_pre,y,ctrs_mat[x_star_pre1],x_star = x[x_star_pre1],K)})
        # cat(C_pre1,"\n")
        # cat("tau",(j*dt),"x_star: ",x[x_star_pre1], "err", output$err,"err1",output$err_order0, "err2", output$err_order1, "\n")
      }
      # cat("tau",(j*dt),"x_star: ",x[x_star_pre1],  "err0", output$err_order0, "\n")
    }
  }
  x_star_pre = x_star_pre1
  coeff_pre  = coeff_pre1
  # coeff = output$coeff
  cat("tau",(j*dt),"x_star: ",x[x_star_pre], "err", output$err,"err1",output$err_order0, "err2", output$err_order1, "\n")
  # C_pre <-  append(sapply(x[1:x_star], function(x){(C(output$coeff,x,ctrs)) } ) , if(x_star<N){RHS_xlarge(x[(x_star+1):N], K)}else{NULL} )
}
time_lagged = proc.time() - ptm; cat("Time Used: ", time_lagged)
C_pre = sapply(xs, function(y){C(coeff_pre,y,ctrs_mat[x_star_pre,], x[x_star_pre], K)})
plot_compare_a2(xs,  coeff_pre, x_star_pre, tau=T, K)
