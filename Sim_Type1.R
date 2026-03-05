
install.packages("purrr")
library("purrr")

source("Ecomp-Rfns.R")


# Set parameters
n = 1000 # sample size
m = 1000 # no of iterations
avec = 0.95 # probability level
delay = 10

# eps = (1:(n+delay)) / (n+delay) # linear trend
# 
# # calculate VaR
# VaRtrue = qnorm(avec) * (1 + eps) # true VaR
# VaR10a = VaRtrue * (1 + .1) # over-report 10%
# VaR10b = VaRtrue * (1 - .1) # under-report 10%
# VaR20a = VaRtrue * (1 + .2) # over-report 20%
# VaR20b = VaRtrue * (1 - .2) # under-report 20%


eps = rdunif(n+delay, 5, -5) / 20 # random noice
# calculate VaR
VaRtrue = rep(qnorm(avec), n+delay) # true VaR
VaReps = VaRtrue + eps # true with noice
VaR20a = VaRtrue * (1 + .2) + eps  # over-report 10%
VaR20b = VaRtrue * (1 - .2) + eps # under-report 10%
VaR50a = VaRtrue * (1 + .5) + eps # over-report 20%
VaR50b = VaRtrue * (1 - .5) + eps # under-report 20%
VaRout = rbind(VaR50b, VaR20b, VaReps, VaR20a, VaR50a, VaRtrue)

# Simulate normal distribution
set.seed(233)
simnorm = matrix(rnorm((n+delay) * m), nrow = m, ncol = n+delay)


# iid + linear trend
# y = sweep(simnorm, 2, (1+eps), "*")
y = simnorm

# plot
plot(1:(n+delay), y[2,], type="l")



#---------------------------------------#
#
#      Comparative E-backtest
#
#---------------------------------------#

final_c = .5 # tuning parameter c for calculation of betting processes

nm=6 ; n=length(VaRtrue)

# Count number of rejections for different thresholds
numrej2 <- numrej5 <- numrej10 <- matrix(0, nrow = nm, ncol = nm)

M <- max(abs(y))


for (iter in 1:m) {
  
  if (iter %% 50 == 0) {
    print(iter)
  }
  
# Matrices to store score values
smatVaR1 <-smatVaR0 <- matrix(nrow=nm,ncol=n)
# matrix to store maximum and final e-values
egrel.mat.VaR.max <- egrel.mat.VaR <- matrix(nrow = nm, ncol = nm)



##------------------------------------------------------------------------------------------------##


for(i in 1:nm)
{
  smatVaR1[i,] <- sfVaR(r=VaRout[i,],x=y[iter,],a=avec,h=1) # 1-homogeneous scoring functions
  smatVaR0[i,] <- sfVaR(r=VaRout[i,],x=y[iter,],a=avec,h=0) # 0-homogeneous scoring functions
}

# Calculate the difference of VaR-score functions between methods

egrel.mat.VaR.max <- egrel.mat.VaR <- matrix(nrow = nm, ncol = nm) # matrix to store maximum and final e-values

options(warn = 0)

  for (i in 1:nm) {
    for (j in (1:nm)[-i]) {
      
      ##'@Notice: please adjust all the score functions including "sdiff.temp" & "sdiff.grel"
      
      sdiff.temp <- smatVaR1[i,] - smatVaR1[j,]
      #sdiff.temp <- smatVaR0[i,] - smatVaR0[j,]
      
      sf.code = 1 #score function = 1 or 0
      
      
      #M <- max(abs(y),abs(VaRout1[c(i,j),]))
      
      # E value calculation
      e.temp <- 1
      for (l in 2:n) {
        
        r = VaRout[i,l]
        r.star = VaRout[j,l]
        
        # gamma: upper bound of lambda
        if (sf.code == 1) {
          gamma.temp <- (1-avec)*(r-r.star) - (r<=r.star)*(min(r,-M) - min(r.star,-M)) - (r>r.star)*(min(r,M) - min(r.star,M))
        }else{
          gamma.temp <- (1-avec)*(log(r)-log(r.star)) - (r>r.star)*(log((min(max(r.star,M),r))/(r.star)))
        }
        
        ## New method
        c <- final_c
        if (gamma.temp>=0) {gamma <- Inf} else {gamma <- -c/gamma.temp} 
        
        
        
        ## GREL method
        
        sdiff.grel <- sfVaR(r=rep(r,l-1), x=y[iter, 1:(l-1)], a=avec, h=sf.code) - sfVaR(r=rep(r.star, l-1), x=y[iter, 1:(l-1)], a=avec, h=sf.code)
        
        
        ## New method
        if(sum(sdiff.grel) == 0) {lambda.grel.temp <- 0} else {
          lambda.grel.temp <- (sum(sdiff.grel))/(sum(sdiff.grel^2))}
        lambda.grel.temp <- min(max(lambda.grel.temp,0),gamma)
        
        
        e.temp <- e.temp*(1 + lambda.grel.temp*sdiff.temp[l])
      }
      
      egrel.mat.VaR.max[i,j] = max(e.temp)
    }}

egrel.mat.VaR <- t(egrel.mat.VaR.max)
numrej2 <- numrej2 + (egrel.mat.VaR >= 2) * 1
numrej5 <- numrej5 + (egrel.mat.VaR >= 5) * 1
numrej10 <- numrej10 + (egrel.mat.VaR >= 10) * 1
}

rej2 <- numrej2 / m
rej5 <- numrej5 / m
rej10 <- numrej10 / m

rej2
rej5
rej10