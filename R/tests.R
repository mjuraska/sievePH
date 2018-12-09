### TESTING CODE 

library(survival)
library(np)
dataDir <- "T:/vaccine/rtss_malaria_sieve/Stephanie's work/AMP"
dataB <- read.csv(file.path(dataDir, "catnap_vrc01_neut_b.csv"))
dataC <- read.csv(file.path(dataDir, "catnap_vrc01_neut_c.csv"))
dataBandC <- read.csv(file.path(dataDir, "catnap_vrc01_neut_all.csv"))
codeDir <- "T:/vaccine/rtss_malaria_sieve/Stephanie's work/AMP/sievePower"
source(file.path(codeDir, "markHR.R"))
source(file.path(codeDir, "covEst.R"))
source(file.path(codeDir, "markDensSim.R"))
Np <- 900
np <- 34
data <- dataB
randomRatio = 2
taumax=80
lambdaT <- (log(1 - (1 + 0.1*Np/np)*(np/Np)))/(-taumax*(1 + 0.1*Np/np))
lambdaC <- 0.1*Np/np*lambdaT
VE=0.7
IC <- "IC50"
mark <- data[,paste0(tolower(IC),".geometric.mean.imputed.log10")]
densbw <- npudensbw(~ mark, ckertype="epanechnikov")

dens <- npudens(densbw)
coeff <- getAlphaBetaGamma(VE, dens, "mark")
alpha=coeff$alpha
beta=coeff$beta
gamma=coeff$gamma
varName="mark"

Z <- c(rep(0, Np), rep(1, randomRatio*Np))        # treatment group
T0 <- rexp(Np, lambdaT)                           # failure times for placebo
T1 <- rexp(randomRatio*Np, lambdaT*exp(gamma))    # failure times for vaccine
T <- c(T0,T1)                                     # failure times
C <- rexp((randomRatio+1)*Np, lambdaC)            # censoring times
X <- pmin(T,C, taumax)                            # observed time: minimum of failure, censoring, and study time
d <- ifelse(T <= pmin(C,taumax),1,0)              # failure indicator (0 if censored)
nInf0 <- sum(d*(1-Z))                             # number of infected in placebo group
nInf1 <- sum(d*Z)                                 # number of infected in vaccine group

Xpoints <- seq(-3,3,len=25000)                    # fine grid of points ranging from -3 to 3 to be sampled from
prob0 <- f0(Xpoints, dens, varName)               # sampling probability for placebos, using nonparametric density estimates
prob1 <- f1(Xpoints, dens, varName, alpha, beta)  # sampling probabiliy for vaccinees, using nonparametric density estimates 
set.seed(1)
V0 <- sample(Xpoints, size=Np, prob=prob0, replace=TRUE)        # sample with replacement with probability prob0 to simulate mark in placebos 
set.seed(2)
V1 <- sample(Xpoints, size=randomRatio*Np, prob=prob1, replace=TRUE)  # sample with replacement with probability prob1 to simulate mark in vaccinees
V <- c(V0, V1)                                    # mark variable
V <- ifelse(V < log10(0.00076), log10(0.00076), ifelse(V <= log10(50), V, log10(50)))  # mark variable with extreme values censored

dRatio <- densRatio(V[d==1],Z[d==1])  

