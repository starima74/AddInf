\name{MMSE2}
\alias{MMSE2}
\title{
Minimum mean squared error estimation with additional information.
}
\description{
This function implements statistical estimation with the use of additional information. This additional information can come from multiple indepedent data sources. Each source is reported by a mean estimate and its variance. In contrast with the \code{MMSE} function where the dataset is formally used, the \code{MMSE2} uses only a variance-covariance matrix \code{K}, and it also assumes a that the vector of biases \code{c} is known.
}
\usage{
  MMSE(K, c, Add.Info, thetahat, thetahat.length, eig.cutoff = 0.9)
}
\arguments{
  \item{K}{
a variance covariance matrix. See the example below}
  \item{c}{
 a vector of biases 
}
  \item{Add.Inf}{
a data frame with additional information, which consists of a columm of means (\code{Means}) and a columm of variances (\code{Vars}). The vector of indicators of possible biases \code{Biases} uses 0 to define unbiased additional information and 1 for possible bias.
}
\item{thetahat}{
 an estimated value for a parameter of interest
}  
\item{betahat}{
 an estimated value for an auxiliary parameter
}
  \item{eig.cutoff}{
a percent cutoff for eigen values to approximate a matrix inverse
}
}
\value{
  \item{Theta.Est}{An estimate of the parameter of interest without the use of addtional information}
  \item{Theta.Est.Var}{A variance covariance matrix of \code{Theta.Est}}
  \item{Theta.Hat}{A new estimate of the parameter of interest with the use of addtional information}
  \item{Theta.Hat.Var}{A variance covariance matrix of \code{Theta.Hat}}
}
\references{
Tarima S.S. and Dmitriev Y.G. Statistical estimation with possibly incorrect model assumptions. The Bulletin of Tomsk State University: control, computing, informatics 2009, 3(8):87-99

Tarima S.S., Vexler A., Singh S., Robust mean estimation under a possibly incorrect log-normality assumption. Communications in Statistics - Simulation and Computation, 42: 316-326, 2013

Sergey Tarima and Kadam Patel. The Association between Gestational Age and 3rd Grade Standardized Reading Score. MCW Biostatistics technical report #68. https://www.mcw.edu/-/media/MCW/Departments/Biostatistics
}
\author{
Sergey Tarima (sergey.s.tarima@gmail.com) and Roman Vygon ()
}
\examples{

### A variance covariance matrix of adjuated and unadjusted
### estimates of differences in reading scores between children 
### born at the gestational age of 37 weeks and full term 
### birth children (39-41 weeks)

K <- matrix(c(0.0467408718, 0.0499675797, 0.0499675797, 0.0768166235),2,2)

### an estimate of the parameter of interest (the adjusted 
### estimate of differences (in percents) in reading scores
### between children born at the gestational age of 37 weeks 
### and full term children
thetahat <- 0.173

### an estimate of the auxiliary parameter (the unadjusted 
### estimate of differences in reading scores between children 
### born at the gestational age of 37 weeks and full term children
betahat <- -0.952

### (Empty) Additional Information (on the auxiliary parameter)
Add.Info.Means <- list()
Add.Info.Vars <- list()
Add.Info.Biases <- list()

### the additional data source says the unadjuted difference 
### is -0.74 with a standard error = 0.098 (variance = 0.009604)
### this information may be biased as it came from 
### a different population

Add.Info.Means[[1]] <- -0.74
Add.Info.Vars[[1]] <- 0.009604
Add.Info.Biases[[1]] <- 1 # biased additional information

### create the data frame where the additional information will be saved
Add.Info <- data.frame(Means = rep(NA,1), 
                       Vars = rep(NA,1), 
                       Biases = rep(NA,1))

### Adding the additional information to the data frame
Add.Info$Means = Add.Info.Means
Add.Info$Vars = Add.Info.Vars
Add.Info$Biases = Add.Info.Biases

### estimation based on both data sources 
res <- MMSE2(K, -0.05, Add.Info, thetahat, betahat,
            eig.cutoff = 1)

### the new estimate and its confidence interval
lo  <- res$Theta.Est - 1.96*sqrt(res$Theta.Est.MSE)
est <- res$Theta.Est
hi  <- res$Theta.Est + 1.96*sqrt(res$Theta.Est.MSE)
c(lo, est, hi)

### the empirical (no additional information) estimate 
### and its confidence interval
theta_lo  <- res$Theta.Hat - 1.96*sqrt(res$Theta.Hat.Var)
theta_est <- res$Theta.Hat
theta_hi  <- res$Theta.Hat + 1.96*sqrt(res$Theta.Hat.Var)
c(theta_lo, theta_est, theta_hi)

### a sensitivity analysis for a range of bais values 

ra <- -10:10/11
n <- length(ra) 
CIs <- data.frame(LO=rep(NA,n),EST=rep(NA,n),HI=rep(NA,n))
s=0
for(delta in ra) { 
 s=s+1
 res <- MMSE2(K, delta, Add.Info, thetahat, betahat,
            eig.cutoff = 1)
 lo  <- res$Theta.Est - 1.96*sqrt(res$Theta.Est.MSE)
 est <- res$Theta.Est
 hi  <- res$Theta.Est + 1.96*sqrt(res$Theta.Est.MSE)
 CIs[s,] <- c(lo, est, hi)
}

plot(ra,CIs$EST,ylim=c(-0.5,1),type="l",lwd=3,
    ylab="Estimated difference and 95% CI",xlab="bias")
lines(ra,CIs$LO)
lines(ra,CIs$HI)
lines(ra,rep(theta_lo,n),lty=2)
lines(ra,rep(theta_hi,n),lty=2)
lines(ra,rep(theta_est,n),lty=2,lwd=3)
lines(ra,rep(0,n),lty=3,lwd=1)
}
