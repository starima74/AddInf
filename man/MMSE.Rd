\name{MMSE}
\alias{MMSE}
\title{
Minimum mean squared error estimation with additional information.
}
\description{
This function implements statistical estimation with the use of additional information. This additional information can come from multiple indepedent data sources. Each source is reported by a mean estimate and its variance. Then, this additional information is used together with an estimating procedure \code{theta.f} applied to the dataset \code{dd}. A new estimator minimizes mean squared error.
}
\usage{
MMSE(dd, theta.f, Add.Inf, nboots = 500, eig.cutoff = 0.9, ...)
}
\arguments{
  \item{dd}{
a data frame containing the variables used by function \code{theta.f}
}
  \item{theta.f}{
a function applied to the dataset \code{dd} to estimate a parameter of interest
}
  \item{Add.Inf}{
a data frame which consist of a columm of means (\code{Means}), a columm of variances (\code{Vars}), a vector of functions (\code{Functions}) and a vector of indicator of possible bias parameters (\code{Biases}). Functions column consists of functions used to calculate means in the additional data source; will be applied to the dataset \code{dd}. A indicator of possible bias parameter equal to 0 defines unbiased additional information, 1 corresponds to possible presence of bias.
}
  \item{nboots}{
a number of boostrap resamples to calculate variances and covariances
}
  \item{eig.cutoff}{
a percent cutoff for eigen values to approximate a matrix inverse
}
  \item{\dots}{
additional parameters. These parameters can be used by function \code{theta.f} or by additional functions.
}
}
\value{
  \item{Theta.Est}{The estimate of the parameter of interest without the use of addtional information}
  \item{Theta.Est.Var}{Variance covariance matrix of \code{Theta.Est}}
  \item{Theta.Hat}{New estimate of the parameter of interest with the use of addtional information}
  \item{Theta.Hat.Var}{Variance covariance matrix of \code{Theta.Hat}}
}
\references{
Tarima S.S. and Dmitriev Y.G. Statistical estimation with possibly incorrect model assumptions. The Bulletin of Tomsk State University: control, computing, informatics 2009, 3(8):87-99

Tarima S.S., Vexler A., Singh S., Robust mean estimation under a possibly incorrect log-normality assumption. Communications in Statistics - Simulation and Computation, 42: 316-326, 2013
}
\author{
Sergey Tarima (sergey.s.tarima@gmail.com)
}
\examples{
### data (an available dataset of patient ages)

dd <- data.frame(age = c(66,59,61,76,35,59,37,47,71,63,
                         41,56,49,58,48,79,87,50,57,58,32,
                         52,57,61,62,42,67,77,66,79))
                         
### the objective is to estimate proportion of pateints 
### younder than 40, 50, 60 and 70 years of age      
times1 <- c(40,50,60,70)

### this function estimates these proportions on a dataset d, 
### defined in its first argument
theta.f <- function(d,times) {
  ecdf(c(d)[[1]])(times)
}

### (Empty) Additional Information
Add.Info.Means <- list()
Add.Info.Vars <- list()
Add.Info.Functions <- list()
Add.Info.Biases <- list()

#### the first additional data source says that average of patients
### is 58.58 years of age with variance = 0.09, this information is possibly biased
Add.Info.Means[[1]] <- 58.58
Add.Info.Vars[[1]] <- 0.09
Add.Info.Functions[[1]] <- function(d) mean(d[[1]],na.rm = TRUE)
Add.Info.Biases[[1]] <- 1 # biased additional information

### create the data frame where the additional information will be saved
Add.Info <- data.frame(Means = rep(NA,1), 
                       Vars = rep(NA,1), 
                       Functions = rep(NA,1),
                       Biases = rep(NA,1))

### Adding the additional information to the data frame
Add.Info$Means = Add.Info.Means
Add.Info$Vars = Add.Info.Vars
Add.Info$Functions = Add.Info.Functions
Add.Info$Biases = Add.Info.Biases

### using both available dataset and additional information for estimation
res <- MMSE(dd, theta.f, Add.Info, nboots = 500, 
               eig.cutoff = 1, times = times1)

### calculating new estimates and pointwise confidence intervals
lo  <- res$Theta.Est - 1.96*sqrt(diag(res$Theta.Est.Var))
cdf <- res$Theta.Est
hi  <- res$Theta.Est + 1.96*sqrt(diag(res$Theta.Est.Var))

### calculating empirical (no additional information) estimates 
### and pointwise confidence intervals
theta_lo  <- res$Theta.Hat - 1.96*sqrt(diag(res$Theta.Hat.Var))
theta_cdf <- res$Theta.Hat
theta_hi  <- res$Theta.Hat + 1.96*sqrt(diag(res$Theta.Hat.Var))

### plot
plot(cdf ~ times1, type = "l", ylab="Proportions", 
          main = paste("Means and 95 percent confidence intervals\n", 
                   "with (dashed) and without (solid) additional information"), 
          xlab="Age", lty = 2, ylim= c(0,1))
lines(lo ~ times1, type = "l", lty = 2)
lines(hi ~ times1, type = "l", lty = 2)
points(cdf ~ times1)
points(lo ~ times1)
points(hi ~ times1)

lines(theta_cdf ~ times1, type = "l")
lines(theta_lo ~ times1, type = "l")
lines(theta_hi ~ times1, type = "l")
points(theta_cdf ~ times1)
points(theta_lo ~ times1)
points(theta_hi ~ times1)
}
