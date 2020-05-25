# info about distribution
mu = 0
sigma = 1

gendata <- function(n) {
  data = mean(ifelse(runif(n) < .5, -1, 1))
}

nvals = c(10,100,1000, 10000)
data = t(sapply(1:10000, function(y) sapply(nvals, function(x) gendata(x))))

# part a: the deviation converges to 0 as n grows. 
matplot(log10(nvals), data[1,]-mu, type = "l")

# part b: for any epsilon, the probability of the sample mean straying
# from the expected mean by more than epsilon always tends to 0 
# as n tends to infinity. By taking the sample mean over N (10000) independent
# trials, we get a good estimate of this probability (again by LLN), i.e.
# \frac{1}{N}\sum_{i=1}^N 1_{|X^i_n-\mu| > epsilon}
# tends to P(|X^i_n-\mu| > epsilon) as N gets large.
epsilon = c(.5, .1, .05)
indicator = t(sapply(epsilon, function(x) sapply(1:4, function(i) length(which(abs(data[,i] - mu) > x))/nrow(data))))
matplot(log10(nvals), indicator[1,], type="l")
matplot(log10(nvals), indicator[2,], type="l")
matplot(log10(nvals), indicator[3,], type="l")

# part c: as n grows, the histograms and Q-Q plots looks "closer" to normal.
# in particular, the histogram converges in distribution to a standard normal
# distribution, i.e. the cdfs tend toward a standard normal cdf.
hist(sqrt(nvals[1])*(data[,1]-mu)/sigma)
hist(sqrt(nvals[2])*(data[,2]-mu)/sigma)
hist(sqrt(nvals[3])*(data[,3]-mu)/sigma)
hist(sqrt(nvals[4])*(data[,4]-mu)/sigma)

# the QQ plot illustrates this convergence as well - the deviation from
# the reference line (at 45 degrees) tends to get smaller as n is larger
# so the two distributions are more and more similar as n grows. See
# http://www.itl.nist.gov/div898/handbook/eda/section3/qqplot.htm for more.
qqnorm(sqrt(nvals[1])*(data[,1]-mu)/sigma)
qqnorm(sqrt(nvals[2])*(data[,2]-mu)/sigma)
qqnorm(sqrt(nvals[3])*(data[,3]-mu)/sigma)
qqnorm(sqrt(nvals[4])*(data[,4]-mu)/sigma)

# part d: we note that the sample distribution in the statement of CLT 
# does NOT converge in PROBABILITY to the standard normal distribution.
# This is captured by the following plot which shows that regardless of n,
# the probability that the random variable sqrt(n)*(Xbar^i_n-mu)/sigma is at least
# a distance epsilon = 0.001 from a standard normal RV, Y_i, is very close to 1
# regardless of n. Almost 100% of the normalized sample data points stray from 
# i.i.d. standard normal Y_i regardless of n.
normals = rnorm(10000)
epsilon = .001
indicator2 = sapply(1:4, function(i) length(which(abs(sqrt(nvals[i])*(data[,i]-mu)/sigma - normals) > epsilon))/nrow(data))
matplot(log10(nvals), indicator2, type = "l")

# Q2 part d: note that since theta = 2, we have the cdf is
# F(x) = int_1^\infty(y^(-2)dy) = (-1/y from 1 to x) = 1 - 1/x
# thus F^(-1)(x) = 1 / (1-x). This allows us to generate i.i.d. X_i ~ P_theta
# in line 65 below.
CItest <- function(){
  n = 100
  theta = 2
  z = 1.96
  q1samples = 1/(1-runif(n))
  theta.hat = n/sum(log(q1samples))+1
  left = theta.hat-z/sqrt(n)*(theta.hat-1)
  right = theta.hat+z/sqrt(n)*(theta.hat-1)
  ifelse(theta > left & theta < right,1,0)
}
N = 10000
mean(sapply(1:N, function(i) CItest())) # should be around .95
