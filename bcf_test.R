set.seed(1)

p <- 3 # two control variables and one effect moderator
n <- 1000 # number of observations
n_burn <- 2000
n_sim <- 1000

x <- matrix(rnorm(n*(p-1)), nrow=n)
x <- cbind(x, x[,2] + rnorm(n))
weights <- abs(rnorm(n))

# create targeted selection, whereby a practice's likelihood of joining the intervention (pi) 
# is related to their expected outcome (mu)
mu <- -1*(x[,1]>(x[,2])) + 1*(x[,1]<(x[,2])) - 0.1

# generate treatment variable
pi <- pnorm(mu)
z <- rbinom(n,1,pi)

# tau is the true treatment effect. It varies across practices as a function of
# X3 and X2
tau <- 1/(1 + exp(-x[,3])) + x[,2]/10

# generate the expected response using mu, tau and z
y_noiseless <- mu + tau*z

# set the noise level relative to the expected mean function of Y
sigma <- diff(range(mu + tau*pi))/8

# draw the response variable with additive error
y <- y_noiseless + sigma*rnorm(n)/sqrt(weights)