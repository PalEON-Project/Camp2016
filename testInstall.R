# test that R/RStudio/JAGS install worked properly

library(neotoma)

meta <- get_dataset(loc = c(-95, 43, -85, 48), datasettype = "pollen")

print(meta[[1]])

cat("You should have seen some information about Lake Allie.\n")

library(R2jags)

mu <- 1
n <- 10
y <- rnorm(n, mu)
out <- jags(data = list(y=y), model.file = 'test.bug', 
            parameters.to.save = 'mu')

cat("You should have seen some output about 'Compiling model graph' and 'Initializing model' as well as a progress bar.\n")
