set.seed(42)
n <- 5000
x1 <- rnorm(n)
x2 <- runif(n, -1, 1)
beta_true <- c(2, -1, 0.5)
errors <- rt(n, df = 3) # heavy-tailed
y <- as.numeric(cbind(1, x1, x2) %*% beta_true + errors)

df <- data.frame(y = y, x1 = x1, x2 = x2)
result <- GA.QR(
    formula = y ~ x1 + x2,
    data = df,
    tau = 0.5,
    pop_size = 100,
    tournament_size = 10,
    elitism_rate = 0.5,
    max_generations = 1000,
    mutation_prob = 0.5,
    mutation_strength = 0.9,
    convergence_generations = 100,
    seed = 42
)

# result$vcov
# summary(result)
# plot(result, type = "convergence")
# plot(result, type = "residuals")

# result$se
# qrSOTA <- quantreg::rq(y ~ x1 + x2, tau = 0.5, data = df)#$coefficients
# summary(qrSOTA, se = "ker") #


# data.frame(
#     `Vraies val.` = beta_true,
#     `Par quantreg` = qrSOTA,
#     `Par GA` = result$coefficients,
#     `Err. absolue` = beta_true - result$coefficients,
#     `Err. relative` = (beta_true - result$coefficients) / beta_true * 100,
#     check.names = FALSE
# )


# profilage
library(microbenchmark)
aaaa <- microbenchmark(
    result <- GA.QR(
        formula = y ~ x1 + x2,
        data = df,
        tau = 0.5,
        pop_size = 200,
        max_generations = 500,
        mutation_strength = 0.25,
        convergence_generations = 60,
        seed = 42
    ),
    times = 10
)
print(aaaa)
