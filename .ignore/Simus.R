# library(quantreg)
library(dplyr)
library(tidyr)
library(kableExtra)
library(knitr)
library(ggplot2)

source("R/quantGA.R")

data.gen <- function(n, taus, beta, seed = 42) {
  set.seed(seed)
  X1 <- runif(n, 0, 10)
  X2 <- rbinom(n, 1, 1/2)
  qshifts <- qnorm(taus)
  Ymat <- sapply(seq_along(taus), function(i) {
    eps <- rnorm(n) - qshifts[i]
    beta[1] + beta[2]*X1 + beta[3] * X2 + eps
  })
  colnames(Ymat) <- paste0("Y_", taus*100, "th")
  data.frame(Ymat, X1 = X1, X2 = X2)
}

fit_tau <- function(df, tau) {
  y_col <- paste0("Y_", tau*100, "th")
  GA.QR(
    formula = as.formula(paste(y_col, "~ X1 + X2")),
    data = df,
    tau = tau,
    pop_size = 100,
    tournament_size = 10,
    elitism_rate = 0.5,
    max_generations = 1000,
    mutation_prob = 0.5,
    mutation_strength = 0.9,
    convergence_generations = 100,
    seed = NULL
  )
}

n      <- 500
taus   <- c(0.1, 0.25, 0.5, 0.75, 0.9)
beta   <- c(1, 2, -3)
n_sims <- 5



allSims <- lapply(1:n_sims, function(i) {
  cat(sprintf("Simulation %d/%d\n", i, n_sims))
  df <- data.gen(n, taus, beta, seed = i)
  
  resList <- list()
  
  # Fit pour chaque tau
  for (tau_val in taus) {
    tryCatch({
      fit_result <- fit_tau(df, tau_val)
      coefNames <- names(fit_result$coefficients)
      
      for (k in coefNames) {
        estVal <- fit_result$coefficients[k]
        estVar <- fit_result$se[k]
        trueBeta <- beta[match(k, c("(Intercept)", "X1", "X2"))]
        
        resList[[length(resList) + 1]] <- data.frame(
          sim = i,
          tau = tau_val,
          coeff = k,
          estimate = estVal,
          estimate_var = estVar,
          true_val = trueBeta,
          stringsAsFactors = FALSE
        )
      }
    }, error = function(e) {
      cat(sprintf("Erreur pour simulation %d, tau=%g : %s\n", 
                  i, tau_val, e$message))
      for (coefJ in c("(Intercept)", "X1", "X2")) {
         resList[[length(resList) + 1]] <- data.frame(
          sim = i,
          tau = tau_val,
          coeff = coefJ,
          estimate = NA,
          estimate_var = NA,
          true_val = beta[match(coefJ, c("(Intercept)", "X1", "X2"))],
          stringsAsFactors = FALSE
        )
      }
    })
  }
  
  bind_rows(resList)
}) %>% bind_rows()

tab <- allSims %>%
  filter(!is.na(estimate)) %>% 
  group_by(tau, coeff) %>%
  summarise(
    `Vraie valeur`            = unique(true_val),
    `Estimation moyenne`      = mean(estimate, na.rm = TRUE),
    `Biais absolu`            = `Estimation moyenne` - `Vraie valeur`,
    `Biais relatif (%)`       = 100*(`Biais absolu` / `Vraie valeur`),
    `Ecart-type asymptotique` = mean(estimate_var, na.rm = TRUE),
    `Ecart-type empirique`    = sd(estimate, na.rm = TRUE),
    `Taux de couverture`      = {
      se_emp <- sd(estimate, na.rm = TRUE)
      lo <- `Estimation moyenne` - 1.96*se_emp
      hi <- `Estimation moyenne` + 1.96*se_emp
      mean(true_val >= lo & true_val <= hi, na.rm = TRUE)
    },
    `Successful simulations` = sum(!is.na(estimate)),
    .groups="drop"
  ) %>%
  arrange(factor(coeff, levels=c("(Intercept)","X1", "X2")), tau)




kable(
  tab,
  digits   = 3,
  col.names = c(
    "tau", "Parameter",
    "True value", "Mean estimate",
    "Absolute bias", "Relative bias (%)", "Asymptotic sd",
    "Empirical sd", "Coverage rate",
    "Successful simulations"
  ),
  caption = "Summary of metrics by coefficient and quantile tau"
) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))






summary_results <- allSims %>%
  filter(!is.na(estimate)) %>%
  mutate(k = coeff)

summary_stats <- data.frame(
  k = c("(Intercept)", "X1", "X2"),
  true_beta = beta
)

ggplot(summary_results, aes(x = factor(tau), y = estimate)) +
  geom_boxplot(aes(fill = k), alpha = 0.7, outlier.shape = NA) +
  geom_jitter(aes(color = k), width = 0.1, alpha = 0.2) +
  geom_hline(data = summary_stats, aes(yintercept = true_beta, color = k), linetype = "dashed", linewidth = 1) +
  facet_wrap(~ k, scales = "free_y") +  labs(title = "Monte Carlo simulations for quantile regression",
       x = "Quantiles",
       y = "Estimated coefficients",
       fill = "Coefficient",
       color = "True value") +
  theme_minimal() +
  scale_color_manual(values = c("red", "blue", "green"))
