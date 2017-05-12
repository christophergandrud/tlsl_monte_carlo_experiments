# TLSL Scenario 2 --------------------------------------------------------------
simpleSetup::library_install(pkgs)
theme_set(theme_bw())

s2_under_list <- list()

s2_under_nolag_list <- list()

s2_over_list <- list()
set.seed(seed)

for (u in 1:nsims) {
    x1 <- sample(x = c(0, 1), size = N, replace = TRUE)
    x2 <- vector()
    for (g in 1:n_indiv) {
        x2 <- c(x2, arima.sim(list(ar = AR), n = t_per_indiv))
    }
    epsilon <- rnorm(N, 0, 1)
    #    epsilon <- vector() # Autoregressive error term, results largely similar
    #    for (g in 1:n_indiv) {
    #        epsilon <- c(epsilon, arima.sim(list(ar = AR),
    #                                        n = t_per_indiv))
    #    }
    location_df <- location_builder_continuous(n_indiv = n_indiv,
                                    t_per_indiv = t_per_indiv)

    # Generate response
    y <- alpha + b1*x1 + b2*x2 + epsilon

    # Create global monadic spatial weight
    comb <- data.frame(id = i, t = t, y = y,
                       x1 = x1, x2 = x2,
                       location = location_df$location)

    sw <- spatialWeights::monadic_spatial_weights(
        comb, id_var = 'id', time_var = 't',
        y_var = 'y', location_var = 'location',
        weight_name = 'wy', mc_cores = num_cores)
    # Lag weight
    sw <- sw %>% arrange(id, t) %>% group_by(id) %>%
        mutate(lag_wy = dplyr::lag(wy, order_by = id))

    sw <- merge(sw, comb)

    # burn in
    sw <- subset(sw, t != 1:burnin)

    # Estimate models
    s2_under <- lm(y ~ x1 + lag_wy, data = sw)

    s2_under_nolag <- lm(y ~ x1 + wy, data = sw)

    s2_over <- lm(y ~ x1 + x2 + lag_wy, data = sw)

    # Save estimates
    s2_under_list <- results_combiner(s2_under_list, s2_under)
    s2_over_list <- results_combiner(s2_over_list, s2_over)

}

# Find the mean squared error for DGP variable coefficients
s2_under_list[['mse']] <- mse(s2_under_list, c('x1'), 'b1', b1)
s2_over_list[['mse']] <- mse(s2_over_list, c('x1', 'x2'),
                            c('b1', 'phi'), c(b1, AR))

# Save simulations -------------------------------------------------------------
save(s2_over_list, s2_under_list, file = 'mc_results/scenario2.rda')


# Plot the results (UNDER) ----------
s2_p_under <- p_plot(s2_under_list, 'lag_wy', 'Scenario 2 (underestimate)')
s2_coef_under <- coef_plot(s2_under_list, 'Scenario 2 (underestimate)')

# Plot the results (OVER) ----------
s2_p_over <- p_plot(s2_over_list, 'lag_wy', 'Scenario 2 (overestimate)')
s2_coef_over <- coef_plot(s2_over_list, 'Scenario 2 (overestimate)')

pdf(file = 'mc_figures/scenario2_plots.pdf', width = 12, height = 12)
gridExtra::grid.arrange(s2_p_under, s2_coef_under, s2_p_over, s2_coef_over,
                        ncol = 2)
dev.off()
