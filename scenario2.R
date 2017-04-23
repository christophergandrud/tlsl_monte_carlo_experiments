# TLSL Scenario 2 --------------------------------------------------------------
simpleSetup::library_install(pkgs)
theme_set(theme_bw())

s2_under_list <- list()
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

    # Estimate models
    s2_under <- lm(y ~ x1 + lag_wy, data = sw)
    s2_over <- lm(y ~ x1 + x2 + lag_wy, data = sw)

    # Save estimates
    s2_under_list <- results_combiner(s2_under_list, s2_under)
    s2_over_list <- results_combiner(s2_over_list, s2_over)
}

# Plot the results (UNDER) ----------
ps_df_u <- extract_element(s2_under_list, 'pvalue', 'lag_wy')

s2_p_under <- ggplot(ps_df_u, aes(variable, value)) +
    geom_boxplot() +
    geom_point(alpha = 0.2, position = 'jitter') +
    geom_hline(yintercept = 0.05, linetype = 'dashed', color = 'red', size = 1) +
    geom_hline(yintercept = 0.1, linetype = 'dotted', color = 'red', size = 1) +
    scale_y_continuous(breaks = c(0, 0.05, 0.1, 0.2, 0.5, 1), limits = c(0, 1)) +
    coord_flip() +
    ylab('p-value of temporally-lagged spatial lag') + xlab('') +
    ggtitle('Scenario 2 (underestimate)')

# Plot coefficients
coef2_interval <- slim_coefs(s2_under_list)

s2_coef_under <- ggplot(coef2_interval, aes(variable, qi_median, ymin = qi_min,
                                            ymax = qi_max)) +
    geom_pointrange() +
    geom_hline(yintercept = 2, linetype = 'dotted') +
    geom_hline(yintercept = 0, colour = 'red') +
    ylab('Simulated Coefficients\n') + xlab('\nVariable') +
    ggtitle('Scenario 2 (underestimate)')

# Plot the results (OVER) ----------
ps_df_o <- extract_element(s2_over_list, 'pvalue', 'lag_wy')

s2_p_over <- ggplot(ps_df_o, aes(variable, value)) +
    geom_boxplot() +
    geom_point(alpha = 0.2, position = 'jitter') +
    geom_hline(yintercept = 0.05, linetype = 'dashed', color = 'red') +
    geom_hline(yintercept = 0.1, linetype = 'dotted', color = 'red') +
    scale_y_continuous(breaks = c(0, 0.05, 0.1, 0.2, 0.5, 1), limits = c(0, 1)) +
    coord_flip() +
    ylab('p-value of temporally-lagged spatial lag') + xlab('') +
    ggtitle('Scenario 2 (overestimate)')

# Plot coefficients
coef2_interval <- slim_coefs(s2_over_list)

s2_coef_over <- ggplot(coef2_interval, aes(variable, qi_median, ymin = qi_min,
                                           ymax = qi_max)) +
    geom_pointrange() +
    geom_hline(yintercept = c(2, 3), linetype = 'dotted') +
    geom_hline(yintercept = 0, colour = 'red') +
    ylab('Coefficient Estimate\n') + xlab('\nVariable') +
    ggtitle('Scenario 2 (overestimate)')

pdf(file = 'mc_figures/scenario2_plots.pdf', width = 12, height = 12)
gridExtra::grid.arrange(s2_p_under, s2_coef_under, s2_p_over, s2_coef_over,
                        ncol = 2)
dev.off()
