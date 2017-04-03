# TLSL Scenario 1 --------------------------------------------------------------
simpleSetup::library_install(pkgs)
theme_set(theme_bw())

s1_over_list <- list()
s1_under_list <- list()

for (u in 1:nsims) {
    set.seed(u)
    x1 <- sample(x = c(0, 1), size = N, replace = TRUE)
    x2 <- runif(n = N, min = 0, max = 10)
    epsilon <- rnorm(N, 0, 1)
    location <- location_builder(n_indiv = n_indiv, t_per_indiv = t_per_indiv)


    # Generate response
    y <- alpha + b1*x1 + b2*x2 + epsilon

    # Create global monadic spatial weight
    comb <- data.frame(id = i, t = t, y = y,
                       x1 = x1, x2 = x2,
                       location = as.factor(location))
    sw <- spatialWeights::monadic_spatial_weights(
                            comb, id_var = 'id', time_var = 't',
                            y_var = 'y', location_var = 'location',
                            weight_name = 'wy', mc_cores = num_cores)
    # Lag weight
    sw <- sw %>% arrange(id, t) %>% group_by(id) %>%
                mutate(lag_wy = dplyr::lag(wy, order_by = id))

    sw <- merge(sw, comb)

    # Estimate models
    s1_over <- lm(y ~ x1 + x2 + lag_wy, data = sw)
    s1_under <- lm(y ~ x1 + lag_wy, data = sw)

    # Save estimates
    s1_over_list <- results_combiner(s1_over_list, s1_over)
    s1_under_list <- results_combiner(s1_under_list, s1_under)
}

# Plot the results (underestimate)
ps_df_u <- extract_element(s1_under_list, 'pvalue', 'lag_wy')

s1_p_under <- ggplot(ps_df_u, aes(value)) +
    geom_density() +
    scale_x_continuous(breaks = c(0, 0.05, 0.1, 0.2, 0.5, 1),
                       limits = c(0, 1)) +
    geom_vline(xintercept = 0.05, linetype = 'dashed') +
    xlab('\np-value of temporally-lagged spatial lag') +
    ggtitle('Scenario 1 (underestimate)')

# Plot coefficients
coef1_interval_u <- slim_coefs(s1_under_list)

s1_coef_under <- ggplot(coef1_interval_u, aes(variable, qi_median, ymin = qi_min,
                           ymax = qi_max)) +
    geom_pointrange() +
    geom_hline(yintercept = 2, linetype = 'dotted') +
    geom_hline(yintercept = 0, colour = 'red') +
    ylab('Simulated Coefficients') + xlab('\nVariable') +
    ggtitle('Scenario 1 (underestimate)')

# Plot the results (overestimate)
ps_df_o <- extract_element(s1_over_list, 'pvalue', 'lag_wy')

s1_p_over <- ggplot(ps_df_o, aes(value)) +
    geom_density() +
    scale_x_continuous(breaks = c(0, 0.05, 0.1, 0.2, 0.5, 1),
                       limits = c(0, 1)) +
    geom_vline(xintercept = 0.05, linetype = 'dashed') +
    xlab('\np-value of temporally-lagged spatial lag') +
    ggtitle('Scenario 1 (overestimate)')

# Plot coefficients
coef1_interval_o <- slim_coefs(s1_over_list)

s1_coef_over <- ggplot(coef1_interval_o, aes(variable, qi_median, ymin = qi_min,
                                              ymax = qi_max)) +
    geom_pointrange() +
    geom_hline(yintercept = c(2, 3), linetype = 'dotted') +
    geom_hline(yintercept = 0, colour = 'red') +
    ylab('Simulated Coefficients') + xlab('\nVariable') +
    ggtitle('Scenario 1 (overestimate)')

pdf(file = 'mc_figures/scenario1_plots.pdf', width = 12, height = 6)
    gridExtra::grid.arrange(s1_p_under, s1_coef_under, s1_p_over, s1_coef_over,
                            ncol = 2)
dev.off()
