# TLSL Scenario 3 --------------------------------------------------------------
simpleSetup::library_install(pkgs)
theme_set(theme_bw())

s3_under_list <- list()
s3_over_list <- list()
s3_under_location_list <- list()

for (u in 1:nsims) {
    set.seed(u)

    x1 <- sample(x = c(0, 1), size = N, replace = TRUE)

    location <- sample(regions, n_indiv, replace = TRUE)
    W <- as.matrix(dist(location, method = 'binary'))
    diag(W) <- 0
    x2_df <- data.frame()
    G <- vector()
    for (k in 1:t_per_indiv) {
        g <- runif(n = obs_per_time, min = 0, max = 10)
        X2 <- c(G, (W %*% g) + rnorm(obs_per_time, 0, 1))
        temp <- data.frame(i = 1:n_indiv, t = k, X2 = X2, location = location)
        x2_df <- rbind(x2_df, temp)
    }
    x2_df <- x2_df[order(x2_df$i, x2_df$t), ]

    epsilon <- rnorm(N, 0, 1)

    # Generate response
    y <- alpha + b1*x1 + b2*x2_df$X2 + epsilon

    # Create global monadic spatial weight
    comb <- data.frame(id = i, t = t, y = y,
                       x1 = x1, x2 = x2_df$X2,
                       location = as.factor(x2_df$location))
    sw <- spatialWeights::monadic_spatial_weights(
                            comb, id_var = 'id', time_var = 't',
                            y_var = 'y', location_var = 'location',
                            weight_name = 'wy',
                            mc_cores = num_cores)
    # Lag weight
    sw <- sw %>% arrange(id, t) %>% group_by(id) %>%
                mutate(lag_wy = dplyr::lag(wy, order_by = id))

   sw <- merge(sw, comb)

    # Estimate models
    s3_under <- lm(y ~ x1 + lag_wy, data = sw)
    s3_over <- lm(y ~ x1 + x2 + lag_wy, data = sw)

    sw$location <- as.factor(sw$location)
    s3_location <- lm(y ~ x1 + location, data = sw)

    # Save estimates
    s3_under_list <- results_combiner(s3_under_list, s3_under)
    s3_over_list <- results_combiner(s3_over_list, s3_over)
    s3_under_location_list <- results_combiner(s3_under_location_list, s3_location)
}

# Plot the results lag p-value (UNDER) -----------------------------------------
ps_df_u <- extract_element(s3_under_list, 'pvalue', 'lag_wy')

s3_p_under <- ggplot(ps_df_u, aes(variable, value)) +
    geom_boxplot() +
#    scale_x_continuous(breaks = c(0, 0.05, 0.1, 0.2, 0.5, 1), limits = c(0, 1)) +
#    geom_vline(xintercept = 0.05, linetype = 'dashed') +
#    xlab('\np-value of temporally-lagged spatial lag') +
    ylab('p-value of temporally-lagged spatial lag') + xlab('') +
    ggtitle('Scenario 3 (underestimate)')

# Plot coefficients
coef3_interval_u <- slim_coefs(s3_under_list)

s3_coef_under <- ggplot(coef3_interval_u, aes(variable, qi_median,
                                              ymin = qi_min, ymax = qi_max)) +
    geom_pointrange() +
    geom_hline(yintercept = c(2), linetype = 'dotted') +
    geom_hline(yintercept = 0, colour = 'red') +
    ylab('Coefficient Estimate\n') + xlab('\nVariable') +
    ggtitle('Scenario 3 (underestimate)')

# Plot the results lag p-value (UNDER, location dummy rather than lag) ---------
ps_df_ul <- extract_element(s3_under_location_list, 'pvalue',
                           c('location1', 'location2', 'location3', 'location4'))

s3_p_underl <- ggplot(ps_df_ul, aes(value, group = variable,
                                    colour = variable)) +
    geom_density() +
    scale_x_continuous(breaks = c(0, 0.05, 0.1, 0.2, 0.5, 1), limits = c(0, 1)) +
    geom_vline(xintercept = 0.05, linetype = 'dashed') +
    xlab('\np-value of temporally-lagged spatial lag') +
    ggtitle('Scenario 3 (underestimate, location dummies)')

# Plot coefficients
coef3_interval_ul <- slim_coefs(s3_under_location_list)

s3_coef_underl <- ggplot(coef3_interval_ul, aes(variable, qi_median,
                                              ymin = qi_min, ymax = qi_max)) +
    geom_pointrange() +
    geom_hline(yintercept = c(2), linetype = 'dotted') +
    geom_hline(yintercept = 0, colour = 'red') +
    ylab('Coefficient Estimate\n') + xlab('\nVariable') +
    ggtitle('Scenario 3 (underestimate, location dummies)')

# Plot the results lag p-value (OVER) ------------------------------------------
ps_df_o <- extract_element(s3_over_list, 'pvalue', 'lag_wy')

s3_p_over <- ggplot(ps_df_o, aes(value)) +
    geom_density() +
    scale_x_continuous(breaks = c(0, 0.05, 0.1, 0.2, 0.5, 1), limits = c(0, 1)) +
    geom_vline(xintercept = 0.05, linetype = 'dashed') +
    xlab('\np-value of temporally-lagged spatial lag') +
    ggtitle('Scenario 3 (overestimate)')

# Plot coefficients
coef3_interval_o <- slim_coefs(s3_over_list)

s3_coef_over <- ggplot(coef3_interval_o, aes(variable, qi_median, ymin = qi_min,
                                             ymax = qi_max)) +
    geom_pointrange() +
    geom_hline(yintercept = c(2, 3), linetype = 'dotted') +
    geom_hline(yintercept = 0, colour = 'red') +
    ylab('Coefficient Estimate\n') + xlab('\nVariable') +
    ggtitle('Scenario 3 (overestimate)')

pdf(file = 'mc_figures/scenario3_plots.pdf', width = 12, height = 18)
    gridExtra::grid.arrange(
                            s3_p_under,
                            s3_coef_under,
                            s3_p_underl, s3_coef_underl,
                            s3_p_over, s3_coef_over,
                            ncol = 2)
dev.off()
