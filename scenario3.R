# TLSL Scenario 3 --------------------------------------------------------------
simpleSetup::library_install(pkgs)
theme_set(theme_bw())

s3_under_list <- list()

s3_under_loc_list <- list()

s3_over_list <- list()
set.seed(seed)

for (u in 1:nsims) {
    x1 <- sample(x = c(0, 1), size = N, replace = TRUE)

    x2_df <- x2_spatial_builder(tu = t_per_indiv)
#    message("Original X2 Moran's I per year")
#    for (l in 1:t_per_indiv) {
#        sub <- subset(x2_df, t == l)$X2
#        print(sprintf('%s: %s', l, format.pval(Moran.I(sub, W)$p.value)))
#    }
    epsilon <- rnorm(N, 0, 1)

    # Generate response
    y <- alpha + b1*x1 + b2*x2_df$X2 + epsilon

    # Create global monadic spatial weight for y
    comb <- data.frame(id = i, t = t, y = y,
                       x1 = x1, x2 = x2_df$X2,
                       location = x2_df$location)
    sw <- spatialWeights::monadic_spatial_weights(
        comb, id_var = 'id', time_var = 't',
        y_var = 'y', location_var = 'location',
        weight_name = 'wy', mc_cores = num_cores)
    # Lag weight
    sw <- sw %>% arrange(id, t) %>% group_by(id) %>%
        mutate(lag_wy = dplyr::lag(wy, order_by = id))

    sw <- merge(sw, comb)

    # burnin
    sw <- subset(sw, t != 1:burnin) %>% arrange(id, t)

    # Estimate models
    s3_under <- lm(y ~ x1 + lag_wy, data = sw)
    s3_over <- lm(y ~ x1 + x2 + lag_wy, data = sw)

    s3_under_loc <- lm(y ~ x1 + location, data = sw)

    # Save estimates
    s3_under_list <- results_combiner(s3_under_list, s3_under)
    s3_over_list <- results_combiner(s3_over_list, s3_over)

    s3_under_loc_list <- results_combiner(s3_under_loc_list, s3_under_loc)
}

# Plot the results lag p-value (UNDER) -----------------------------------------
ps_df_u <- extract_element(s3_under_list, 'pvalue', 'lag_wy')

s3_p_under <- ggplot(ps_df_u, aes(variable, value)) +
    geom_boxplot() +
    geom_point(alpha = 0.2, position = 'jitter') +
    geom_hline(yintercept = 0.05, linetype = 'dashed', color = 'red', size = 1) +
    geom_hline(yintercept = 0.1, linetype = 'dotted', color = 'red', size = 1) +
    scale_y_continuous(breaks = c(0, 0.05, 0.1, 0.2, 0.5, 1), limits = c(0, 1)) +
    coord_flip() +
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
    ggtitle('Scenario 3 (continuous, underestimate)')

# Plot the results lag p-value (location, UNDER) -----------------------------------------
ps_df_uloc <- extract_element(s3_under_loc_list, 'pvalue', 'location')

s3_p_underloc <- ggplot(ps_df_uloc, aes(variable, value)) +
    geom_boxplot() +
    geom_point(alpha = 0.2, position = 'jitter') +
    geom_hline(yintercept = 0.05, linetype = 'dashed', color = 'red', size = 1) +
    geom_hline(yintercept = 0.1, linetype = 'dotted', color = 'red', size = 1) +
    scale_y_continuous(breaks = c(0, 0.05, 0.1, 0.2, 0.5, 1), limits = c(0, 1)) +
    coord_flip() +
    ylab('p-value of temporally-lagged spatial lag') + xlab('') +
    ggtitle('Scenario 3 (raw location, underestimate)')

# Plot coefficients
coef3_interval_uloc <- slim_coefs(s3_under_loc_list)

s3_coef_underloc <- ggplot(coef3_interval_uloc, aes(variable, qi_median,
                                              ymin = qi_min, ymax = qi_max)) +
    geom_pointrange() +
    geom_hline(yintercept = c(2), linetype = 'dotted') +
    geom_hline(yintercept = 0, colour = 'red') +
    ylab('Coefficient Estimate\n') + xlab('\nVariable') +
    ggtitle('Scenario 3 (continuous, underestimate)')

# Plot the results lag p-value (OVER) ------------------------------------------
ps_df_o <- extract_element(s3_over_list, 'pvalue', 'lag_wy')

s3_p_over <- ggplot(ps_df_o, aes(variable, value)) +
    geom_boxplot() +
    geom_point(alpha = 0.2, position = 'jitter') +
    geom_hline(yintercept = 0.05, linetype = 'dashed', color = 'red', size = 1) +
    geom_hline(yintercept = 0.1, linetype = 'dotted', color = 'red', size = 1) +
    scale_y_continuous(breaks = c(0, 0.05, 0.1, 0.2, 0.5, 1), limits = c(0, 1)) +
    coord_flip() +
    ylab('p-value of temporally-lagged spatial lag') + xlab('') +
    ggtitle('Scenario 3 (overestimate)')

# Plot coefficients
coef3_interval_o <- slim_coefs(s3_over_list)

s3_coef_over <- ggplot(coef3_interval_o, aes(variable, qi_median, ymin = qi_min,
                                             ymax = qi_max)) +
    geom_pointrange() +
    geom_hline(yintercept = c(2, 3), linetype = 'dotted') +
    geom_hline(yintercept = 0, colour = 'red') +
    ylab('Coefficient Estimate\n') + xlab('\nVariable') +
    ggtitle('Scenario 3 (continuous,overestimate)')

pdf(file = 'mc_figures/scenario3_plots.pdf', width = 18, height = 18)
gridExtra::grid.arrange(
    s3_p_under, s3_coef_under,
    s3_p_underloc, s3_coef_underloc,
    s3_p_over, s3_coef_over,
    ncol = 2)
dev.off()
