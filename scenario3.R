# TLSL Scenario 3 --------------------------------------------------------------
simpleSetup::library_install(pkgs)
theme_set(theme_bw())

s3_under_list <- list()
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
    y <- alpha + b1*x1 + rho*x2_df$X2 + epsilon

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

    # Save estimates
    s3_under_list <- results_combiner(s3_under_list, s3_under)
    s3_over_list <- results_combiner(s3_over_list, s3_over)
}

# Find the root mean squared error for DGP variable coefficients
s3_under_list[['rmse']] <- rmse(s3_under_list, c('x1'), 'b1', b1)
s3_over_list[['rmse']] <- rmse(s3_over_list, c('x1', 'x2'),
                            c('b1', 'rho'), c(b1, rho))

# Save simulations -------------------------------------------------------------
save(s3_over_list, s3_under_list, file = 'mc_results/scenario3.rda')

# Plot the results (UNDER) ----------
s3_p_under <- p_plot(s3_under_list, 'lag_wy', 'Scenario 3 (underestimate)')
s3_coef_under <- coef_plot(s3_under_list, 'Scenario 3 (underestimate)',
                           yintercepts = 2)

# Plot the results (OVER) ----------
s3_p_over <- p_plot(s3_over_list, 'lag_wy', 'Scenario 3 (overestimate)')
s3_coef_over <- coef_plot(s3_over_list, 'Scenario 3 (overestimate)',
                          yintercepts = c(0.001, 2))

pdf(file = 'mc_figures/scenario3_plots.pdf', width = 12, height = 12)
gridExtra::grid.arrange(
    s3_p_under, s3_coef_under,
    s3_p_over, s3_coef_over,
    ncol = 2)
dev.off()
