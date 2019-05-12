# TLSL Scenario 1 --------------------------------------------------------------
xfun::pkg_attach2(pkgs)
theme_set(theme_minimal())

s1_over_list <- list()
s1_under_list <- list()
set.seed(seed)

for (u in 1:nsims) {
    x1 <- sample(x = c(0, 1), size = N, replace = TRUE)
    x2 <- rnorm(n = N, 0, 1)
    epsilon <- rnorm(N, 0, 1)
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
    s1_over <- lm(y ~ x1 + x2 + lag_wy, data = sw)
    s1_under <- lm(y ~ x1 + lag_wy, data = sw)

    # Save estimates
    s1_over_list <- results_combiner(s1_over_list, s1_over)
    s1_under_list <- results_combiner(s1_under_list, s1_under)
}

# Find the root mean squared error for DGP variable coefficients
s1_under_list[['rmse']] <- rmse(s1_under_list, c('x1'), 'b1', b1)
s1_over_list[['rmse']] <- rmse(s1_over_list, c('x1', 'x2'),
                             c('b1', 'b2'), c(b1, b2))

# Save simulations -------------------------------------------------------------
save(s1_over_list, s1_under_list, file = 'mc_results/scenario1.rda')
