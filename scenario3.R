# TLSL Scenario 3 --------------------------------------------------------------
simpleSetup::library_install(pkgs)
theme_set(theme_bw())

s3_under_list <- list()
s3_over_list <- list()
set.seed(seed)

for (u in 1:nsims) {
    x1 <- sample(x = c(0, 1), size = N, replace = TRUE)
    x2_df <- x2_spatial_builder(tu = t_per_indiv)
    epsilon <- rnorm(N, 0, 1)

    # Generate response
    y <- alpha + b1*x1 + theta_wz*x2_df$X2 + epsilon

    # Create global monadic spatial weight for y
    comb <- data.frame(id = i, t = t, y = y,
                       x1 = x1, x2 = x2_df$X2,
                       location = x2_df$location)
    sw <- spatialWeights::monadic_spatial_weights(
        comb, id_var = 'id', time_var = 't',
        y_var = 'y', location_var = 'location',
        weight_name = 'wy', mc_cores = num_cores)
    sw$t <- as.integer(sw$t)
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
                            c('b1', 'theta_wz'), c(b1, theta_wz))

# Save simulations -------------------------------------------------------------
save(s3_over_list, s3_under_list, file = 'mc_results/scenario3.rda')
