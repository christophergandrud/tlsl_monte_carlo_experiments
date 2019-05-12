# TLSL Scenario 3 (range of theta_wz) -----------------------------------------------

xfun::pkg_attach2(pkgs)
theme_set(theme_minimal())
set.seed(seed)

theta_wz_range <- c(0.00001, 0.0001, 0.001, 0.01, 0.01)

s3_theta_wz_range_under_list <- list()
s3_theta_wz_range_over_list <- list()

one_run_theta_wz <- function(n = nsims_less, theta_wz, under) {
    out <- list()
    for (u in 1:n) {
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
        if (under)
            model_results <- lm(y ~ x1 + lag_wy, data = sw)
        else
            model_results <- lm(y ~ x1 + x2 + lag_wy, data = sw)

        out <- results_combiner(out, model_results)
    }
    return(out)
}

for (r in theta_wz_range) {
    message(r)
    r_char <- as.character(r)
    s3_theta_wz_range_under_list[[r_char]] <- one_run_theta_wz(theta_wz = r, under = TRUE)
    s3_theta_wz_range_under_list[[r_char]][['rmse']] <- rmse(s3_theta_wz_range_under_list[[r_char]],
                                                        c('x1'), 'b1', b1)
}

for (r in theta_wz_range) {
    message(r)
    r_char <- as.character(r)
    s3_theta_wz_range_over_list[[r_char]] <- one_run_theta_wz(theta_wz = r, under = FALSE)
    s3_theta_wz_range_over_list[[r_char]][['rmse']] <- rmse(s3_theta_wz_range_over_list[[r_char]],
                                       c('x1', 'x2'), c('b1', 'theta_wz'), c(b1, theta_wz))
}


# Save simulations -------------------------------------------------------------
save(s3_theta_wz_range_under_list, s3_theta_wz_range_over_list,
    file = 'mc_results/scenario3_range_theta_wz.rda')
