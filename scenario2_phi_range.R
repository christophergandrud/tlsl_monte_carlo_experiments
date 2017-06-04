# TLSL Scenario 2 (range of phi) -----------------------------------------------

simpleSetup::library_install(pkgs)
theme_set(theme_bw())
set.seed(seed)

phi_range <- seq(0.1, 0.9, by = 0.2)

s2_phi_range_under_list <- list()
s2_phi_range_over_list <- list()

nsims_less <- 100

one_run_phi <- function(n = nsims_less, phi, under) {
    out <- list()
    for (u in 1:n) {

        x1 <- sample(x = c(0, 1), size = N, replace = TRUE)
        x2 <- as.vector(replicate(n_indiv,
                                  arima.sim(list(ar = phi), n = t_per_indiv)))

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
        sw$t <- as.integer(sw$t)
        # Lag weight
        sw <- sw %>% arrange(id, t) %>% group_by(id) %>%
            mutate(lag_wy = dplyr::lag(wy, order_by = id))

        sw <- merge(sw, comb)

        # burn in
        sw <- subset(sw, t != 1:burnin)

        # Estimate models
        if (under)
            model_results <- lm(y ~ x1 + lag_wy, data = sw)
        else
            model_results <- lm(y ~ x1 + x2 + lag_wy, data = sw)

        out <- results_combiner(out, model_results)
    }
    return(out)
}

for (p in phi_range) {
    message(p)
    p_char <- as.character(p)
    s2_phi_range_under_list[[p_char]] <- one_run_phi(phi = p, under = TRUE)
    s2_phi_range_under_list[[p_char]][['rmse']] <- rmse(s2_phi_range_under_list[[p_char]],
                                                        c('x1'), 'b1', b1)
}

for (p in phi_range) {
    message(p)
    p_char <- as.character(p)
    s2_phi_range_over_list[[p_char]] <- one_run_phi(phi = p, under = FALSE)
    s2_phi_range_over_list[[p_char]][['rmse']] <- rmse(s2_phi_range_over_list[[p_char]],
                                       c('x1', 'x2'), c('b1', 'b2'), c(b1, p))
}


# Save simulations -------------------------------------------------------------
save(s2_phi_range_under_list, s2_phi_range_over_list,
    file = 'mc_results/scenario2_range_phi.rda')
