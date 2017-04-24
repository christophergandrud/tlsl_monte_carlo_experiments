# TLSL Scenario 4 --------------------------------------------------------------
simpleSetup::library_install(pkgs)
theme_set(theme_bw())

s4_under_list <- list()
set.seed(seed)

for (u in 1:nsims) {
    tu <- t_per_indiv + 1
    comb <- data.frame()
    for (n in 1:n_indiv) {
        epsilon <- rnorm(tu, 0, 1)
        x1 <- sample(x = c(0, 1), size = tu, replace = TRUE)
        x2 <- runif(n = tu, min = 0, max = 1)

        y <- numeric(length(epsilon))
        yinit <- rnorm(1, 0, 1)

        for(l in 1:length(y)){
            if(l==1) y[l] <- alpha + b1*x1[l] + b2*x2[l] + phi * yinit + epsilon[l]
            if(l > 1) y[l] <- alpha + b1*x1[l] + b2*x2[l] + phi * y[l-1] + epsilon[l]
        }
        temp <- data.frame(id = n, t = 1:t_per_indiv, y = y[-1],
                           x1 = x1[-1], x2 = x2[-1],
                           ytm1 = y[-length(y)])
        comb <- rbind(comb, temp)
    }
    location_df <- location_builder_continuous(n_indiv = n_indiv,
                                                t_per_indiv = t_per_indiv)
    comb <- cbind(comb, location = location_df$location)

    # Create global monadic spatial weight
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
    s4_under <- lm(y ~ x1 + x2 + lag_wy, data = sw)

#    s4_over <- lm(y ~ x1 + x2 + ytm1, data = sw)

    # Save estimates
    s4_under_list <- results_combiner(s4_under_list, s4_under)
}

# Plot the results
ps_df_u <- extract_element(s4_under_list, 'pvalue', 'lag_wy')

s4_p <- ggplot(ps_df_u, aes(variable, value)) +
    geom_boxplot() +
    geom_point(alpha = 0.2, position = 'jitter') +
    geom_hline(yintercept = 0.05, linetype = 'dashed', color = 'red', size = 1) +
    geom_hline(yintercept = 0.1, linetype = 'dotted', color = 'red', size = 1) +
    scale_y_continuous(breaks = c(0, 0.05, 0.1, 0.2, 0.5, 1), limits = c(0, 1)) +
    coord_flip() +
    ylab('p-value of temporally-lagged spatial lag') + xlab('') +
    ggtitle('Scenario 4 (mischaracterised)')

# Plot coefficients
coef4_interval <- slim_coefs(s4_under_list)

s4_coef <- ggplot(coef4_interval, aes(variable, qi_median, ymin = qi_min, ymax = qi_max)) +
    geom_pointrange() +
    geom_hline(yintercept = c(2, 3, 4), linetype = 'dotted') +
    geom_hline(yintercept = 0, colour = 'red') +
    ylab('Coefficient Estimate\n') + xlab('\nVariable') +
    ggtitle('Scenario 4 (mischaracterised)')

pdf(file = 'mc_figures/scenario4_plots.pdf', width = 12, height = 6)
    gridExtra::grid.arrange(s4_p, s4_coef, ncol = 2)
dev.off()
