# Find Moran's I for TLSL in scenarios 2 and 3

# Moran's I false discovery rate
fdr_fun <- function(x) (sum(x < 0.05)/length(x))

# Scenario 2 -------------------------------------------------------------------
xfun::pkg_attach2(pkgs)
theme_set(theme_minimal())
set.seed(seed)
mbreaks = c(0, 0.05, 0.1, 0.3, 0.5)


phi_range <- seq(0.1, 0.9, by = 0.2)
s2_moransi <- data.frame()

one_run_phi_moran <- function(n = nsims_less, phi) {
    out <- data.frame()
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
            weight_name = 'wy', mc_cores = 1, morans_i = 'table')
        sw$t <- as.integer(sw$t)
        sw$morans_i_p_value <- as.numeric(sw$morans_i_p_value)
        sw$phi <- phi

    }
    out <- rbind(out, sw)
    return(out)
}

for (p in phi_range) {
    message(p)
    p_char <- as.character(p)
    s2_moransi <- rbind(s2_moransi, one_run_phi_moran(phi = p))
}

s2_morans_fdr <- s2_moransi %>% group_by(phi) %>%
    summarise(morans_i_fdr = fdr_fun(morans_i_p_value))

s2_m_fdr_plot <- ggplot(s2_morans_fdr, aes(phi, morans_i_fdr)) +
    geom_line() +
    scale_y_continuous(limits = c(0, 0.5), breaks = mbreaks) +
    geom_hline(yintercept = 0.05, linetype = 'dotted') +
    ylab('FDR\n') +
    xlab(expression(phi)) +
    ggtitle('Scenario 2')

# Scenario 3 -------------------------------------------------------------------
theta_range <- c(0.00001, 0.0001, 0.001, 0.01, 0.01)
s3_moransi <- data.frame()

one_run_theta_moran <- function(n = nsims_less, theta) {
    out <- data.frame()
    for (u in 1:n) {
        x1 <- sample(x = c(0, 1), size = N, replace = TRUE)
        x2_df <- x2_spatial_builder(tu = t_per_indiv)
        epsilon <- rnorm(N, 0, 1)

        # Generate response
        y <- alpha + b1*x1 + theta*x2_df$X2 + epsilon

        # Create global monadic spatial weight for y
        comb <- data.frame(id = i, t = t, y = y,
                           x1 = x1, x2 = x2_df$X2,
                           location = x2_df$location)
        sw <- spatialWeights::monadic_spatial_weights(
            comb, id_var = 'id', time_var = 't',
            y_var = 'y', location_var = 'location',
            weight_name = 'wy', mc_cores = 1, morans_i = 'table')
        sw$t <- as.integer(sw$t)
        sw$morans_i_p_value <- as.numeric(sw$morans_i_p_value)
        sw$theta <- theta

    }
    out <- rbind(out, sw)
    return(out)
}

for (r in theta_range) {
    message(r)
    p_char <- as.character(r)
    s3_moransi <- rbind(s3_moransi, one_run_theta_moran(theta = r))
}

# Moran's I false discovery rate
fdr_fun <- function(x) (sum(x < 0.05)/length(x))

s3_morans_fdr <- s3_moransi %>% group_by(theta) %>%
    summarise(morans_i_fdr = fdr_fun(morans_i_p_value))

save(s2_morans_fdr, s3_morans_fdr, file = 'mc_results/morans_i_fdr.rda')

s3_m_fdr_plot <- ggplot(s3_morans_fdr, aes(theta, morans_i_fdr)) +
    geom_line() +
    scale_y_continuous(limits = c(0, 0.5), breaks = mbreaks) +
    geom_hline(yintercept = 0.05, linetype = 'dotted') +
    ylab('FDR\n') + xlab(expression(theta[WZ])) +
    ggtitle('Scenario 3')


pdf(file = 'mc_figures/morans_i_figures.pdf', width = 9, height = 4.5)
    ggarrange(s2_m_fdr_plot, s3_m_fdr_plot, ncol = 2)
dev.off()
