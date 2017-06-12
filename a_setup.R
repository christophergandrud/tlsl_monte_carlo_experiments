# ------------------------------------------------------------------------------
# Set up MC Experiments
# Christopher Gandrud
# MIT License
# ------------------------------------------------------------------------------

# Install and load required packages
library(simpleSetup)
pkgs <- c('devtools', 'spatialWeights', 'car', 'dplyr', 'Matrix', 'ggplot2',
          'gridExtra', 'ape', 'tibble', 'xtable')

if (!('spatialWeights' %in% installed.packages()[,1]))
    devtools::install_github('christophergandrud/spatialWeights')

simpleSetup::library_install(pkgs[-1])

# Set working directory
dirs <- c('~/Dropbox/dynsimRwriteUp/tlsl/tlsl_monte_carlo_experiments/',
          '/nfs/home/C/cgandrud/git_repositories/tlsl_monte_carlo_experiments')
set_valid_wd(dirs)

# Use BW plotting theme
theme_set(theme_bw())

# Set seed
seed <- 1234

# Number of cores
num_cores <- 7

# Number of simulations
nsims = 100

nsims_less <- 100

# Burn in
burnin <- 2

# Number of time points per unit
t_per_indiv = 100 + burnin

# Number of units
n_indiv <- 100

# Number of simulated observations
N = n_indiv * t_per_indiv

# Time points
t <- rep(x = 1:t_per_indiv, n_indiv)

# Individuals
i <- rep(1:n_indiv, times = 1, each = t_per_indiv)

# Observations per time point
obs_per_time <- N/t_per_indiv

# Parameters
alpha = 1
b1 = 2
b2 = 3
theta_wz = 0.001
AR <- phi <- 0.6

#### FUNCTIONS
# Function to randomly create a location variable # NOT USED
location_builder_discrete <- function(n_indiv = n_indiv,
                             t_per_indiv = t_per_indiv,
                             n_regions = nregions) {
    location <- vector()
    for (j in 1:n_indiv) {
        temp <- rep(sample(0:n_regions, 1), t_per_indiv)
        location <- c(location, temp)
    }
    return(location)
}

location_builder_continuous <- function(n_indiv = n_indiv,
                                        t_per_indiv = t_per_indiv) {
    library(dplyr)
    location_1 <- rnorm(n = n_indiv, mean = 0, sd = 1)
    if (!is.null(t_per_indiv)) {
        location <- rep(location_1, times = t_per_indiv)
        id <- rep(1:n_indiv, times = t_per_indiv)
    }
    else {
        location <- location_1
        id <- 1:n_indiv
    }
    location_df <- data.frame(id = id, location = location) %>% arrange(id)
    return(location_df)
}

# X2 spatially clustered builder
x2_spatial_builder <- function(tu) {
    # Continuous location variable
    location <- rnorm(n = n_indiv, mean = 0, sd = 1)
    W <- as.matrix(dist(location))
    diag(W) <- 0
    x2_df <- data.frame()
    Z <- vector()
    for (k in 1:tu) {
        g <- runif(n = obs_per_time, min = 0, max = 1)
        X2 <- c(Z, (colSums(W * g) + rnorm(obs_per_time, 0, 1)))
        temp <- data.frame(i = 1:n_indiv, t = k, X2 = X2, location = location)
        x2_df <- rbind(x2_df, temp)
    }
    x2_df <- x2_df[order(x2_df$i, x2_df$t), ]
    return(x2_df)
}

# Function to extract standard errors from fitted model objects
se <- function(x) summary(x)[["coefficients"]][, 2]

# Function to extract p-values from fitted model objects
pv <- function(x) {
    full <- summary(x)[["coefficients"]]
    p <- full[, 'Pr(>|t|)']
    return(p)
}

# Combine resluts from each simulation
results_combiner <- function(l, m) {
    l[['coefs']] <- c(l[['coefs']], coef(m))
    l[['se']] <- c(l[['se']], se(m))
    l[['pvalue']] <- c(l[['pvalue']], pv(m))
    l[['vif']] <- c(l[['vif']], vif(m))

    # Need Moran's I
    # Need QI
    return(l)
}

extract_element <- function(results, type, var) {
    all <- results[[type]]
    sub <- cbind(as.data.frame(all), vars = attributes(all)$names)
    if (!missing(var))
        sub <- subset(sub, vars == var)
    names(sub) <- c('value', 'variable')
    return(sub)
}

qi_slimmer <- function(df, scenario_var, qi_var){
    qi_ <- scenario_ <- NULL

    names(df)[names(df) == qi_var] <- 'qi_'
    names(df)[names(df) == scenario_var] <- 'scenario_'

    if (!(names(df)[[ncol(df)]] == 'scenario_'))
        df$scenario_ <- interaction(df[, 1:(ncol(df)-1)], drop = TRUE)

    df_out <- df %>% group_by(scenario_) %>%
        summarise(qi_min = min(qi_),
                  qi_median = median(qi_),
                  qi_max = max(qi_)
        ) %>%
        data.frame

    names(df_out)[names(df_out) == 'qi_'] <- qi_var
    names(df_out)[names(df_out) == 'scenario_'] <- scenario_var
    return(df_out)
}

slim_coefs <- function(results) {
    coefs_raw <- extract_element(results, 'coefs')
    coef_interval <- qi_slimmer(coefs_raw, scenario_var = 'variable',
                                qi_var = 'value')
    coef_interval <- subset(coef_interval, variable != '(Intercept)')
    return(coef_interval)
}

# Root mean squared error
rmse <- function(results, vars, param_labels, p) {
    coefs <- extract_element(results, type = 'coefs')
    rmse_fun__ <- function(phat, p) sqrt(mean((phat - p)^2))
    rmse_df <- data.frame()
    for (i in vars) {
        position <- grep(i, vars)
        temp <- coefs[coefs$variable == i, ]
        rmse_value <- rmse_fun__(temp$value, p = p[position])
        rmse_df <- rbind(rmse_df, data.frame(variable = param_labels[position],
                                           rmse = rmse_value))
    }
    return(rmse_df)
}

# Plotters
p_plot <- function(results, var, title) {
    extracted <- extract_element(results, 'pvalue', var)

    p <- ggplot(extracted, aes(variable, value)) +
        geom_boxplot() +
        geom_point(alpha = 0.2, position = 'jitter') +
        geom_hline(yintercept = 0.05, linetype = 'dashed', color = 'red', size = 1) +
        geom_hline(yintercept = 0.1, linetype = 'dotted', color = 'red', size = 1) +
        scale_y_continuous(breaks = c(0, 0.05, 0.1, 0.2, 0.5, 1), limits = c(0, 1)) +
        coord_flip() +
        ylab('p-value of temporally-lagged spatial lag') + xlab('') +
        ggtitle(title)
    return(p)
}

coef_plot <- function(results, title, yintercepts) {
    extracted <- slim_coefs(results)

    p <- ggplot(extracted, aes(variable, qi_median, ymin = qi_min,
                                               ymax = qi_max)) +
        geom_pointrange() +
        geom_hline(yintercept = 0, colour = 'red') +
        ylab('Coefficient Estimate\n') + xlab('\nVariable') +
        ggtitle(title)
    if (missing(yintercepts))
        yi <- c(2, 3)
    else
        yi = yintercepts
    p <- p + geom_hline(yintercept = yi, linetype = 'dotted')

    return(p)
}
