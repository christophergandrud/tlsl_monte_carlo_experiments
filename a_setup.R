# ------------------------------------------------------------------------------
# Set up MC Experiments
# Christopher Gandrud
# MIT License
# ------------------------------------------------------------------------------

setwd('~/Dropbox/dynsimRwriteUp/tlsl/tlsl_monte_carlo_experiments/')

library(simpleSetup)
pkgs <- c('spatialWeights', 'car', 'dplyr', 'Matrix', 'ggplot2', 'gridExtra',
          'ape')
simpleSetup::library_install(pkgs)

theme_set(theme_bw())

# Set seed
seed <- 1234

# Number of cores
num_cores <- 7

#### TRUE VALUES
# Number of simulations
nsims = 1000

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
phi = 4

# AR
AR = 0.6

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
    G <- vector()
    for (k in 1:tu) {
        g <- runif(n = obs_per_time, min = 0, max = 1)
        X2 <- c(G, colSums(W * g) + rnorm(obs_per_time, 0, 1))
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
    if (!missing(var)) sub <- subset(sub, vars == var)
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
