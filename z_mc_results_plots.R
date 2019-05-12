# ------------------------------------------------------------------------------
# Plot MC Experiment Results
# Christopher Gandrud
# MIT License
# ------------------------------------------------------------------------------

library(xfun)
pkg_attach2("dplyr", "ggplot2", "ggpubr", "xtable")
theme_set(theme_minimal())

# Load MC results
rm_lists <- (Filter( function(x) 'list' %in% class( get(x) ), ls() ))
rm(rm_lists)
results_files <- list.files('mc_results')
lapply(sprintf('mc_results/%s', results_files), load, .GlobalEnv)
rm(results_files)
all_results <- Filter( function(x) 'list' %in% class( get(x) ), ls() )
no_range_results <- all_results[!(all_results %in%
                                 c('s2_phi_range_over_list',
                                   's2_phi_range_under_list',
                                   's3_theta_wz_range_over_list',
                                   's3_theta_wz_range_under_list'))]

# False discovery rate formula
fdr_fun <- function(x) (sum(x < 0.05)/length(x))

# Y-axis breas
y_breaks <- c(0, 0.05, 0.25, 0.5, 0.75, 1)

# Plot p-values for TLSL -------------------------------------------------------
pvalues_df <- data.frame()
for (i in no_range_results) {
    message(i)
    ptemp <- extract_element(eval(parse(text=paste(i))), 'pvalue', 'lag_wy')
    ptemp$scenario <- i
    pvalues_df <- rbind(pvalues_df, ptemp)
}

pvalues_df$under_l <- grepl('under', pvalues_df$scenario)

pvalues_df$Type[pvalues_df$under_l] <- "Under-estimated"
pvalues_df$Type[pvalues_df$under_l == FALSE] <- "Over-estimated"

p_labels <- c('Scenario 1 (over)', 'Scenario 1 (under)',
              'Scenario 2 (over)', 'Scenario 2 (under)',
              'Scenario 3 (over)', 'Scenario 3 (under)',
              'Scenario 4 (over)', 'Scenario 4 (under)',
              'Scenario 5 (over)', 'Scenario 5 (under)')
pvalues_df$scenario <- factor(pvalues_df$scenario, labels = p_labels)

no_range_fdr <- pvalues_df %>% group_by(scenario) %>%
    summarise(fdr = fdr_fun(value))

is.even <- function(x) x %% 2 == 0

under <- no_range_fdr[is.even(1:nrow(no_range_fdr)), ]
over <- no_range_fdr[!is.even(1:nrow(no_range_fdr)), ]

no_range_fdr <- data.frame(Scenario = 1:5, Under = under[, 2],
                           Over = over[, 2])
names(no_range_fdr) <- c('Scenario No.', 'FDR (under)', 'FDR (over)')

print(xtable(no_range_fdr), file = 'mc_tables/no_range_fdr.tex',
      floating = FALSE, row.names = FALSE)


# Scenario 2 (omitted autoregressive covariate) ----- --------------------------
# Coefficient Bias
load('mc_results/scenario2_range_phi.rda')
rmse_s2_range_df <- data.frame()
for (i in names(s2_phi_range_under_list)) {
    message(i)
    rmse_temp <- s2_phi_range_under_list[[i]][['rmse']]
    rmse_temp <- cbind(rmse_temp, data.frame(phi = as.numeric(i),
                                             type = 'Scenario 2 under estimate'))
    rmse_s2_range_df <- rbind(rmse_s2_range_df, rmse_temp)
}

for (i in names(s2_phi_range_over_list)) {
    message(i)
    rmse_temp <- s2_phi_range_over_list[[i]][['rmse']]
    rmse_temp <- cbind(rmse_temp, data.frame(phi = as.numeric(i),
                                             type = 'Scenario 2 over estimate'))
    rmse_s2_range_df <- rbind(rmse_s2_range_df, rmse_temp)
}

rmse_s2_range <- ggplot(rmse_s2_range_df, aes(phi, rmse, group = variable, 
                                              linetype = variable)) +
    facet_wrap(~type) +
    geom_line() +
    scale_linetype(name = "", labels = expression(hat(beta)[1], hat(beta)[2])) +
    geom_hline(yintercept = 0, linetype = 'dotted') +
#    scale_y_continuous(limits = c(0, 2.5)) +
    scale_x_continuous(breaks = as.numeric(names(s2_phi_range_under_list))) +
    ylab('') + xlab(expression(phi))

# False Discovery Rate (underestimated)
fdr_s2_range_df1 <- data.frame()
for (i in names(s2_phi_range_under_list)) {
    message(i)
    pvalue <- s2_phi_range_under_list[[i]][['pvalue']]
    pvalue <- pvalue[names(pvalue) == 'lag_wy']
    fdr_temp <- cbind(pvalue, data.frame(phi = as.numeric(i),
                                         type = 'Under estimate'))
    fdr_s2_range_df1 <- rbind(fdr_s2_range_df1, fdr_temp)
}
fdr_s2_range_df1 <- fdr_s2_range_df1 %>% group_by(phi, type) %>%
    summarise(fdr = fdr_fun(pvalue))


# False Discovery Rate (over estimated)
fdr_s2_range_df2 <- data.frame()
for (i in names(s2_phi_range_over_list)) {
    message(i)
    pvalue <- s2_phi_range_over_list[[i]][['pvalue']]
    pvalue <- pvalue[names(pvalue) == 'lag_wy']
    fdr_temp <- cbind(pvalue, data.frame(phi = as.numeric(i),
                                         type = 'Over estimate'))
    fdr_s2_range_df2 <- rbind(fdr_s2_range_df2, fdr_temp)
}
fdr_s2_range_df2 <- fdr_s2_range_df2 %>% group_by(phi, type) %>%
    summarise(fdr = fdr_fun(pvalue))

fdr_s2_range_df <- rbind(fdr_s2_range_df1, fdr_s2_range_df2)

# fdr_s2_range <- ggplot(fdr_s2_range_df, aes(phi, fdr, group = type, 
#                                             linetype = type)) +
#     geom_line(size = 1) +
#     scale_linetype(name = "") +
#     geom_hline(yintercept = 0.05, linetype = 'dotted') +
#     scale_y_continuous(limits = c(0, 1)) +
#     scale_x_continuous(breaks = as.numeric(names(s2_phi_range_under_list))) +
#     ggtitle("Scenario 2") +
#     ylab('') + 
#     xlab(expression(phi))


# Scenario 3 (omitted spatially clustered covariate) ---------------------------
# Coefficient Bias
load('mc_results/scenario3_range_theta_wz.rda')
rmse_s3_range_df <- data.frame()
for (i in names(s3_theta_wz_range_under_list)) {
    message(i)
    rmse_temp <- s3_theta_wz_range_under_list[[i]][['rmse']]
    rmse_temp <- cbind(rmse_temp, data.frame(theta_wz = as.numeric(i),
                                             type = 'Scenario 3 under estimate'))
    rmse_s3_range_df <- rbind(rmse_s3_range_df, rmse_temp)
}

for (i in names(s3_theta_wz_range_over_list)) {
    message(i)
    rmse_temp <- s3_theta_wz_range_over_list[[i]][['rmse']]
    rmse_temp <- cbind(rmse_temp, data.frame(theta_wz = as.numeric(i),
                                             type = 'Scenario 3 over estimate'))
    rmse_s3_range_df <- rbind(rmse_s3_range_df, rmse_temp)
}

rmse_s3_range_df$theta_wz_log <- log(rmse_s3_range_df$theta_wz)

rmse_s3_range <- ggplot(rmse_s3_range_df, aes(theta_wz_log, rmse, group = variable,
                                              linetype = variable)) +
    facet_wrap(~type) +
    geom_line() +
    geom_hline(yintercept = 0, linetype = 'dotted') +
    scale_linetype(name = "", labels = expression(hat(beta)[1], hat(theta)[wz])) +
    scale_x_continuous(breaks = unique(rmse_s3_range_df$theta_wz_log),
                       labels = as.numeric(names(s3_theta_wz_range_under_list))) +
 #   scale_y_continuous(limits = c(0, 2.5)) +
    ylab('') + xlab(expression(paste(theta[wz], ' (log spaced scale)')))


# False Discovery Rate (underestimated)
fdr_s3_range_df1 <- data.frame()
for (i in names(s3_theta_wz_range_under_list)) {
    message(i)
    pvalue <- s3_theta_wz_range_under_list[[i]][['pvalue']]
    pvalue <- pvalue[names(pvalue) == 'lag_wy']
    fdr_temp <- cbind(pvalue, data.frame(theta_wz = as.numeric(i),
                                         type = 'Under estimate'))
    fdr_s3_range_df1 <- rbind(fdr_s3_range_df1, fdr_temp)
}
fdr_s3_range_df1 <- fdr_s3_range_df1 %>% group_by(theta_wz, type) %>%
    summarise(fdr = fdr_fun(pvalue))


# False Discovery Rate (over estimated)
fdr_s3_range_df2 <- data.frame()
for (i in names(s3_theta_wz_range_over_list)) {
    message(i)
    pvalue <- s3_theta_wz_range_over_list[[i]][['pvalue']]
    pvalue <- pvalue[names(pvalue) == 'lag_wy']
    fdr_temp <- cbind(pvalue, data.frame(theta_wz = as.numeric(i),
                                         type = 'Over estimate'))
    fdr_s3_range_df2 <- rbind(fdr_s3_range_df2, fdr_temp)
}
fdr_s3_range_df2 <- fdr_s3_range_df2 %>% group_by(theta_wz, type) %>%
    summarise(fdr = fdr_fun(pvalue))

fdr_s3_range_df <- rbind(fdr_s3_range_df1, fdr_s3_range_df2)
fdr_s3_range_df$theta_wz_log <- log(fdr_s3_range_df$theta_wz)

# fdr_s3_range <- ggplot(fdr_s3_range_df, aes(theta_wz_log, fdr, group = type,
#                                             linetype = type)) +
#     geom_line(size = 1) +
#     scale_linetype(name = "") +
#     geom_hline(yintercept = 0.05, linetype = 'dotted') +
#     scale_x_continuous(breaks = unique(rmse_s3_range_df$theta_wz_log),
#                        labels = as.numeric(names(s3_theta_wz_range_under_list))) +
#     scale_y_continuous(limits = c(0, 1)) +
#     ggtitle("Scenario 3") +
#     ylab('') + 
#     xlab(expression(theta[WZ]))

# Add Moran's I false discovery rate ------
load('mc_results/morans_i_fdr.rda') # created in scen2_scen3_morans_i.R

# Scenario 2
# Clean up to merge
s2_morans_fdr$type <- "Moran's I"
s2_morans_fdr <- s2_morans_fdr %>% 
  dplyr::rename(fdr = morans_i_fdr)
fdr_s2_range_df <- bind_rows(fdr_s2_range_df, s2_morans_fdr)
fdr_s2_sub_df <- subset(fdr_s2_range_df, !(type %in% "Over estimate"))

fdr_s2_range <- ggplot(fdr_s2_sub_df, aes(phi, fdr, group = type, 
                                            linetype = type)) +
    geom_line(size = 1) +
    scale_linetype(name = "") +
    geom_hline(yintercept = 0.05, linetype = 'dotted') +
    scale_y_continuous(limits = c(0, 1),
                       breaks = y_breaks) +
    scale_x_continuous(breaks = as.numeric(names(s2_phi_range_under_list))) +
    ggtitle("Scenario 2") +
    ylab('') + 
    xlab(expression(phi))

# Scenario 3 
# Clean up to merge
s3_morans_fdr$type <- "Moran's I"
s3_morans_fdr <- s3_morans_fdr %>% 
                    dplyr::rename(theta_wz = theta) %>%
                    dplyr::rename(fdr = morans_i_fdr)
s3_morans_fdr$theta_wz_log <- log(s3_morans_fdr$theta_wz)
fdr_s3_range_df <- bind_rows(fdr_s3_range_df, s3_morans_fdr)
fdr_s3_sub_df <- subset(fdr_s3_range_df, !(type %in% "Over estimate"))

# Plot
fdr_s3_range <- ggplot(fdr_s3_sub_df, aes(theta_wz_log, fdr, group = type,
                                            linetype = type)) +
    geom_line(size = 1) +
    scale_linetype(name = "") +
    geom_hline(yintercept = 0.05, linetype = 'dotted') +
    scale_x_continuous(breaks = unique(rmse_s3_range_df$theta_wz_log),
                       labels = as.numeric(names(s3_theta_wz_range_under_list))) +
    scale_y_continuous(limits = c(0, 1), 
                       breaks = y_breaks) +
    ggtitle("Scenario 3") +
    ylab('') + 
    xlab(expression(theta[WZ]))

# Plot FDRs for scenarios 2 and 3 ----
p_fdr_combined <-ggarrange(fdr_s2_range, fdr_s3_range,
                           common.legend = TRUE, legend = "bottom", ncol = 2)
p_fdr_combined <- annotate_figure(p_fdr_combined, left = "TLSL False Discovery Rate")
pdf('mc_figures/fdr_scen_2_3.pdf', width = 10, height = 6)
  p_fdr_combined
dev.off()

# Plot bias for scenarios 2 and 3 ----
p_rmse_combined <-ggarrange(rmse_s2_range, rmse_s3_range, ncol = 2)
p_rmse_combined <- annotate_figure(p_rmse_combined, left = "Bias (RMSE)")

pdf('mc_figures/rmse_scen_2_3.pdf', width = 12, height = 7)
  p_rmse_combined
dev.off()

