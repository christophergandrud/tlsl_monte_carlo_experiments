# ------------------------------------------------------------------------------
# Plot MC Experiment Results
# Christopher Gandrud
# MIT License
# ------------------------------------------------------------------------------

library(ggplot2)
theme_set(theme_bw())

# Load MC results
rm_lists <- (Filter( function(x) 'list' %in% class( get(x) ), ls() ))
rm(rm_lists)
results_files <- list.files('mc_results')
lapply(sprintf('mc_results/%s', results_files), load, .GlobalEnv)
rm(results_files)
all_results <- Filter( function(x) 'list' %in% class( get(x) ), ls() )

# Plot p-values for TLSL -------------------------------------------------------
pvalues_df <- data.frame()
for (i in all_results) {
    message(i)
    ptemp <- extract_element(eval(parse(text=paste(i))), 'pvalue', 'lag_wy')
    ptemp$scenario <- i
    pvalues_df <- rbind(pvalues_df, ptemp)
}

pvalues_df$under_l <- grepl('under', pvalues_df$scenario)

pvalues_df$Type[pvalues_df$under_l] <- "Under Estimated"
pvalues_df$Type[pvalues_df$under_l == FALSE] <- "Over Estimated"

p_labels <- c('Scenario 1 (over)', 'Scenario 1 (under)',
              'Scenario 2 (over)', 'Scenario 2 (under)',
              'Scenario 3 (over)', 'Scenario 3 (under)',
              'Scenario 4 (over)', 'Scenario 4 (under)',
              'Scenario 5 (over)', 'Scenario 5 (under)')
pvalues_df$scenario <- factor(pvalues_df$scenario, labels = p_labels)

ggplot(pvalues_df, aes(scenario, value, group = scenario,
                       color = Type)) +
    geom_boxplot() +
    geom_point(alpha = 0.2, position = 'jitter') +
    geom_hline(yintercept = 0.05, linetype = 'dashed', color = 'red', size = 1) +
    geom_hline(yintercept = 0.1, linetype = 'dotted', color = 'red', size = 1) +
    scale_y_continuous(breaks = c(0, 0.05, 0.1, 0.2, 0.5, 1), limits = c(0, 1)) +
    scale_color_manual(values =  c('#bdbdbd', '#636363')) +
    coord_flip() +
    ylab('p-values for temporally-lagged spatial lag') + xlab('')

ggsave('mc_figures/mc_lagwy_pvalues.pdf', height = 12, width = 10)


# Coefficient Bias ----- -------------------------------------------------------
mse_df <- data.frame()
for (i in all_results) {
    message(i)
    mtemp <- eval(parse(text=paste(i)))[['mse']]
    mtemp$scenario <- i
    mse_df <- rbind(mse_df, mtemp)
}

mse_df$under <- grepl('under', mse_df$scenario)

ggplot(mse_df, aes(scenario, mse, group = variable, colour = variable,
                   shape = variable)) +
    facet_grid(~under) +
    geom_point()
