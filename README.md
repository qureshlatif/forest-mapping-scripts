# forest-mapping-scripts
 Scripts for mapping forest restoration benefits for birds
 (Archived scripts not cataloged)

# Data compilation
00-Data_wrangle.R - Compiles data for analysis to produce Data_compiled.RDATA (requires access to Bird Conservancy internal database; not necessary to run if this file has already been compiled).
01-Priority_spp_lists.R - Generates habitat specialization index and categorizes species based on strength of association with ponderosa pine forest. Generates "Spp_list_detected_&_categorized.csv", so not necessary to run if this file has already been generated.

# Analysis
01-Analyze_community_occupancy.R - Implements model fit (not necessary to run if the R object file "mod" has already been generated).

# Results
02-Tabulate_parameters.R - Tabulates summaries of species-specific occupancy model parameter estimates for review.
03-Map_predictions_grid.R - Generates spatial model predictions mapping avian responses to management at the grid cell resolution
03-Map_predictions_point.R - Generates spatial model predictions mapping avian responses to management at the point resolution
03-Plot_community_metrics.R - Plots showing community metrics in relation to model covariates.
03-Plot_effects.R - Plots species-specific covariate effects on occupancy.
03-Plot_QuadCanCov_relations.R - Plots species occupancy in relation to point-level canopy cover where quadratic relations were statistically supported.
03-Plot_species_occupancy.R - Plots individual species occupancy in relation to covariates for review
04-Correlate_catchment_metrics.R - Generates Pearson's correlation coefficients relating avian and non-avian ecosystem metrics (not reported in manuscript)
04-Plot_community_means.R - Plots landscape-wide mean predictions for avian community metrics
04-Summarize_predictions.R - Tabulates summaries of spatial predictions of avian responses to management
05-Plot_gains_and_losses.R - Plots summaries of the extent of statistically supported projected changes in avian community metrics with management.

# Sourced scripts
Data_processing_community_occupancy.R - Sourced file that organizes avian detection data and covariates for analysis
model.jags - JAGS script representing avian community occupancy model