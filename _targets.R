source("R/packages.R")
source("R/functions.R")
source("R/traitDependent_functions.R")
options(timeout = max(10000, getOption("timeout")))
options(clustermq.scheduler="multiprocess")

list(
  # Download and load files
  tar_target(tree, download_tree(), format="file"),
  
  # Input files which cannot be downloaded automatically
  tar_target(climate_map, "data/Beck_KG_V1/Beck_KG_V1_present_0p0083.tif", format="file"),
  tar_target(confidence_map, "data/Beck_KG_V1/Beck_KG_V1_present_conf_0p0083.tif", format="file"),
  tar_target(redlist_original, "data/MAMMALS/MAMMALS.shp", format="file"),

  tar_target(redlist_preprocessed, preprocessing_redlist(redlist_original)), # Preprocessing
  tar_target(redlist_climate_original, merge_redlist_with_climate(redlist_preprocessed, climate_map, confidence_map)), # Merge IUCN Red List data with Köppen-Geiger climate classification (This step will take about a day)
  tar_target(redlist_climate, summarize_area_by_species(redlist_climate_original)), # Summarize area of climate class by species names
  tar_target(shannon_index, calculate_shannon_index(redlist_climate)), # Calculate the Shannon index
  tar_target(threshold, calculate_threshold(shannon_index)), # Set the threshold for the classification of biome specialists and generalists
  tar_target(character, annotate_specialist_generalist(shannon_index, threshold)), # Classify species into biome specialists and generalists
  tar_target(tree_binomial, extract_species_from_tree(tree)), # Extract binomial names from the phylogenetic tree

  tar_target(sampling_step, 100000),
  tar_target(removed_step, 499),
  tar_target(sampling_interval, 500),

  # "Mammalia", "Rodentia", "Primates", "Cetartiodactyla", "Chiroptera", "Eulipotyphla", "Marsupialia"
  tar_target(taxa_woCarnivora, c("Mammalia", "Rodentia", "Primates", "Cetartiodactyla", "Chiroptera", "Eulipotyphla", "Marsupialia")),
  tar_target(redlist_binomial_woCarnivora, extract_species_from_redlist(character, taxa_woCarnivora), pattern=map(taxa_woCarnivora)),
  tar_target(redlist_binomial_in_tree_woCarnivora, map_redlist_to_tree(redlist_binomial_woCarnivora, tree_binomial), pattern=map(redlist_binomial_woCarnivora)),
  tar_target(pruned_tree_woCarnivora, prune_tree(tree, redlist_binomial_in_tree_woCarnivora), pattern=map(redlist_binomial_in_tree_woCarnivora)),
  tar_target(sampling_fraction_woCarnivora, calculate_sampling_fraction(character, redlist_binomial_woCarnivora, redlist_binomial_in_tree_woCarnivora), pattern=map(redlist_binomial_woCarnivora, redlist_binomial_in_tree_woCarnivora)),
  tar_target(bisse_result_woCarnivora, run_bisse(character, pruned_tree_woCarnivora, redlist_binomial_in_tree_woCarnivora, sampling_fraction_woCarnivora, iteration=sampling_step), pattern=map(pruned_tree_woCarnivora, redlist_binomial_in_tree_woCarnivora, sampling_fraction_woCarnivora)),
  tar_target(shannon_index_tree_figure_woCarnivora, plot_shannon_index(shannon_index, threshold, redlist_binomial_woCarnivora, redlist_binomial_in_tree_woCarnivora, taxa_woCarnivora), pattern=map(redlist_binomial_woCarnivora, redlist_binomial_in_tree_woCarnivora, taxa_woCarnivora)),
  tar_target(shannon_index_all_figure_woCarnivora, plot_shannon_index_all(shannon_index, threshold, redlist_binomial_woCarnivora, taxa_woCarnivora), pattern=map(redlist_binomial_woCarnivora, taxa_woCarnivora)),
  tar_target(likelihood_figure_woCarnivora, plot_likelihood(bisse_result_woCarnivora, taxa_woCarnivora, removed_step=removed_step), pattern=map(bisse_result_woCarnivora, taxa_woCarnivora)),
  tar_target(autocorrelation_figure_woCarnivora, plot_autocorrelation(bisse_result_woCarnivora, taxa_woCarnivora, removed_step=removed_step, max_interval=sampling_interval*2, threshold=sampling_interval), pattern=map(bisse_result_woCarnivora, taxa_woCarnivora)),
  tar_target(bisse_figure_woCarnivora, plot_bisse(bisse_result_woCarnivora, taxa_woCarnivora, burn_in=removed_step, interval=sampling_interval), pattern=map(bisse_result_woCarnivora, taxa_woCarnivora)),
  tar_target(pruned_tree_figure_woCarnivora, plot_tree(pruned_tree_woCarnivora, character, redlist_binomial_in_tree_woCarnivora, taxa_woCarnivora), pattern=map(pruned_tree_woCarnivora, redlist_binomial_in_tree_woCarnivora, taxa_woCarnivora)),

  # "Carnivora" (I skipped ML estimation in Carnivora to avoid warnings about the likelihood calculation failure)
  tar_target(taxa_Carnivora, "Carnivora"),
  tar_target(redlist_binomial_Carnivora, extract_species_from_redlist(character, taxa_Carnivora)),
  tar_target(redlist_binomial_in_tree_Carnivora, map_redlist_to_tree(redlist_binomial_Carnivora, tree_binomial)),
  tar_target(pruned_tree_Carnivora, prune_tree(tree, redlist_binomial_in_tree_Carnivora)),
  tar_target(sampling_fraction_Carnivora, calculate_sampling_fraction(character, redlist_binomial_Carnivora, redlist_binomial_in_tree_Carnivora)),
  tar_target(bisse_result_Carnivora, run_bisse_skip_ml(character, pruned_tree_Carnivora, redlist_binomial_in_tree_Carnivora, sampling_fraction_Carnivora, iteration=sampling_step)),
  tar_target(shannon_index_tree_figure_Carnivora, plot_shannon_index(shannon_index, threshold, redlist_binomial_Carnivora, redlist_binomial_in_tree_Carnivora, taxa_Carnivora)),
  tar_target(shannon_index_all_figure_Carnivora, plot_shannon_index_all(shannon_index, threshold, redlist_binomial_Carnivora, taxa_Carnivora)),
  tar_target(likelihood_figure_Carnivora, plot_likelihood(bisse_result_Carnivora, taxa_Carnivora, removed_step=removed_step)),
  tar_target(autocorrelation_figure_Carnivora, plot_autocorrelation(bisse_result_Carnivora, taxa_Carnivora, removed_step=removed_step, max_interval=sampling_interval*2, threshold=sampling_interval)),
  tar_target(bisse_figure_Carnivora, plot_bisse(bisse_result_Carnivora, taxa_Carnivora, burn_in=removed_step, interval=sampling_interval)),
  tar_target(pruned_tree_figure_Carnivora, plot_tree(pruned_tree_Carnivora, character, redlist_binomial_in_tree_Carnivora, taxa_Carnivora)),
  
  # Visualize results from these mammalian lineages at the same time
  # Figure S3
  tar_target(shannon_index_all_summary_woCarnivora, summarize_shannon_index_all(shannon_index, redlist_binomial_woCarnivora, taxa_woCarnivora), pattern=map(redlist_binomial_woCarnivora, redlist_binomial_in_tree_woCarnivora, taxa_woCarnivora)),
  tar_target(shannon_index_all_summary_Carnivora, summarize_shannon_index_all(shannon_index, redlist_binomial_Carnivora, taxa_Carnivora)),
  tar_target(shannon_index_all_summary, bind_rows(shannon_index_all_summary_woCarnivora, shannon_index_all_summary_Carnivora)),
  tar_target(shannon_index_all_summary_figure, plot_shannon_index_all_summary(shannon_index_all_summary, threshold)),

  # Figure S5
  tar_target(shannon_index_tree_summary_woCarnivora, summarize_shannon_index_tree(shannon_index, redlist_binomial_woCarnivora, redlist_binomial_in_tree_woCarnivora, taxa_woCarnivora), pattern=map(redlist_binomial_woCarnivora, redlist_binomial_in_tree_woCarnivora, taxa_woCarnivora)),
  tar_target(shannon_index_tree_summary_Carnivora, summarize_shannon_index_tree(shannon_index, redlist_binomial_Carnivora, redlist_binomial_in_tree_Carnivora, taxa_Carnivora)),
  tar_target(shannon_index_tree_summary, bind_rows(shannon_index_tree_summary_woCarnivora, shannon_index_tree_summary_Carnivora)),
  tar_target(shannon_index_tree_summary_figure, plot_shannon_index_tree_summary(shannon_index_tree_summary, threshold)),

  # Figure S7
  tar_target(likelihood_summary_woCarnivora, summarize_likelihood(bisse_result_woCarnivora, taxa_woCarnivora), pattern=map(bisse_result_woCarnivora, taxa_woCarnivora)),
  tar_target(likelihood_summary_Carnivora, summarize_likelihood(bisse_result_Carnivora, taxa_Carnivora)),
  tar_target(likelihood_summary, bind_rows(likelihood_summary_woCarnivora, likelihood_summary_Carnivora)),
  tar_target(likelihood_summary_figure, plot_likelihood_summary(likelihood_summary, removed_step=removed_step)),

  # Figure S8
  tar_target(autocorrelation_summary_woCarnivora, summarize_autocorrelation(bisse_result_woCarnivora, taxa_woCarnivora, removed_step=removed_step, max_interval=sampling_interval*2), pattern=map(bisse_result_woCarnivora, taxa_woCarnivora)),
  tar_target(autocorrelation_summary_Carnivora, summarize_autocorrelation(bisse_result_Carnivora, taxa_Carnivora, removed_step=removed_step, max_interval=sampling_interval*2)),
  tar_target(autocorrelation_summary, bind_rows(autocorrelation_summary_woCarnivora, autocorrelation_summary_Carnivora)),
  tar_target(autocorrelation_summary_figure, plot_autocorrelation_summary(autocorrelation_summary, threshold=sampling_interval)),

  # Figure S9
  tar_target(filtered_transition_woCarnivora, filter_mcmc_transition(bisse_result_woCarnivora, taxa_woCarnivora, burn_in=removed_step, interval=sampling_interval), pattern=map(bisse_result_woCarnivora, taxa_woCarnivora)),
  tar_target(transition_summary_woCarnivora, summarize_samples(filtered_transition_woCarnivora, taxa_woCarnivora), pattern=map(filtered_transition_woCarnivora, taxa_woCarnivora)),
  tar_target(filtered_transition_Carnivora, filter_mcmc_transition(bisse_result_Carnivora, taxa_Carnivora, burn_in=removed_step, interval=sampling_interval)),
  tar_target(transition_summary_Carnivora, summarize_samples(filtered_transition_Carnivora, taxa_Carnivora)),
  tar_target(filtered_transition, bind_rows(filtered_transition_woCarnivora, filtered_transition_Carnivora)),
  tar_target(transition_summary, bind_rows(transition_summary_woCarnivora, transition_summary_Carnivora)),
  tar_target(transition_figure, plot_transition(filtered_transition, transition_summary)),

  # Figure S10
  tar_target(filtered_transition_diff_woCarnivora, filter_mcmc_transition_diff(bisse_result_woCarnivora, taxa_woCarnivora, burn_in=removed_step, interval=sampling_interval), pattern=map(bisse_result_woCarnivora, taxa_woCarnivora)),
  tar_target(transition_summary_diff_woCarnivora, summarize_samples_diff(filtered_transition_diff_woCarnivora, taxa_woCarnivora), pattern=map(filtered_transition_diff_woCarnivora, taxa_woCarnivora)),
  tar_target(filtered_transition_diff_Carnivora, filter_mcmc_transition_diff(bisse_result_Carnivora, taxa_Carnivora, burn_in=removed_step, interval=sampling_interval)),
  tar_target(transition_summary_diff_Carnivora, summarize_samples_diff(filtered_transition_diff_Carnivora, taxa_Carnivora)),
  tar_target(filtered_transition_diff, bind_rows(filtered_transition_diff_woCarnivora, filtered_transition_diff_Carnivora)),
  tar_target(transition_summary_diff, bind_rows(transition_summary_diff_woCarnivora, transition_summary_diff_Carnivora)),
  tar_target(transition_diff_figure, plot_transition_diff(filtered_transition_diff, transition_summary_diff)),

  # Figure S11
  tar_target(filtered_speciation_woCarnivora, filter_mcmc_speciation(bisse_result_woCarnivora, taxa_woCarnivora, burn_in=removed_step, interval=sampling_interval), pattern=map(bisse_result_woCarnivora, taxa_woCarnivora)),
  tar_target(speciation_summary_woCarnivora, summarize_samples(filtered_speciation_woCarnivora, taxa_woCarnivora), pattern=map(filtered_speciation_woCarnivora, taxa_woCarnivora)),
  tar_target(filtered_speciation_Carnivora, filter_mcmc_speciation(bisse_result_Carnivora, taxa_Carnivora, burn_in=removed_step, interval=sampling_interval)),
  tar_target(speciation_summary_Carnivora, summarize_samples(filtered_speciation_Carnivora, taxa_Carnivora)),
  tar_target(filtered_speciation, bind_rows(filtered_speciation_woCarnivora, filtered_speciation_Carnivora)),
  tar_target(speciation_summary, bind_rows(speciation_summary_woCarnivora, speciation_summary_Carnivora)),
  tar_target(speciation_figure, plot_speciation(filtered_speciation, speciation_summary)),

  # Figure S12
  tar_target(filtered_speciation_diff_woCarnivora, filter_mcmc_speciation_diff(bisse_result_woCarnivora, taxa_woCarnivora, burn_in=removed_step, interval=sampling_interval), pattern=map(bisse_result_woCarnivora, taxa_woCarnivora)),
  tar_target(speciation_summary_diff_woCarnivora, summarize_samples_diff(filtered_speciation_diff_woCarnivora, taxa_woCarnivora), pattern=map(filtered_speciation_diff_woCarnivora, taxa_woCarnivora)),
  tar_target(filtered_speciation_diff_Carnivora, filter_mcmc_speciation_diff(bisse_result_Carnivora, taxa_Carnivora, burn_in=removed_step, interval=sampling_interval)),
  tar_target(speciation_summary_diff_Carnivora, summarize_samples_diff(filtered_speciation_diff_Carnivora, taxa_Carnivora)),
  tar_target(filtered_speciation_diff, bind_rows(filtered_speciation_diff_woCarnivora, filtered_speciation_diff_Carnivora)),
  tar_target(speciation_summary_diff, bind_rows(speciation_summary_diff_woCarnivora, speciation_summary_diff_Carnivora)),
  tar_target(speciation_diff_figure, plot_speciation_diff(filtered_speciation_diff, speciation_summary_diff)),
  
  #
  tar_target(filtered_extinction_woCarnivora, filter_mcmc_extinction(bisse_result_woCarnivora, taxa_woCarnivora, burn_in=removed_step, interval=sampling_interval), pattern=map(bisse_result_woCarnivora, taxa_woCarnivora)),
  tar_target(extinction_summary_woCarnivora, summarize_samples(filtered_extinction_woCarnivora, taxa_woCarnivora), pattern=map(filtered_extinction_woCarnivora, taxa_woCarnivora)),
  tar_target(filtered_extinction_Carnivora, filter_mcmc_extinction(bisse_result_Carnivora, taxa_Carnivora, burn_in=removed_step, interval=sampling_interval)),
  tar_target(extinction_summary_Carnivora, summarize_samples(filtered_extinction_Carnivora, taxa_Carnivora)),
  tar_target(filtered_extinction, bind_rows(filtered_extinction_woCarnivora, filtered_extinction_Carnivora)),
  tar_target(extinction_summary, bind_rows(extinction_summary_woCarnivora, extinction_summary_Carnivora)),
  tar_target(extinction_figure, plot_extinction(filtered_extinction, extinction_summary)),

  #
  tar_target(filtered_extinction_diff_woCarnivora, filter_mcmc_extinction_diff(bisse_result_woCarnivora, taxa_woCarnivora, burn_in=removed_step, interval=sampling_interval), pattern=map(bisse_result_woCarnivora, taxa_woCarnivora)),
  tar_target(extinction_summary_diff_woCarnivora, summarize_samples_diff(filtered_extinction_diff_woCarnivora, taxa_woCarnivora), pattern=map(filtered_extinction_diff_woCarnivora, taxa_woCarnivora)),
  tar_target(filtered_extinction_diff_Carnivora, filter_mcmc_extinction_diff(bisse_result_Carnivora, taxa_Carnivora, burn_in=removed_step, interval=sampling_interval)),
  tar_target(extinction_summary_diff_Carnivora, summarize_samples_diff(filtered_extinction_diff_Carnivora, taxa_Carnivora)),
  tar_target(filtered_extinction_diff, bind_rows(filtered_extinction_diff_woCarnivora, filtered_extinction_diff_Carnivora)),
  tar_target(extinction_summary_diff, bind_rows(extinction_summary_diff_woCarnivora, extinction_summary_diff_Carnivora)),
  tar_target(extinction_diff_figure, plot_extinction_diff(filtered_extinction_diff, extinction_summary_diff)),

  # 
  tar_target(number_summary_woCarnivora, summarize_number(character, redlist_binomial_woCarnivora, redlist_binomial_in_tree_woCarnivora, taxa_woCarnivora), pattern=map(redlist_binomial_woCarnivora, redlist_binomial_in_tree_woCarnivora, taxa_woCarnivora)),
  tar_target(number_summary_Carnivora, summarize_number(character, redlist_binomial_Carnivora, redlist_binomial_in_tree_Carnivora, taxa_Carnivora)),
  tar_target(number_summary, bind_rows(number_summary_woCarnivora, number_summary_Carnivora)),
  tar_target(number_summary_table, output_summarized_number(number_summary)),

  # DR statistic and FiSSE
  tar_target(taxa, c("Mammalia", "Rodentia", "Primates", "Cetartiodactyla", "Chiroptera", "Eulipotyphla", "Marsupialia", "Carnivora")),
  tar_target(redlist_binomial, extract_species_from_redlist(character, taxa), pattern=map(taxa)),
  tar_target(redlist_binomial_in_tree, map_redlist_to_tree(redlist_binomial, tree_binomial), pattern=map(redlist_binomial)),
  tar_target(dr_statistic, calculate_dr_statistic(tree, character, redlist_binomial_in_tree, taxa), pattern=map(redlist_binomial_in_tree, taxa)),
  tar_target(dr_statistic_figure, plot_dr_statistic(dr_statistic, taxa), pattern=map(dr_statistic, taxa)),
  tar_target(fisse_result, run_fisse(tree, character, redlist_binomial_in_tree, taxa), pattern=map(redlist_binomial_in_tree, taxa)),

  # DR statistic from complete trees
  tar_target(dr_stat_complete_trees, "data/tree/DR-SUMMARY_MamPhy_BDvr_Completed_5911sp_topoCons_NDexp_all10k_v2_expanded.txt"),
  tar_target(dr_stat_mean_complete_trees, extract_dr_stat_complete_trees(dr_stat_complete_trees, character, redlist_binomial_in_tree, taxa), pattern=map(redlist_binomial_in_tree, taxa)),
  tar_target(dr_stat_mean_complete_trees_figure, plot_dr_stat_complete_trees(dr_stat_mean_complete_trees, character, taxa), pattern=map(dr_stat_mean_complete_trees, taxa)),
  tar_target(dr_stat_all_figure, plot_dr_stat_all(dr_stat_mean_complete_trees, dr_statistic, character, taxa), pattern=map(dr_stat_mean_complete_trees, dr_statistic, taxa))
)