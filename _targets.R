source("R/packages.R")
source("R/functions.R")
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
  tar_target(bisse_figure_woCarnivora, plot_speciation_transition(bisse_result_woCarnivora, taxa_woCarnivora, burn_in=removed_step, interval=sampling_interval), pattern=map(bisse_result_woCarnivora, taxa_woCarnivora)),
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
  tar_target(bisse_figure_Carnivora, plot_speciation_transition(bisse_result_Carnivora, taxa_Carnivora, burn_in=removed_step, interval=sampling_interval)),
  tar_target(pruned_tree_figure_Carnivora, plot_tree(pruned_tree_Carnivora, character, redlist_binomial_in_tree_Carnivora, taxa_Carnivora)),
  
  # Visualize results from these mammalian lineages at the same time
  # Figure S2
  tar_target(shannon_index_all_summary_woCarnivora, summarize_shannon_index_all(shannon_index, redlist_binomial_woCarnivora, taxa_woCarnivora), pattern=map(redlist_binomial_woCarnivora, redlist_binomial_in_tree_woCarnivora, taxa_woCarnivora)),
  tar_target(shannon_index_all_summary_Carnivora, summarize_shannon_index_all(shannon_index, redlist_binomial_Carnivora, taxa_Carnivora)),
  tar_target(shannon_index_all_summary, bind_rows(shannon_index_all_summary_woCarnivora, shannon_index_all_summary_Carnivora)),
  tar_target(shannon_index_all_summary_figure, plot_shannon_index_all_summary(shannon_index_all_summary, threshold)),

  # Figure S4
  tar_target(shannon_index_tree_summary_woCarnivora, summarize_shannon_index_tree(shannon_index, redlist_binomial_woCarnivora, redlist_binomial_in_tree_woCarnivora, taxa_woCarnivora), pattern=map(redlist_binomial_woCarnivora, redlist_binomial_in_tree_woCarnivora, taxa_woCarnivora)),
  tar_target(shannon_index_tree_summary_Carnivora, summarize_shannon_index_tree(shannon_index, redlist_binomial_Carnivora, redlist_binomial_in_tree_Carnivora, taxa_Carnivora)),
  tar_target(shannon_index_tree_summary, bind_rows(shannon_index_tree_summary_woCarnivora, shannon_index_tree_summary_Carnivora)),
  tar_target(shannon_index_tree_summary_figure, plot_shannon_index_tree_summary(shannon_index_tree_summary, threshold)),

  # Figure S5
  tar_target(likelihood_summary_woCarnivora, summarize_likelihood(bisse_result_woCarnivora, taxa_woCarnivora), pattern=map(bisse_result_woCarnivora, taxa_woCarnivora)),
  tar_target(likelihood_summary_Carnivora, summarize_likelihood(bisse_result_Carnivora, taxa_Carnivora)),
  tar_target(likelihood_summary, bind_rows(likelihood_summary_woCarnivora, likelihood_summary_Carnivora)),
  tar_target(likelihood_summary_figure, plot_likelihood_summary(likelihood_summary, removed_step=removed_step)),

  # Figure S6
  tar_target(autocorrelation_summary_woCarnivora, summarize_autocorrelation(bisse_result_woCarnivora, taxa_woCarnivora, removed_step=removed_step, max_interval=sampling_interval*2), pattern=map(bisse_result_woCarnivora, taxa_woCarnivora)),
  tar_target(autocorrelation_summary_Carnivora, summarize_autocorrelation(bisse_result_Carnivora, taxa_Carnivora, removed_step=removed_step, max_interval=sampling_interval*2)),
  tar_target(autocorrelation_summary, bind_rows(autocorrelation_summary_woCarnivora, autocorrelation_summary_Carnivora)),
  tar_target(autocorrelation_summary_figure, plot_autocorrelation_summary(autocorrelation_summary, threshold=sampling_interval)),

  # 
  tar_target(filtered_speciation_woCarnivora, filter_mcmc_speciation(bisse_result_woCarnivora, taxa_woCarnivora, burn_in=removed_step, interval=sampling_interval), pattern=map(bisse_result_woCarnivora, taxa_woCarnivora)),
  tar_target(speciation_summary_woCarnivora, summarize_samples(filtered_speciation_woCarnivora, taxa_woCarnivora), pattern=map(filtered_speciation_woCarnivora, taxa_woCarnivora)),
  tar_target(filtered_speciation_Carnivora, filter_mcmc_speciation(bisse_result_Carnivora, taxa_Carnivora, burn_in=removed_step, interval=sampling_interval)),
  tar_target(speciation_summary_Carnivora, summarize_samples(filtered_speciation_Carnivora, taxa_Carnivora)),
  tar_target(filtered_speciation, bind_rows(filtered_speciation_woCarnivora, filtered_speciation_Carnivora)),
  tar_target(speciation_summary, bind_rows(speciation_summary_woCarnivora, speciation_summary_Carnivora)),
  tar_target(speciation_figure, plot_speciation(filtered_speciation, speciation_summary)),

  # Figure S7
  tar_target(filtered_speciation_diff_woCarnivora, filter_mcmc_speciation_diff(bisse_result_woCarnivora, taxa_woCarnivora, burn_in=removed_step, interval=sampling_interval), pattern=map(bisse_result_woCarnivora, taxa_woCarnivora)),
  tar_target(speciation_summary_diff_woCarnivora, summarize_samples_diff(filtered_speciation_diff_woCarnivora, taxa_woCarnivora), pattern=map(filtered_speciation_diff_woCarnivora, taxa_woCarnivora)),
  tar_target(filtered_speciation_diff_Carnivora, filter_mcmc_speciation_diff(bisse_result_Carnivora, taxa_Carnivora, burn_in=removed_step, interval=sampling_interval)),
  tar_target(speciation_summary_diff_Carnivora, summarize_samples_diff(filtered_speciation_diff_Carnivora, taxa_Carnivora)),
  tar_target(filtered_speciation_diff, bind_rows(filtered_speciation_diff_woCarnivora, filtered_speciation_diff_Carnivora)),
  tar_target(speciation_summary_diff, bind_rows(speciation_summary_diff_woCarnivora, speciation_summary_diff_Carnivora)),
  tar_target(speciation_diff_figure, plot_speciation_diff(filtered_speciation_diff, speciation_summary_diff)),

  # 
  tar_target(filtered_transition_woCarnivora, filter_mcmc_transition(bisse_result_woCarnivora, taxa_woCarnivora, burn_in=removed_step, interval=sampling_interval), pattern=map(bisse_result_woCarnivora, taxa_woCarnivora)),
  tar_target(transition_summary_woCarnivora, summarize_samples(filtered_transition_woCarnivora, taxa_woCarnivora), pattern=map(filtered_transition_woCarnivora, taxa_woCarnivora)),
  tar_target(filtered_transition_Carnivora, filter_mcmc_transition(bisse_result_Carnivora, taxa_Carnivora, burn_in=removed_step, interval=sampling_interval)),
  tar_target(transition_summary_Carnivora, summarize_samples(filtered_transition_Carnivora, taxa_Carnivora)),
  tar_target(filtered_transition, bind_rows(filtered_transition_woCarnivora, filtered_transition_Carnivora)),
  tar_target(transition_summary, bind_rows(transition_summary_woCarnivora, transition_summary_Carnivora)),
  tar_target(transition_figure, plot_transition(filtered_transition, transition_summary)),

  # Figure S8
  tar_target(filtered_transition_diff_woCarnivora, filter_mcmc_transition_diff(bisse_result_woCarnivora, taxa_woCarnivora, burn_in=removed_step, interval=sampling_interval), pattern=map(bisse_result_woCarnivora, taxa_woCarnivora)),
  tar_target(transition_summary_diff_woCarnivora, summarize_samples_diff(filtered_transition_diff_woCarnivora, taxa_woCarnivora), pattern=map(filtered_transition_diff_woCarnivora, taxa_woCarnivora)),
  tar_target(filtered_transition_diff_Carnivora, filter_mcmc_transition_diff(bisse_result_Carnivora, taxa_Carnivora, burn_in=removed_step, interval=sampling_interval)),
  tar_target(transition_summary_diff_Carnivora, summarize_samples_diff(filtered_transition_diff_Carnivora, taxa_Carnivora)),
  tar_target(filtered_transition_diff, bind_rows(filtered_transition_diff_woCarnivora, filtered_transition_diff_Carnivora)),
  tar_target(transition_summary_diff, bind_rows(transition_summary_diff_woCarnivora, transition_summary_diff_Carnivora)),
  tar_target(transition_diff_figure, plot_transition_diff(filtered_transition_diff, transition_summary_diff)),

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

  # ---
  # Aves
  # Aves (Hackett, MCC):
  # Generate a MCC tree from 10 000 trees sampled in the pseudo-posterior distribution using TreeAnnotator
  # Concatenate tree files:
  ## cat data/Aves/mnt/data/projects/birdphylo/Tree_sets/Stage1_full_data/CombinedTrees/HackettStage1Full_{1..10}.tre > data/Aves/mnt/data/projects/birdphylo/Tree_sets/Stage1_full_data/CombinedTrees/HackettStage1Full_1-10.tre
  # Convert newick to nexus (in R):
  ## library(ape)
  ## tree_path <- "data/Aves/mnt/data/projects/birdphylo/Tree_sets/Stage1_full_data/CombinedTrees/HackettStage1Full_1-10.tre"
  ## treelines <- readLines(tree_path, n = 10000)
  ## tree <- read.tree(text = treelines)
  ## ape::write.nexus(tree, file="data/Aves/mnt/data/projects/birdphylo/Tree_sets/Stage1_full_data/CombinedTrees/HackettStage1Full_1-10.nex")
  # Install BEAST: https://beast.community/installing
  # Edit BEASTv1.10.4/bin/treeannotator and expand heap (-Xms6144m -Xmx6144m)
  # Run TreeAnnotator
  ## treeannotator data/Aves/mnt/data/projects/birdphylo/Tree_sets/Stage1_full_data/CombinedTrees/HackettStage1Full_1-10.nex data/Aves/mnt/data/projects/birdphylo/Tree_sets/Stage1_full_data/CombinedTrees/HackettStage1Full_1-10_MCC.tre

  # Download BOTW.7z from birdlife.org and extract it beforehand.

  tar_target(tree_Aves, "data/tree/HackettStage1Full_1-10_MCC.tre", format="file"),
  tar_target(redlist_original_Aves, "data/Aves/a00000009.gdbtable", format="file"),
  tar_target(redlist_preprocessed_Aves, preprocessing_redlist_Aves(redlist_original_Aves)),
  tar_target(redlist_climate_original_Aves, merge_redlist_with_climate_Aves(redlist_preprocessed_Aves, climate_map, confidence_map)), # This step will take about a week
  tar_target(redlist_climate_Aves, summarize_area_by_species(redlist_climate_original_Aves)),
  tar_target(shannon_index_Aves, calculate_shannon_index(redlist_climate_Aves)),
  tar_target(threshold_Aves, calculate_threshold(shannon_index_Aves)),
  tar_target(character_Aves, annotate_specialist_generalist(shannon_index_Aves, threshold_Aves)),
  tar_target(redlist_binomial_Aves, extract_species_from_redlist_Aves(character_Aves)),
  tar_target(tree_binomial_Aves, extract_species_from_tree_Aves(tree_Aves)),
  tar_target(redlist_binomial_in_tree_Aves, map_redlist_to_tree(redlist_binomial_Aves, tree_binomial_Aves)),
  tar_target(pruned_tree_Aves, prune_tree(tree_Aves, redlist_binomial_in_tree_Aves)),
  tar_target(sampling_fraction_Aves, calculate_sampling_fraction(character_Aves, redlist_binomial_Aves, redlist_binomial_in_tree_Aves)),
  tar_target(bisse_result_Aves, run_bisse(character_Aves, pruned_tree_Aves, redlist_binomial_in_tree_Aves, sampling_fraction_Aves, iteration=sampling_step)),
  tar_target(shannon_index_figure_Aves, plot_shannon_index(shannon_index_Aves, threshold_Aves, redlist_binomial_Aves, redlist_binomial_in_tree_Aves, "Aves")),
  tar_target(likelihood_figure_Aves, plot_likelihood(bisse_result_Aves, "Aves", removed_step=removed_step)),
  tar_target(autocorrelation_figure_Aves, plot_autocorrelation(bisse_result_Aves, "Aves", removed_step=removed_step, max_interval=sampling_interval*2, threshold=sampling_interval)),
  tar_target(bisse_figure_Aves, plot_bisse(bisse_result_Aves, "Aves", burn_in=removed_step, interval=sampling_interval)),
  tar_target(number_summary_table_Aves, summarize_number_Aves(character_Aves, redlist_binomial_Aves, redlist_binomial_in_tree_Aves, "Aves")),
  tar_target(pruned_tree_figure_Aves, plot_tree_Aves(pruned_tree_Aves, character_Aves, redlist_binomial_in_tree_Aves, "Aves")),
  tar_target(bisse_diff_figure_Aves, plot_bisse_diff(bisse_result_Aves, "Aves", burn_in=removed_step, interval=sampling_interval)),

  # ---
  # Typical clades
  # Figure 4a
  tar_target(typical_taxa_woCarnivora, c("Rodentia", "Primates", "Cetartiodactyla")),
  tar_target(filtered_speciation_typical_woCarnivora, filtered_speciation_woCarnivora, pattern=slice(filtered_speciation_woCarnivora, index=c(2, 3, 4))),
  tar_target(speciation_summary_typical_woCarnivora, summarize_samples(filtered_speciation_typical_woCarnivora, typical_taxa_woCarnivora), pattern=map(filtered_speciation_typical_woCarnivora, typical_taxa_woCarnivora)),
  tar_target(filtered_speciation_typical, bind_rows(filtered_speciation_typical_woCarnivora, filtered_speciation_Carnivora)),
  tar_target(speciation_summary_typical, bind_rows(speciation_summary_typical_woCarnivora, speciation_summary_Carnivora)),
  tar_target(speciation_figure_typical, plot_speciation_typical(filtered_speciation_typical, speciation_summary_typical)),

  # Figure 4b
  tar_target(filtered_transition_typical_woCarnivora, filtered_transition_woCarnivora, pattern=slice(filtered_transition_woCarnivora, index=c(2, 3, 4))),
  tar_target(transition_summary_typical_woCarnivora, summarize_samples(filtered_transition_typical_woCarnivora, typical_taxa_woCarnivora), pattern=map(filtered_transition_typical_woCarnivora, typical_taxa_woCarnivora)),
  tar_target(filtered_transition_typical, bind_rows(filtered_transition_typical_woCarnivora, filtered_transition_Carnivora)),
  tar_target(transition_summary_typical, bind_rows(transition_summary_typical_woCarnivora, transition_summary_Carnivora)),
  tar_target(transition_figure_typical, plot_transition_typical(filtered_transition_typical, transition_summary_typical)),

  # Atypical clades
  # Figure 5a
  tar_target(atypical_mammal, c("Chiroptera", "Eulipotyphla", "Marsupialia")),
  tar_target(filtered_speciation_atypical_mammal, filtered_speciation_woCarnivora, pattern=slice(filtered_speciation_woCarnivora, index=c(5, 6, 7))),
  tar_target(speciation_summary_atypical_mammal, summarize_samples(filtered_speciation_atypical_mammal, atypical_mammal), pattern=map(filtered_speciation_atypical_mammal, atypical_mammal)),
  tar_target(filtered_speciation_Aves, filter_mcmc_speciation(bisse_result_Aves, "Aves", burn_in=removed_step, interval=sampling_interval)),
  tar_target(speciation_summary_Aves, summarize_samples(filtered_speciation_Aves, "Aves")),
  tar_target(atypical_taxa, c("Chiroptera", "Eulipotyphla", "Marsupialia", "Aves")),
  tar_target(filtered_speciation_atypical, bind_rows(filtered_speciation_atypical_mammal, filtered_speciation_Aves)),
  tar_target(speciation_summary_atypical, bind_rows(speciation_summary_atypical_mammal, speciation_summary_Aves)),
  tar_target(speciation_figure_atypical, plot_speciation_atypical(filtered_speciation_atypical, speciation_summary_atypical)),

  # Figure 5b
  tar_target(filtered_transition_atypical_mammal, filtered_transition_woCarnivora, pattern=slice(filtered_transition_woCarnivora, index=c(5, 6, 7))),
  tar_target(transition_summary_atypical_mammal, summarize_samples(filtered_transition_atypical_mammal, atypical_mammal), pattern=map(filtered_transition_atypical_mammal, atypical_mammal)),
  tar_target(filtered_transition_Aves, filter_mcmc_transition(bisse_result_Aves, "Aves", burn_in=removed_step, interval=sampling_interval)),
  tar_target(transition_summary_Aves, summarize_samples(filtered_transition_Aves, "Aves")),
  tar_target(filtered_transition_atypical, bind_rows(filtered_transition_atypical_mammal, filtered_transition_Aves)),
  tar_target(transition_summary_atypical, bind_rows(transition_summary_atypical_mammal, transition_summary_Aves)),
  tar_target(transition_figure_atypical, plot_transition_atypical(filtered_transition_atypical, transition_summary_atypical))
)