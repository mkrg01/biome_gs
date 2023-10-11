# Download MCC consensus tree of the DNA-only distributions
download_tree <- function() {
  download.file(
    url = "https://raw.githubusercontent.com/n8upham/MamPhy_v1/master/_DATA/MamPhy_fullPosterior_BDvr_DNAonly_4098sp_topoFree_NDexp_MCC_v2_target.tre",
    destfile = "data/tree/MamPhy_fullPosterior_BDvr_DNAonly_4098sp_topoFree_NDexp_MCC_v2_target.tre"
  )
  return("data/tree/MamPhy_fullPosterior_BDvr_DNAonly_4098sp_topoFree_NDexp_MCC_v2_target.tre")
}


preprocessing_redlist <- function(redlist_path) {
  # Read polygon data
  redlist <- st_read(redlist_path)

  # Use "Extant", "Possibly Extinct", and "Extinct" distribution
  redlist <- filter(redlist, presence==1 | presence==4 | presence==5)

  # Use only "Native" distribution
  redlist <- filter(redlist, origin==1)

  # Use only terrestrial distribution
  redlist <- filter(redlist, terrestial=="true")

  return(redlist)
}


# Merge IUCN Red List polygon data with Köppen-Geiger climate classification map
merge_redlist_with_climate <- function(redlist, climate_map_path, confidence_map_path) {
  # Read raster files
  climate_map <- raster(climate_map_path)
  confidence_map <- raster(confidence_map_path)

  # Use only 100% confidence cells in the map
  kg_conf_binary <- calc(confidence_map, fun=function(x) {x==100})
  kg_100conf <- climate_map * kg_conf_binary

  # Layerize the map for parallelization by exactextractr
  kg_brick <- layerize(kg_100conf)

  # Merge polygon data with climate classes
  redlist <- cbind(redlist, exact_extract(kg_brick, redlist, 'weighted_sum', weights = area(kg_brick), stack_apply = TRUE))
  redlist <- st_drop_geometry(redlist)

  # Change colnames and select some columns
  tmp <- data.frame(binomial="temporary", weighted_sum.X0=0 , weighted_sum.X1=0, weighted_sum.X2=0, weighted_sum.X3=0, weighted_sum.X4=0, weighted_sum.X5=0, weighted_sum.X6=0, weighted_sum.X7=0, weighted_sum.X8=0, weighted_sum.X9=0, weighted_sum.X10=0, weighted_sum.X11=0, weighted_sum.X12=0, weighted_sum.X13=0, weighted_sum.X14=0, weighted_sum.X15=0, weighted_sum.X16=0, weighted_sum.X17=0, weighted_sum.X18=0, weighted_sum.X19=0, weighted_sum.X20=0, weighted_sum.X21=0, weighted_sum.X22=0, weighted_sum.X23=0, weighted_sum.X24=0, weighted_sum.X25=0, weighted_sum.X26=0, weighted_sum.X27=0, weighted_sum.X28=0, weighted_sum.X29=0, weighted_sum.X30=0)
  redlist <- merge(redlist, tmp, all=T) %>%
    as_tibble() %>%
    filter(binomial != "temporary") %>%
    mutate_if(is.numeric, ~tidyr::replace_na(., 0))
  redlist <- redlist %>%
    dplyr::rename(N = weighted_sum.X0, Af = weighted_sum.X1, Am = weighted_sum.X2, Aw = weighted_sum.X3, BWh = weighted_sum.X4, BWk = weighted_sum.X5, BSh = weighted_sum.X6, BSk = weighted_sum.X7, Csa = weighted_sum.X8, Csb = weighted_sum.X9, Csc = weighted_sum.X10, Cwa = weighted_sum.X11, Cwb = weighted_sum.X12, Cwc = weighted_sum.X13, Cfa = weighted_sum.X14, Cfb = weighted_sum.X15, Cfc = weighted_sum.X16, Dsa = weighted_sum.X17, Dsb = weighted_sum.X18, Dsc = weighted_sum.X19, Dsd = weighted_sum.X20, Dwa = weighted_sum.X21, Dwb = weighted_sum.X22, Dwc = weighted_sum.X23, Dwd = weighted_sum.X24, Dfa = weighted_sum.X25, Dfb = weighted_sum.X26, Dfc = weighted_sum.X27, Dfd = weighted_sum.X28, ET = weighted_sum.X29, EF = weighted_sum.X30)
  kg_list <- c("N", "Af", "Am", "Aw", "BWh", "BWk", "BSh", "BSk", "Csa", "Csb", "Csc", "Cwa", "Cwb", "Cwc", "Cfa", "Cfb", "Cfc", "Dsa", "Dsb", "Dsc", "Dsd", "Dwa", "Dwb", "Dwc", "Dwd", "Dfa", "Dfb", "Dfc", "Dfd", "ET", "EF")
  redlist$species <- redlist$binomial
  redlist$order <- redlist$order_
  redlist <- dplyr::select(redlist, c("species", all_of(kg_list)), order)

  return(redlist)
}


summarize_area_by_species <- function(redlist_climate) {
  # Summarize by species
  by_species <- dplyr::group_by(redlist_climate, species, order)
  redlist_climate <- dplyr::summarise_all(by_species, sum)
  redlist_climate <- ungroup(redlist_climate)
  return(redlist_climate)
}


calculate_shannon_index <- function(redlist) {
  # Preprocessing
  kg_list <- c("Af", "Am", "Aw", "BWh", "BWk", "BSh", "BSk", "Csa", "Csb", "Csc", "Cwa", "Cwb", "Cwc", "Cfa", "Cfb", "Cfc", "Dsa", "Dsb", "Dsc", "Dsd", "Dwa", "Dwb", "Dwc", "Dwd", "Dfa", "Dfb", "Dfc", "Dfd", "ET", "EF")
  redlist <- redlist %>%
    mutate(sum = rowSums(redlist[,kg_list])) %>%
    filter(sum > 0)
  
  kg_list <- c("A", "B", "C", "D", "E")
  redlist <- redlist %>%
    mutate(A = Af + Am + Aw, B = BWh + BWk + BSh + BSk, C = Csa + Csb + Csc + Cwa + Cwb + Cwc + Cfa + Cfb + Cfc, D = Dsa + Dsb + Dsc + Dsd + Dwa + Dwb + Dwc + Dwd + Dfa + Dfb + Dfc + Dfd, E = ET + EF) %>%
    dplyr::select(species, all_of(kg_list), sum, order)
  
  # Calculate Shannon index
  redlist <- redlist %>%
    mutate_at(kg_list, function(x) x / redlist$sum) %>%
    mutate_at(kg_list, function(x) {ifelse(x==0, 0, -x * log(x))})
  redlist <- redlist %>%
    mutate(shannon_index = rowSums(redlist[,kg_list])) %>%
    dplyr::select(species, shannon_index, order)
  
  return(redlist)
}


# Calculate the mean of shannon index
calculate_threshold <- function(shannon_index) {
  return(mean(shannon_index$shannon_index))
}


# Annotate biome specialists and generalists
annotate_specialist_generalist <- function(redlist_shannon, threshold) {
  redlist_shannon <- redlist_shannon %>%
    mutate(generalist = shannon_index > threshold) %>%
    mutate_at(vars(generalist), as.numeric) %>%
    dplyr::select(species, generalist, order)
  return(redlist_shannon)
}


# Extract binomial names from a tree
extract_species_from_tree <- function(tree_path) {
  tree <- ape::read.nexus(tree_path)

  # Remove outgroup
  tree <- ape::drop.tip(tree, "_Anolis_carolinensis")

  # Extract binomial name
  spl <- stringr::str_split(tree$tip.label, "_")
  leaf_list_renamed <- c()
  for (i in 1:length(tree$tip.label)) leaf_list_renamed <- append(leaf_list_renamed, paste(spl[i][[1]][1], spl[i][[1]][2]))
  tree$tip.label <- leaf_list_renamed
  return(tree$tip.label)
}


# Extract species names from IUCN Red List
extract_species_from_redlist <- function(redlist, taxon) {
  if (taxon == "Mammalia") {
    return(redlist$species)
  } else if (taxon == "Marsupialia") {
    order_list <- c("PAUCITUBERCULATA", "DIDELPHIMORPHIA", "MICROBIOTHERIA", "NOTORYCTEMORPHIA", "PERAMELEMORPHIA", "DASYUROMORPHIA", "DIPROTODONTIA")
    redlist <- redlist %>%
      filter(order %in% order_list)
    return(redlist$species)
  } else {
    order_list <- toupper(taxon)
    redlist <- redlist %>%
      filter(order %in% order_list)
    return(redlist$species)
  }
}


# Map species in Red List into external nodes in a phylogenetic tree
map_redlist_to_tree <- function(redlist_species, tree_species) {
  redlist_species <- tibble(species = redlist_species) %>%
    filter(species %in% tree_species)
  return(redlist_species$species)
}


# Prune a tree for BiSSE
prune_tree <- function(tree_path, species) {
  tree <- read.nexus(tree_path)
  spl <- stringr::str_split(tree$tip.label, "_")
  leaf_list_renamed <- c()
  for (i in 1:length(tree$tip.label)) leaf_list_renamed <- append(leaf_list_renamed, paste(spl[i][[1]][1], spl[i][[1]][2]))
  tree$tip.label <- gsub(" ", "_", leaf_list_renamed)
  tree <- keep.tip(tree, gsub(" ", "_", species))
  if (is.ultrametric(tree) == F) {
    tree <- force.ultrametric(tree)
  }
  return(tree)
}

calculate_sampling_fraction <- function(character, redlist_species, redlist_species_in_tree) {
  species_all <- character %>% 
    filter(species %in% redlist_species)
  species_tree <- character %>% 
    filter(species %in% redlist_species_in_tree)
  
  specialist_all <- species_all %>% 
    filter(generalist == 0) %>% 
    .$species %>% 
    length()
  generalist_all <- length(species_all$species) - specialist_all
  specialist_tree <- species_tree %>% 
    filter(generalist == 0) %>% 
    .$species %>% 
    length()
  generalist_tree <- length(species_tree$species) - specialist_tree
  
  sampling_fraction_specialist <- specialist_tree / specialist_all
  sampling_fraction_generalist <- generalist_tree / generalist_all
  return(c(sampling_fraction_specialist, sampling_fraction_generalist))
}


run_bisse <- function(character, pruned_tree, redlist_species_in_tree, sampling.f, iteration) {
  # Use only species in the tree
  character <- character %>%
    filter(species %in% redlist_species_in_tree)

  # Run BiSSE
  states <- setNames(character$generalist, gsub(" ", "_", character$species))
  lik <- make.bisse(pruned_tree, states, sampling.f=sampling.f)
  p <- starting.point.bisse(pruned_tree)
  fit <- find.mle(lik, p, control=list(maxit=20000))
  prior <- make.prior.exponential(1 / (2 * (p[1] - p[3])))
  tmp <- mcmc(lik, fit$par, nsteps=100, prior=prior, lower=0, w=.1, print.every=0)
  w <- diff(sapply(tmp[2:7], quantile, c(.05, .95)))
  samples <- mcmc(lik, fit$par, nsteps=iteration, w=w, lower=0, prior=prior, print.every=500)
  
  return(samples)
}


run_bisse_skip_ml <- function(character, pruned_tree, redlist_species_in_tree, sampling.f, iteration) {
  # Use only species in the tree
  character <- character %>% 
    filter(species %in% redlist_species_in_tree)

  # Run BiSSE (I skipped ML estimation to avoid the likelihood calculation failure in Carnivora)
  states <- setNames(character$generalist, gsub(" ", "_", character$species))
  lik <- make.bisse(pruned_tree, states, sampling.f=sampling.f)
  p <- starting.point.bisse(pruned_tree)
  prior <- make.prior.exponential(1 / (2 * (p[1] - p[3])))
  samples <- mcmc(lik, p, nsteps=iteration, w=0.1, lower=0, prior=prior, print.every=500)
  
  return(samples)
}


# Visualize Shannon index
plot_shannon_index <- function(shannon_index, threshold, redlist_species, redlist_species_in_tree, taxon) {
  # Select a particular taxon
  shannon_index <- shannon_index %>% 
    filter(species %in% redlist_species)

  shannon_index_all <- shannon_index %>%
    dplyr::select(shannon_index)
  shannon_index_tree <- shannon_index %>%
    filter(species %in% redlist_species_in_tree) %>%
    dplyr::select(shannon_index)
  
  colnames(shannon_index_all) <- c("All")
  colnames(shannon_index_tree) <- c("Tree")

  shannon_index_all <- shannon_index_all %>%
    gather(key=case, value=shannon_index)
  shannon_index_tree <- shannon_index_tree %>%
    gather(key=case, value=shannon_index)
  
  shannon_index_concat <- bind_rows(shannon_index_all, shannon_index_tree)

  # Plot
  p <- ggplot(shannon_index_concat) + 
    geom_histogram(aes(shannon_index), binwidth=0.05) + 
    geom_vline(aes(xintercept=threshold), color="red", linetype="dashed") +
    facet_wrap(~case, ncol=2) +
    ggtitle(taxon) +
    xlab("Shannon index") +
    ylab("Number of species")
  
  ggsave(paste0("data/output/", taxon, "/shannon_index.pdf"), plot=p, width=5, height=3)
  return(paste0("data/output/", taxon, "/shannon_index.pdf"))
}


plot_shannon_index_all <- function(shannon_index, threshold, redlist_species, taxon) {
  # Select a particular taxon
  shannon_index <- shannon_index %>% 
    filter(species %in% redlist_species)

  shannon_index_all <- shannon_index %>%
    dplyr::select(shannon_index)
  
  colnames(shannon_index_all) <- c("All")

  shannon_index_all <- shannon_index_all %>%
    gather(key=case, value=shannon_index)

  # Plot
  p <- ggplot(shannon_index_all) + 
    geom_histogram(aes(shannon_index), binwidth=0.05) + 
    geom_vline(aes(xintercept=threshold), color="red", linetype="dashed") +
    xlab("Shannon index") +
    ylab("Number of species")
  
  ggsave(paste0("data/output/", taxon, "/shannon_index_all.pdf"), plot=p, width=3, height=3)
  return(paste0("data/output/", taxon, "/shannon_index_all.pdf"))
}


# Concatenate figures of Shannon index for species in the tree
summarize_shannon_index_all <- function(shannon_index, redlist_species, taxa) {
  shannon_index_all <- shannon_index %>% 
    filter(species %in% redlist_species) %>% 
    dplyr::select(shannon_index)
  
  colnames(shannon_index_all) <- c(taxa)
  shannon_index_all <- shannon_index_all %>% 
    gather(key=case, value=shannon_index)
  
  return(shannon_index_all)
}


plot_shannon_index_all_summary <- function(shannon_index_all, threshold) {
  shannon_index_all$case <- factor(shannon_index_all$case, levels=c("Mammalia", "Rodentia", "Primates", "Cetartiodactyla", "Carnivora", "Chiroptera", "Eulipotyphla", "Marsupialia"))
  p <- ggplot(shannon_index_all) + 
    geom_histogram(aes(shannon_index), binwidth=0.05) + 
    geom_vline(aes(xintercept=threshold), color="red", linetype="dashed") +
    facet_wrap(~case, ncol=4, scales="free_y") +
    ggtitle("All") +
    xlab("Shannon index") +
    ylab("Number of species")
  
  ggsave("data/output/Mammalia/shannon_index_all_summary.pdf", plot=p, width=8, height=4)
  return("data/output/Mammalia/shannon_index_all_summary.pdf")
}


summarize_shannon_index_tree <- function(shannon_index, redlist_species, redlist_species_in_tree, taxa) {
  shannon_index_taxon <- shannon_index %>% 
    filter(species %in% redlist_species)
  shannon_index_tree <- shannon_index_taxon %>%
    filter(species %in% redlist_species_in_tree) %>%
    dplyr::select(shannon_index)
  
  colnames(shannon_index_tree) <- c(taxa)
  shannon_index_tree <- shannon_index_tree %>% 
    gather(key=case, value=shannon_index)
  
  return(shannon_index_tree)
}


plot_shannon_index_tree_summary <- function(shannon_index_tree, threshold) {
  shannon_index_tree$case <- factor(shannon_index_tree$case, levels=c("Mammalia", "Rodentia", "Primates", "Cetartiodactyla", "Carnivora", "Chiroptera", "Eulipotyphla", "Marsupialia"))
  p <- ggplot(shannon_index_tree) + 
    geom_histogram(aes(shannon_index), binwidth=0.05) + 
    geom_vline(aes(xintercept=threshold), color="red", linetype="dashed") +
    facet_wrap(~case, ncol=4, scales="free_y") +
    ggtitle("Tree") +
    xlab("Shannon index") +
    ylab("Number of species")
  
  ggsave("data/output/Mammalia/shannon_index_tree_summary.pdf", plot=p, width=8, height=4)
  return("data/output/Mammalia/shannon_index_tree_summary.pdf")
}


# Visualize likelihoods to decide steps to be removed as burn-in period of MCMC
plot_likelihood <- function(bisse_result, taxon, removed_step) {
  burn_in <- tibble(threshold=c(removed_step), taxon=c(taxon))
  p <- ggplot(bisse_result, aes(x=i, y=p)) +
  ggtitle(taxon) +
    xlab("Step") +
    ylab("Log-likelihood") +
    geom_line() +
    geom_vline(data=burn_in, aes(xintercept=threshold), color="red", linetype="dashed")
  
  ggsave(paste0("data/output/", taxon, "/likelihood.pdf"), plot=p, width=8, height=3)
  return(paste0("data/output/", taxon, "/likelihood.pdf"))
}


summarize_likelihood <- function(bisse_result, taxon) {
  likelihood <- bisse_result %>% 
    dplyr::select(i, p) %>% 
    mutate(taxon = taxon) %>% 
    as_tibble()
  return(likelihood)
}


plot_likelihood_summary <- function(likelihood, removed_step) {
  likelihood$taxon <- factor(likelihood$taxon, levels=c("Mammalia", "Rodentia", "Primates", "Cetartiodactyla", "Carnivora", "Chiroptera", "Eulipotyphla", "Marsupialia"))
  burn_in <- tibble(threshold=rep(removed_step, 8), taxon=factor(c("Mammalia", "Rodentia", "Primates", "Cetartiodactyla", "Carnivora", "Chiroptera", "Eulipotyphla", "Marsupialia"), levels=c("Mammalia", "Rodentia", "Primates", "Cetartiodactyla", "Carnivora", "Chiroptera", "Eulipotyphla", "Marsupialia")))

  p <- ggplot(likelihood, aes(x=i, y=p)) +
      facet_wrap(~taxon, ncol=4, scales="free_y") +
      xlab("Step") +
      ylab("Log-likelihood") +
      geom_line() +
      geom_vline(data=burn_in, aes(xintercept=threshold), color="red", linetype="dashed")
  
  ggsave("data/output/Mammalia/likelihood_summary.pdf", plot=p, width=13, height=5)
  return("data/output/Mammalia/likelihood_summary.pdf")
}


# Visualize autocorrelation to decide sampling interval
plot_autocorrelation <- function(bisse_result, taxon, removed_step, max_interval, threshold) {
  # Remove steps as a burn-in period
  bisse_result <- bisse_result[(removed_step + 1):(nrow(bisse_result)-removed_step),]

  # Calculate autocorrelation of the likelihoods in each interval
  x <- bisse_result$p
  mean_x <- mean(x)
  autocorrelation <- c()
  for (k in 1:max_interval) {
    covariance <- 0
    mean_x1 <- mean(x[1:(length(x)-k)])
    mean_x2 <- mean(x[(1+k):length(x)])
    for (i in 1:(length(x) - k)) {
      covariance <- covariance + (x[i] - mean_x1) * (x[i+k] - mean_x2)
    }
    variance1 <- 0
    variance2 <- 0
    for (i in 1:(length(x) - k)) {
      variance1 <- variance1 + (x[i] - mean_x1) * (x[i] - mean_x1)
      variance2 <- variance2 + (x[i+k] - mean_x2) * (x[i+k] - mean_x2)
    }
    autocorrelation <- append(autocorrelation, covariance / (sqrt(variance1) * sqrt(variance2)))
  }

  interval <- tibble(threshold=c(threshold), taxon=c(taxon))
  p <- ggplot(tibble(autocorrelation), aes(seq_along(autocorrelation), autocorrelation)) +
    geom_bar(stat="identity") +
    ggtitle(taxon) +
    xlab("Sampling interval") +
    ylab("Autocorrelation") +
    geom_vline(data=interval, aes(xintercept=threshold), color="red", linetype="dashed")
  
  ggsave(paste0("data/output/", taxon, "/autocorrelation.pdf"), plot=p, width=3.5, height=3.5)
  return(paste0("data/output/", taxon, "/autocorrelation.pdf"))
}

summarize_autocorrelation <- function(bisse_result, taxon, removed_step, max_interval) {
  # Remove steps as a burn-in period
  bisse_result <- bisse_result[(removed_step + 1):(nrow(bisse_result)-removed_step),]

  # Calculate autocorrelation of the likelihoods in each interval
  x <- bisse_result$p
  mean_x <- mean(x)
  autocorrelation <- c()
  for (k in 1:max_interval) {
    covariance <- 0
    mean_x1 <- mean(x[1:(length(x)-k)])
    mean_x2 <- mean(x[(1+k):length(x)])
    for (i in 1:(length(x) - k)) {
      covariance <- covariance + (x[i] - mean_x1) * (x[i+k] - mean_x2)
    }
    variance1 <- 0
    variance2 <- 0
    for (i in 1:(length(x) - k)) {
      variance1 <- variance1 + (x[i] - mean_x1) * (x[i] - mean_x1)
      variance2 <- variance2 + (x[i+k] - mean_x2) * (x[i+k] - mean_x2)
    }
    autocorrelation <- append(autocorrelation, covariance / (sqrt(variance1) * sqrt(variance2)))
  }
  
  autocorrelation <- tibble(ac=autocorrelation)
  autocorrelation <- autocorrelation %>% 
    mutate(taxon = taxon) %>% 
    mutate(i = as.numeric(row.names(autocorrelation)))
  
  return(autocorrelation)
}

plot_autocorrelation_summary <- function(autocorrelation, threshold) {
  autocorrelation$taxon <- factor(autocorrelation$taxon, levels=c("Mammalia", "Rodentia", "Primates", "Cetartiodactyla", "Carnivora", "Chiroptera", "Eulipotyphla", "Marsupialia"))
  interval <- tibble(threshold=rep(threshold, 8), taxon=factor(c("Mammalia", "Rodentia", "Primates", "Cetartiodactyla", "Carnivora", "Chiroptera", "Eulipotyphla", "Marsupialia"), levels=c("Mammalia", "Rodentia", "Primates", "Cetartiodactyla", "Carnivora", "Chiroptera", "Eulipotyphla", "Marsupialia")))

  p <- ggplot(autocorrelation, aes(i, ac)) +
    facet_wrap(~taxon, ncol=4, scales="fixed") +
    xlab("Sampling interval") +
    ylab("Autocorrelation") +
    geom_bar(stat="identity") +
    geom_vline(data=interval, aes(xintercept=threshold), color="red", linetype="dashed")
  
  ggsave(paste0("data/output/Mammalia/autocorrelation_summary.pdf"), plot=p, width=10, height=5)
  return(paste0("data/output/Mammalia/autocorrelation_summary.pdf"))
}


# Visualize a tree with species names
plot_tree <- function(pruned_tree, character, redlist_species_in_tree, taxon) {
  character <- character %>% 
    filter(species %in% redlist_species_in_tree) %>%
    mutate(type = ifelse(generalist==0, "Specialist", "Generalist")) %>% 
    mutate(pos = 1) %>% 
    dplyr::select(species, type, pos)
  character$type <- factor(character$type, levels=c("Specialist", "Generalist"))
  character$species <- gsub(" ", "_", character$species)
  
  p <- ggtree(pruned_tree, layout="fan", open.angle=10, size=0.5) +
    geom_tiplab(size=0.7) + 
    geom_fruit(
      data=character,
      geom=geom_tile,
      mapping=aes(y=species, x=pos, fill=type),
      offset=0.006,
      pwidth=0.04
    ) + 
    scale_fill_manual(values=c("#004165", "#eaab00"))

  ggsave(paste0("data/output/", taxon, "/pruned_tree.pdf"), plot=p, width=40, height=40)
  return(paste0("data/output/", taxon, "/pruned_tree.pdf"))
}


# Visualize estimated rates of speciation, extinction, and transition in BiSSE model
plot_bisse <- function(bisse_result, taxon, burn_in, interval) {
  # Drop burn-in
  bisse_result <- bisse_result[burn_in+1:(nrow(bisse_result)-burn_in),]
  # Sampling
  samples <- filter(bisse_result, row_number()%%interval == 1)

  # Arrange dataframe
  speciation <- samples %>%
    dplyr::select(lambda0, lambda1) %>%
    rename("Specialist"=lambda0, "Generalist"=lambda1) %>%
    gather(key=Character, value=estimates)
  speciation$case <- rep("Speciation rate", length(speciation$estimates))

  extinction <- samples %>%
    dplyr::select(mu0, mu1) %>%
    rename("Specialist"=mu0, "Generalist"=mu1) %>%
    gather(key=Character, value=estimates)
  extinction$case <- rep("Extinction rate", length(extinction$estimates))

  transition <- samples %>%
    dplyr::select(q01, q10) %>%
    rename("Specialist"=q01, "Generalist"=q10) %>%
    gather(key=Character, value=estimates)
  transition$case <- rep("Transition rate", length(transition$estimates))

  cc <- speciation %>%
    bind_rows(extinction) %>%
    bind_rows(transition)
  
  # Mean
  cc_mean <- tibble(Character="Specialist", estimates=mean(samples$lambda0), case="Speciation rate") %>%
    bind_rows(tibble(Character="Generalist", estimates=mean(samples$lambda1), case="Speciation rate")) %>%
    bind_rows(tibble(Character="Specialist", estimates=mean(samples$mu0), case="Extinction rate")) %>%
    bind_rows(tibble(Character="Generalist", estimates=mean(samples$mu1), case="Extinction rate")) %>%
    bind_rows(tibble(Character="Specialist", estimates=mean(samples$q01), case="Transition rate")) %>%
    bind_rows(tibble(Character="Generalist", estimates=mean(samples$q10), case="Transition rate"))

  # Reorder
  cc$case <- factor(cc$case, levels=c("Speciation rate", "Transition rate", "Extinction rate"))
  cc$Character <- factor(cc$Character, levels=c("Specialist", "Generalist"))
  cc_mean$case <- factor(cc_mean$case, levels=c("Speciation rate", "Transition rate", "Extinction rate"))
  cc_mean$Character <- factor(cc_mean$Character, levels=c("Specialist", "Generalist"))

  # Rename columns
  cc <- cc %>% rename(State=Character)
  cc_mean <- cc_mean %>% rename(State=Character)

  # Plot
  p <- ggplot(cc) +
    geom_histogram(aes(estimates, fill=State), alpha=0.5, position="identity", binwidth=0.01) +
    geom_vline(data=cc_mean, aes(xintercept=estimates, color=State), linetype="dashed") +
    scale_fill_manual(values=c("#004165", "#eaab00")) +
    scale_color_manual(values=c("#004165", "#eaab00")) +
    facet_wrap(~case, ncol=3, scales="fixed") +
    ggtitle(str_to_title(taxon)) +
    xlab("Parameter estimate") +
    ylab("Number of estimates")
  
  # Save
  ggsave(paste0("data/output/", taxon, "/bisse.pdf"), plot=p, width=7, height=2.4)

  tibble(State="Specialist", case="Speciation rate", mean=mean(samples$lambda0), sd=sd(samples$lambda0)) %>%
    bind_rows(tibble(State="Generalist", case="Speciation rate", mean=mean(samples$lambda1), sd=sd(samples$lambda1))) %>%
    bind_rows(tibble(State="Specialist", case="Extinction rate", mean=mean(samples$mu0), sd=sd(samples$mu0))) %>%
    bind_rows(tibble(State="Generalist", case="Extinction rate", mean=mean(samples$mu1), sd=sd(samples$mu1))) %>%
    bind_rows(tibble(State="Specialist", case="Transition rate", mean=mean(samples$q01), sd=sd(samples$q01))) %>%
    bind_rows(tibble(State="Generalist", case="Transition rate", mean=mean(samples$q10), sd=sd(samples$q10))) %>%
    write_tsv(paste0("data/output/", taxon, "/bisse.tsv"))
  
  return(paste0("data/output/", taxon, "/bisse.tsv"))
}


# Visualize estimated rates of speciation and transition in BiSSE model
plot_speciation_transition <- function(bisse_result, taxon, burn_in, interval) {
  # Drop burn-in
  bisse_result <- bisse_result[burn_in+1:(nrow(bisse_result)-burn_in),]
  # Sampling
  samples <- filter(bisse_result, row_number()%%interval == 1)

  # Arrange dataframe
  speciation <- samples %>%
    dplyr::select(lambda0, lambda1) %>%
    rename("Specialist"=lambda0, "Generalist"=lambda1) %>%
    gather(key=Character, value=estimates)
  speciation$case <- rep("Speciation rate", length(speciation$estimates))

  transition <- samples %>%
    dplyr::select(q01, q10) %>%
    rename("Specialist"=q01, "Generalist"=q10) %>%
    gather(key=Character, value=estimates)
  transition$case <- rep("Transition rate", length(transition$estimates))

  cc <- speciation %>%
    bind_rows(transition)
  
  # Mean
  cc_mean <- tibble(Character="Specialist", estimates=mean(samples$lambda0), case="Speciation rate") %>%
    bind_rows(tibble(Character="Generalist", estimates=mean(samples$lambda1), case="Speciation rate")) %>%
    bind_rows(tibble(Character="Specialist", estimates=mean(samples$q01), case="Transition rate")) %>%
    bind_rows(tibble(Character="Generalist", estimates=mean(samples$q10), case="Transition rate"))

  # Reorder
  cc$case <- factor(cc$case, levels=c("Speciation rate", "Transition rate"))
  cc$Character <- factor(cc$Character, levels=c("Specialist", "Generalist"))
  cc_mean$case <- factor(cc_mean$case, levels=c("Speciation rate", "Transition rate"))
  cc_mean$Character <- factor(cc_mean$Character, levels=c("Specialist", "Generalist"))

  # Rename columns
  cc <- cc %>% rename(State=Character)
  cc_mean <- cc_mean %>% rename(State=Character)

  # Plot
  p <- ggplot(cc) +
    geom_histogram(aes(estimates, fill=State), alpha=0.5, position="identity", binwidth=0.01) +
    geom_vline(data=cc_mean, aes(xintercept=estimates, color=State), linetype="dashed") +
    scale_fill_manual(values=c("#004165", "#eaab00")) +
    scale_color_manual(values=c("#004165", "#eaab00")) +
    facet_wrap(~case, ncol=3, scales="fixed") +
    ggtitle(str_to_title(taxon)) +
    xlab("Parameter estimate") +
    ylab("Number of estimates")
  
  # Save
  ggsave(paste0("data/output/", taxon, "/bisse.pdf"), plot=p, width=5, height=2.4)

  tibble(State="Specialist", case="Speciation rate", mean=mean(samples$lambda0), sd=sd(samples$lambda0)) %>%
    bind_rows(tibble(State="Generalist", case="Speciation rate", mean=mean(samples$lambda1), sd=sd(samples$lambda1))) %>%
    bind_rows(tibble(State="Specialist", case="Transition rate", mean=mean(samples$q01), sd=sd(samples$q01))) %>%
    bind_rows(tibble(State="Generalist", case="Transition rate", mean=mean(samples$q10), sd=sd(samples$q10))) %>%
    write_tsv(paste0("data/output/", taxon, "/bisse.tsv"))
  
  return(paste0("data/output/", taxon, "/bisse.tsv"))
}

plot_bisse_diff <- function(bisse_result, taxon, burn_in, interval) {
  # Drop burn-in
  bisse_result <- bisse_result[burn_in+1:(nrow(bisse_result)-burn_in),]
  # Sampling
  samples <- filter(bisse_result, row_number()%%interval == 1)

  samples$speciation <- samples$lambda0 - samples$lambda1
  samples$transition <- samples$q01 - samples$q10
  samples$extinction <- samples$mu0 - samples$mu1

  # Arrange dataframe
  speciation <- samples %>%
    dplyr::select(speciation) %>%
    rename(estimates=speciation)
  speciation$case <- rep("Speciation rate", length(speciation$estimates))

  transition <- samples %>%
    dplyr::select(transition) %>%
    rename(estimates=transition)
  transition$case <- rep("Transition rate", length(transition$estimates))

  extinction <- samples %>%
    dplyr::select(extinction) %>%
    rename(estimates=extinction)
  extinction$case <- rep("Extinction rate", length(extinction$estimates))

  cc <- speciation %>%
    bind_rows(transition) %>%
    bind_rows(extinction)

  # Probability direction
  cc_pd <- tibble(p_direction=formatC(max(sum(speciation$estimates > 0) / length(speciation$estimates), sum(speciation$estimates < 0) / length(speciation$estimates)), digits=2, format="f"), case="Speciation rate") %>%
    bind_rows(tibble(p_direction=formatC(max(sum(transition$estimates > 0) / length(transition$estimates), sum(transition$estimates < 0) / length(transition$estimates)), digits=2, format="f"), case="Transition rate")) %>%
    bind_rows(tibble(p_direction=formatC(max(sum(extinction$estimates > 0) / length(extinction$estimates), sum(extinction$estimates < 0) / length(extinction$estimates)), digits=2, format="f"), case="Extinction rate"))

  # Reorder
  cc$case <- factor(cc$case, levels=c("Speciation rate", "Transition rate", "Extinction rate"))
  cc_pd$case <- factor(cc_pd$case, levels=c("Speciation rate", "Transition rate", "Extinction rate"))

  cc_pd$zero <- 0

  # Plot
  p <- ggplot(cc) +
    geom_histogram(aes(estimates), alpha=0.5, position="identity", binwidth=0.01) +
    geom_vline(data=cc_pd, aes(xintercept=zero)) +
    geom_label(data=cc_pd, aes(label=paste0("pd = ", p_direction)), x=Inf, y=Inf, hjust=1.05,vjust=1.2) +
    facet_wrap(~case, ncol=3, scales="fixed") +
    ggtitle(str_to_title(taxon)) +
    xlab("Parameter estimate (specialist - generalist)") +
    ylab("Number of estimates")
  
  # Save
  ggsave(paste0("data/output/", taxon, "/bisse_diff.pdf"), plot=p, width=7, height=2.4)
  return(paste0("data/output/", taxon, "/bisse_diff.pdf"))
}


# Visualize estimated speciation rates in BiSSE model
filter_mcmc_speciation <- function(bisse_result, taxa, burn_in, interval) {
  # Drop burn-in
  bisse_result <- bisse_result[burn_in+1:(nrow(bisse_result)-burn_in),]
  # Sampling
  samples <- filter(bisse_result, row_number()%%interval == 1)

  speciation <- samples %>%
    dplyr::select(lambda0, lambda1) %>%
    rename("Specialist"=lambda0, "Generalist"=lambda1) %>%
    gather(key=Character, value=estimates) %>% 
    mutate(taxon = taxa) %>% 
    as_tibble()
  return(speciation)
}


filter_mcmc_speciation_diff <- function(bisse_result, taxa, burn_in, interval) {
  # Drop burn-in
  bisse_result <- bisse_result[burn_in+1:(nrow(bisse_result)-burn_in),]
  # Sampling
  samples <- filter(bisse_result, row_number()%%interval == 1)

  samples$diff <- samples$lambda0 - samples$lambda1
  speciation <- samples %>%
    dplyr::select(diff) %>%
    mutate(taxon = taxa) %>% 
    as_tibble()
  return(speciation)
}


filter_mcmc_transition <- function(bisse_result, taxa, burn_in, interval) {
  # Drop burn-in
  bisse_result <- bisse_result[burn_in+1:(nrow(bisse_result)-burn_in),]
  # Sampling
  samples <- filter(bisse_result, row_number()%%interval == 1)

  transition <- samples %>%
    dplyr::select(q01, q10) %>%
    rename("Specialist"=q01, "Generalist"=q10) %>%
    gather(key=Character, value=estimates) %>% 
    mutate(taxon = taxa) %>% 
    as_tibble()
  return(transition)
}


filter_mcmc_transition_diff <- function(bisse_result, taxa, burn_in, interval) {
  # Drop burn-in
  bisse_result <- bisse_result[burn_in+1:(nrow(bisse_result)-burn_in),]
  # Sampling
  samples <- filter(bisse_result, row_number()%%interval == 1)

  samples$diff <- samples$q01 - samples$q10
  transition <- samples %>%
    dplyr::select(diff) %>%
    mutate(taxon = taxa) %>% 
    as_tibble()
  return(transition)
}


filter_mcmc_extinction <- function(bisse_result, taxa, burn_in, interval) {
  # Drop burn-in
  bisse_result <- bisse_result[burn_in+1:(nrow(bisse_result)-burn_in),]
  # Sampling
  samples <- filter(bisse_result, row_number()%%interval == 1)

  extinction <- samples %>%
    dplyr::select(mu0, mu1) %>%
    rename("Specialist"=mu0, "Generalist"=mu1) %>%
    gather(key=Character, value=estimates) %>% 
    mutate(taxon = taxa) %>% 
    as_tibble()
  return(extinction)
}


filter_mcmc_extinction_diff <- function(bisse_result, taxa, burn_in, interval) {
  # Drop burn-in
  bisse_result <- bisse_result[burn_in+1:(nrow(bisse_result)-burn_in),]
  # Sampling
  samples <- filter(bisse_result, row_number()%%interval == 1)

  samples$diff <- samples$mu0 - samples$mu1
  extinction <- samples %>%
    dplyr::select(diff) %>%
    mutate(taxon = taxa) %>% 
    as_tibble()
  return(extinction)
}


summarize_samples <- function(filtered_samples, taxa) {
  summary <- filtered_samples %>% 
    group_by(Character, taxon) %>%
    summarise_all(mean) %>% 
    ungroup()
  return(summary)
}


summarize_samples_diff <- function(filtered_samples, taxa) {
  summary <- filtered_samples %>% 
    group_by(taxon) %>%
    summarise_all(mean) %>% 
    ungroup()
  summary$p_direction <- formatC(max(sum(filtered_samples$diff > 0) / length(filtered_samples$diff), sum(filtered_samples$diff < 0) / length(filtered_samples$diff)), digits=2, format="f")
  return(summary)
}


plot_speciation <- function(speciation, speciation_summary) {
  speciation$Character <- factor(speciation$Character, levels=c("Specialist", "Generalist"))
  speciation_summary$Character <- factor(speciation_summary$Character, levels=c("Specialist", "Generalist"))
  speciation$taxon <- factor(speciation$taxon, levels=c("Mammalia", "Rodentia", "Primates", "Cetartiodactyla", "Carnivora", "Chiroptera", "Eulipotyphla", "Marsupialia"))
  speciation_summary$taxon <- factor(speciation_summary$taxon, levels=c("Mammalia", "Rodentia", "Primates", "Cetartiodactyla", "Carnivora", "Chiroptera", "Eulipotyphla", "Marsupialia"))
  speciation <- speciation %>% rename(State=Character)
  speciation_summary <- speciation_summary %>% rename(State=Character)

  p <- ggplot(speciation) +
    geom_histogram(aes(estimates, fill=State), alpha=0.5, position="identity", binwidth=0.01) +
    geom_vline(data=speciation_summary, aes(xintercept=estimates, color=State), linetype="dashed") +
    scale_fill_manual(values=c("#004165", "#eaab00")) +
    scale_color_manual(values=c("#004165", "#eaab00")) +
    facet_wrap(~taxon, ncol=4, scales="fixed") +
    xlab("Speciation rate") +
    ylab("Number of estimates")
  
  ggsave("data/output/Mammalia/speciation.pdf", plot=p, width=9.5, height=4)
  return("data/output/Mammalia/speciation.pdf")
}


plot_speciation_typical <- function(speciation, speciation_summary) {
  speciation$Character <- factor(speciation$Character, levels=c("Specialist", "Generalist"))
  speciation_summary$Character <- factor(speciation_summary$Character, levels=c("Specialist", "Generalist"))
  speciation$taxon <- factor(speciation$taxon, levels=c("Rodentia", "Primates", "Cetartiodactyla", "Carnivora"))
  speciation_summary$taxon <- factor(speciation_summary$taxon, levels=c("Rodentia", "Primates", "Cetartiodactyla", "Carnivora"))
  speciation <- speciation %>% rename(State=Character)
  speciation_summary <- speciation_summary %>% rename(State=Character)

  p <- ggplot(speciation) +
    geom_histogram(aes(estimates, fill=State), alpha=0.5, position="identity", binwidth=0.01) +
    geom_vline(data=speciation_summary, aes(xintercept=estimates, color=State), linetype="dashed") +
    scale_fill_manual(values=c("#004165", "#eaab00")) +
    scale_color_manual(values=c("#004165", "#eaab00")) +
    facet_wrap(~taxon, ncol=4, scales="fixed") +
    xlab("Speciation rate") +
    ylab("Number of estimates")
  
  ggsave("data/output/Mammalia/speciation_typical.pdf", plot=p, width=9.5, height=2.4)
  return("data/output/Mammalia/speciation_typical.pdf")
}


plot_speciation_atypical <- function(speciation, speciation_summary) {
  speciation$Character <- factor(speciation$Character, levels=c("Specialist", "Generalist"))
  speciation_summary$Character <- factor(speciation_summary$Character, levels=c("Specialist", "Generalist"))
  speciation$taxon <- factor(speciation$taxon, levels=c("Chiroptera", "Eulipotyphla", "Marsupialia", "Aves"))
  speciation_summary$taxon <- factor(speciation_summary$taxon, levels=c("Chiroptera", "Eulipotyphla", "Marsupialia", "Aves"))
  speciation <- speciation %>% rename(State=Character)
  speciation_summary <- speciation_summary %>% rename(State=Character)

  p <- ggplot(speciation) +
    geom_histogram(aes(estimates, fill=State), alpha=0.5, position="identity", binwidth=0.01) +
    geom_vline(data=speciation_summary, aes(xintercept=estimates, color=State), linetype="dashed") +
    scale_fill_manual(values=c("#004165", "#eaab00")) +
    scale_color_manual(values=c("#004165", "#eaab00")) +
    facet_wrap(~taxon, ncol=4, scales="fixed") +
    xlab("Speciation rate") +
    ylab("Number of estimates")
  
  ggsave("data/output/Mammalia/speciation_atypical.pdf", plot=p, width=9.5, height=2.4)
  return("data/output/Mammalia/speciation_atypical.pdf")
}


plot_speciation_diff <- function(speciation, speciation_summary) {
  speciation$taxon <- factor(speciation$taxon, levels=c("Mammalia", "Rodentia", "Primates", "Cetartiodactyla", "Carnivora", "Chiroptera", "Eulipotyphla", "Marsupialia"))
  speciation_summary$taxon <- factor(speciation_summary$taxon, levels=c("Mammalia", "Rodentia", "Primates", "Cetartiodactyla", "Carnivora", "Chiroptera", "Eulipotyphla", "Marsupialia"))
  speciation_summary$zero <- 0

  p <- ggplot(speciation) +
    geom_histogram(aes(diff), alpha=0.5, position="identity", binwidth=0.01) +
    # geom_vline(data=speciation_summary, aes(xintercept=diff), linetype="dashed") +
    geom_vline(data=speciation_summary, aes(xintercept=zero)) +
    geom_label(data=speciation_summary, aes(label=paste0("pd = ", p_direction)), x=Inf, y=Inf, hjust=1.05,vjust=1.2) +
    facet_wrap(~taxon, ncol=4, scales="fixed") +
    xlab("(Speciation rate of specialist) - (Speciation rate of generalist)") +
    ylab("Number of estimates")
  
  ggsave("data/output/Mammalia/speciation_diff.pdf", plot=p, width=9.5, height=4)
  return("data/output/Mammalia/speciation_diff.pdf")
}


plot_transition <- function(transition, transition_summary) {
  transition$Character <- factor(transition$Character, levels=c("Specialist", "Generalist"))
  transition_summary$Character <- factor(transition_summary$Character, levels=c("Specialist", "Generalist"))
  transition$taxon <- factor(transition$taxon, levels=c("Mammalia", "Rodentia", "Primates", "Cetartiodactyla", "Carnivora", "Chiroptera", "Eulipotyphla", "Marsupialia"))
  transition_summary$taxon <- factor(transition_summary$taxon, levels=c("Mammalia", "Rodentia", "Primates", "Cetartiodactyla", "Carnivora", "Chiroptera", "Eulipotyphla", "Marsupialia"))
  transition <- transition %>% rename(State=Character)
  transition_summary <- transition_summary %>% rename(State=Character)

  p <- ggplot(transition) +
    geom_histogram(aes(estimates, fill=State), alpha=0.5, position="identity", binwidth=0.01) +
    geom_vline(data=transition_summary, aes(xintercept=estimates, color=State), linetype="dashed") +
    scale_fill_manual(values=c("#004165", "#eaab00")) +
    scale_color_manual(values=c("#004165", "#eaab00")) +
    facet_wrap(~taxon, ncol=4, scales="fixed") +
    xlab("Transition rate") +
    ylab("Number of estimates")
  
  ggsave("data/output/Mammalia/transition.pdf", plot=p, width=9.5, height=4)
  return("data/output/Mammalia/transition.pdf")
}


plot_transition_typical <- function(transition, transition_summary) {
  transition$Character <- factor(transition$Character, levels=c("Specialist", "Generalist"))
  transition_summary$Character <- factor(transition_summary$Character, levels=c("Specialist", "Generalist"))
  transition$taxon <- factor(transition$taxon, levels=c("Rodentia", "Primates", "Cetartiodactyla", "Carnivora"))
  transition_summary$taxon <- factor(transition_summary$taxon, levels=c("Rodentia", "Primates", "Cetartiodactyla", "Carnivora"))
  transition <- transition %>% rename(State=Character)
  transition_summary <- transition_summary %>% rename(State=Character)

  p <- ggplot(transition) +
    geom_histogram(aes(estimates, fill=State), alpha=0.5, position="identity", binwidth=0.01) +
    geom_vline(data=transition_summary, aes(xintercept=estimates, color=State), linetype="dashed") +
    scale_fill_manual(values=c("#004165", "#eaab00")) +
    scale_color_manual(values=c("#004165", "#eaab00")) +
    facet_wrap(~taxon, ncol=4, scales="fixed") +
    xlab("Transition rate") +
    ylab("Number of estimates")
  
  ggsave("data/output/Mammalia/transition_typical.pdf", plot=p, width=9.5, height=2.4)
  return("data/output/Mammalia/transition_typical.pdf")
}


plot_transition_atypical <- function(transition, transition_summary) {
  transition$Character <- factor(transition$Character, levels=c("Specialist", "Generalist"))
  transition_summary$Character <- factor(transition_summary$Character, levels=c("Specialist", "Generalist"))
  transition$taxon <- factor(transition$taxon, levels=c("Chiroptera", "Eulipotyphla", "Marsupialia", "Aves"))
  transition_summary$taxon <- factor(transition_summary$taxon, levels=c("Chiroptera", "Eulipotyphla", "Marsupialia", "Aves"))
  transition <- transition %>% rename(State=Character)
  transition_summary <- transition_summary %>% rename(State=Character)

  p <- ggplot(transition) +
    geom_histogram(aes(estimates, fill=State), alpha=0.5, position="identity", binwidth=0.01) +
    geom_vline(data=transition_summary, aes(xintercept=estimates, color=State), linetype="dashed") +
    scale_fill_manual(values=c("#004165", "#eaab00")) +
    scale_color_manual(values=c("#004165", "#eaab00")) +
    facet_wrap(~taxon, ncol=4, scales="fixed") +
    xlab("Transition rate") +
    ylab("Number of estimates")
  
  ggsave("data/output/Mammalia/transition_atypical.pdf", plot=p, width=9.5, height=2.4)
  return("data/output/Mammalia/transition_atypical.pdf")
}


plot_transition_diff <- function(transition, transition_summary) {
  transition$taxon <- factor(transition$taxon, levels=c("Mammalia", "Rodentia", "Primates", "Cetartiodactyla", "Carnivora", "Chiroptera", "Eulipotyphla", "Marsupialia"))
  transition_summary$taxon <- factor(transition_summary$taxon, levels=c("Mammalia", "Rodentia", "Primates", "Cetartiodactyla", "Carnivora", "Chiroptera", "Eulipotyphla", "Marsupialia"))
  transition_summary$zero <- 0

  p <- ggplot(transition) +
    geom_histogram(aes(diff), alpha=0.5, position="identity", binwidth=0.01) +
    # geom_vline(data=transition_summary, aes(xintercept=diff), linetype="dashed") +
    geom_vline(data=transition_summary, aes(xintercept=zero)) +
    geom_label(data=transition_summary, aes(label=paste0("pd = ", p_direction)), x=Inf, y=Inf, hjust=1.05,vjust=1.2) +
    facet_wrap(~taxon, ncol=4, scales="fixed") +
    xlab("(Transition rate from specialist to generalist) - (Transition rate from generalist to specialist)") +
    ylab("Number of estimates")
  
  ggsave("data/output/Mammalia/transition_diff.pdf", plot=p, width=9.5, height=4)
  return("data/output/Mammalia/transition_diff.pdf")
}


plot_extinction <- function(extinction, extinction_summary) {
  extinction$Character <- factor(extinction$Character, levels=c("Specialist", "Generalist"))
  extinction_summary$Character <- factor(extinction_summary$Character, levels=c("Specialist", "Generalist"))
  extinction$taxon <- factor(extinction$taxon, levels=c("Mammalia", "Rodentia", "Primates", "Cetartiodactyla", "Carnivora", "Chiroptera", "Eulipotyphla", "Marsupialia"))
  extinction_summary$taxon <- factor(extinction_summary$taxon, levels=c("Mammalia", "Rodentia", "Primates", "Cetartiodactyla", "Carnivora", "Chiroptera", "Eulipotyphla", "Marsupialia"))
  extinction <- extinction %>% rename(State=Character)
  extinction_summary <- extinction_summary %>% rename(State=Character)

  p <- ggplot(extinction) +
    geom_histogram(aes(estimates, fill=State), alpha=0.5, position="identity", binwidth=0.01) +
    geom_vline(data=extinction_summary, aes(xintercept=estimates, color=State), linetype="dashed") +
    scale_fill_manual(values=c("#004165", "#eaab00")) +
    scale_color_manual(values=c("#004165", "#eaab00")) +
    facet_wrap(~taxon, ncol=4, scales="fixed") +
    xlab("Extinction rate") +
    ylab("Number of estimates")
  
  ggsave("data/output/Mammalia/extinction.pdf", plot=p, width=9.5, height=4)
  return("data/output/Mammalia/extinction.pdf")
}


plot_extinction_diff <- function(extinction, extinction_summary) {
  extinction$taxon <- factor(extinction$taxon, levels=c("Mammalia", "Rodentia", "Primates", "Cetartiodactyla", "Carnivora", "Chiroptera", "Eulipotyphla", "Marsupialia"))
  extinction_summary$taxon <- factor(extinction_summary$taxon, levels=c("Mammalia", "Rodentia", "Primates", "Cetartiodactyla", "Carnivora", "Chiroptera", "Eulipotyphla", "Marsupialia"))
  extinction_summary$zero <- 0

  p <- ggplot(extinction) +
    geom_histogram(aes(diff), alpha=0.5, position="identity", binwidth=0.01) +
    # geom_vline(data=extinction_summary, aes(xintercept=diff), linetype="dashed") +
    geom_vline(data=extinction_summary, aes(xintercept=zero)) +
    geom_label(data=extinction_summary, aes(label=paste0("pd = ", p_direction)), x=-Inf, y=Inf, hjust=-0.05,vjust=1.2) +
    facet_wrap(~taxon, ncol=4, scales="fixed") +
    xlab("(Extinction rate of specialist) - (Extinction rate of generalist)") +
    ylab("Number of estimates")
  
  ggsave("data/output/Mammalia/extinction_diff.pdf", plot=p, width=9.5, height=4)
  return("data/output/Mammalia/extinction_diff.pdf")
}


summarize_number <- function(character, redlist_species, redlist_species_in_tree, taxon) {
  species_all <- character %>% 
    filter(species %in% redlist_species)
  species_tree <- character %>% 
    filter(species %in% redlist_species_in_tree)
  
  specialist_all <- species_all %>% 
    filter(generalist == 0) %>% 
    .$species %>% 
    length()
  generalist_all <- length(species_all$species) - specialist_all
  specialist_tree <- species_tree %>% 
    filter(generalist == 0) %>% 
    .$species %>% 
    length()
  generalist_tree <- length(species_tree$species) - specialist_tree

  sampling_fraction_specialist <- specialist_tree / specialist_all
  sampling_fraction_generalist <- generalist_tree / generalist_all

  df <- tibble(Taxon = taxon, Specialist_all = specialist_all, Generalist_all = generalist_all, Specialist_tree = specialist_tree, Generalist_tree = generalist_tree, Specialist_sampling.f = sampling_fraction_specialist, Generalist_sampling.f = sampling_fraction_generalist)
  return(df)
}


output_summarized_number <- function(df) {
  df$Taxon <- factor(df$Taxon, levels=c("Mammalia", "Rodentia", "Primates", "Cetartiodactyla", "Carnivora", "Chiroptera", "Eulipotyphla", "Marsupialia"))
  df %>% 
    write_tsv("data/output/Mammalia/number_summary.tsv")
}



# --- Aves

preprocessing_redlist_Aves <- function(redlist_path) {
  # Read polygon data
  redlist <- st_read(redlist_path)

  # Resolve MULTISURFACE
  redlist <- st_cast(redlist, "MULTIPOLYGON")

  # Use "Extant", "Possibly Extinct", and "Extinct" distribution
  redlist <- filter(redlist, presence==1 | presence==4 | presence==5)

  # Use only "Native" distribution
  redlist <- filter(redlist, origin==1)

  return(redlist)
}


# Merge IUCN Red List polygon data with Köppen-Geiger climate classification map
merge_redlist_with_climate_Aves <- function(redlist, climate_map_path, confidence_map_path) {
  # Read raster files
  climate_map <- raster(climate_map_path)
  confidence_map <- raster(confidence_map_path)

  # Use only 100% confidence cells in the map
  kg_conf_binary <- calc(confidence_map, fun=function(x) {x==100})
  kg_100conf <- climate_map * kg_conf_binary

  # Layerize the map for parallelization in exactextractr
  kg_brick <- layerize(kg_100conf)

  # Merge polygon data with climate classes
  redlist <- cbind(redlist, exact_extract(kg_brick, redlist, 'weighted_sum', weights = area(kg_brick), stack_apply = TRUE))
  redlist <- st_drop_geometry(redlist)

  # Change colnames and select some columns
  tmp <- data.frame(binomial="temporary", weighted_sum.X0=0 , weighted_sum.X1=0, weighted_sum.X2=0, weighted_sum.X3=0, weighted_sum.X4=0, weighted_sum.X5=0, weighted_sum.X6=0, weighted_sum.X7=0, weighted_sum.X8=0, weighted_sum.X9=0, weighted_sum.X10=0, weighted_sum.X11=0, weighted_sum.X12=0, weighted_sum.X13=0, weighted_sum.X14=0, weighted_sum.X15=0, weighted_sum.X16=0, weighted_sum.X17=0, weighted_sum.X18=0, weighted_sum.X19=0, weighted_sum.X20=0, weighted_sum.X21=0, weighted_sum.X22=0, weighted_sum.X23=0, weighted_sum.X24=0, weighted_sum.X25=0, weighted_sum.X26=0, weighted_sum.X27=0, weighted_sum.X28=0, weighted_sum.X29=0, weighted_sum.X30=0)
  redlist <- merge(redlist, tmp, all=T) %>%
    as_tibble() %>%
    filter(binomial != "temporary") %>%
    mutate_if(is.numeric, ~tidyr::replace_na(., 0))
  redlist <- redlist %>%
    dplyr::rename(N = weighted_sum.X0, Af = weighted_sum.X1, Am = weighted_sum.X2, Aw = weighted_sum.X3, BWh = weighted_sum.X4, BWk = weighted_sum.X5, BSh = weighted_sum.X6, BSk = weighted_sum.X7, Csa = weighted_sum.X8, Csb = weighted_sum.X9, Csc = weighted_sum.X10, Cwa = weighted_sum.X11, Cwb = weighted_sum.X12, Cwc = weighted_sum.X13, Cfa = weighted_sum.X14, Cfb = weighted_sum.X15, Cfc = weighted_sum.X16, Dsa = weighted_sum.X17, Dsb = weighted_sum.X18, Dsc = weighted_sum.X19, Dsd = weighted_sum.X20, Dwa = weighted_sum.X21, Dwb = weighted_sum.X22, Dwc = weighted_sum.X23, Dwd = weighted_sum.X24, Dfa = weighted_sum.X25, Dfb = weighted_sum.X26, Dfc = weighted_sum.X27, Dfd = weighted_sum.X28, ET = weighted_sum.X29, EF = weighted_sum.X30)
  kg_list <- c("N", "Af", "Am", "Aw", "BWh", "BWk", "BSh", "BSk", "Csa", "Csb", "Csc", "Cwa", "Cwb", "Cwc", "Cfa", "Cfb", "Cfc", "Dsa", "Dsb", "Dsc", "Dsd", "Dwa", "Dwb", "Dwc", "Dwd", "Dfa", "Dfb", "Dfc", "Dfd", "ET", "EF")
  redlist$species <- redlist$binomial
  redlist$order <- ""
  redlist <- dplyr::select(redlist, c("species", all_of(kg_list)), order)

  return(redlist)
}


extract_species_from_redlist_Aves <- function(redlist) {
  return(redlist$species)
}


extract_species_from_tree_Aves <- function(tree_path) {
  tree <- ape::read.nexus(tree_path)
  return(gsub("_", " ", tree$tip.label))
}


# Visualize a tree with species names
plot_tree_Aves <- function(pruned_tree, character, redlist_species_in_tree, taxon) {
  character <- character %>% 
    filter(species %in% redlist_species_in_tree) %>%
    mutate(type = ifelse(generalist==0, "Specialist", "Generalist")) %>% 
    mutate(pos = 1) %>% 
    dplyr::select(species, type, pos)
  character$type <- factor(character$type, levels=c("Specialist", "Generalist"))
  character$species <- gsub(" ", "_", character$species)
  
  p <- ggtree(pruned_tree, layout="fan", open.angle=10, size=0.5) +
    geom_tiplab(size=0.6) + 
    geom_fruit(
      data=character,
      geom=geom_tile,
      mapping=aes(y=species, x=pos, fill=type),
      offset=0.025,
      pwidth=0.0005
    ) + 
    scale_fill_manual(values=c("#004165", "#eaab00"))

  ggsave(paste0("data/output/", taxon, "/pruned_tree.pdf"), plot=p, width=70, height=70, limitsize=FALSE)
  return(paste0("data/output/", taxon, "/pruned_tree.pdf"))
}


summarize_number_Aves <- function(character, redlist_species, redlist_species_in_tree, taxon) {
  species_all <- character %>% 
    filter(species %in% redlist_species)
  species_tree <- character %>% 
    filter(species %in% redlist_species_in_tree)
  
  specialist_all <- species_all %>% 
    filter(generalist == 0) %>% 
    .$species %>% 
    length()
  generalist_all <- length(species_all$species) - specialist_all
  specialist_tree <- species_tree %>% 
    filter(generalist == 0) %>% 
    .$species %>% 
    length()
  generalist_tree <- length(species_tree$species) - specialist_tree

  sampling_fraction_specialist <- specialist_tree / specialist_all
  sampling_fraction_generalist <- generalist_tree / generalist_all

  df <- tibble(Taxon = taxon, Specialist_all = specialist_all, Generalist_all = generalist_all, Specialist_tree = specialist_tree, Generalist_tree = generalist_tree, Specialist_sampling.f = sampling_fraction_specialist, Generalist_sampling.f = sampling_fraction_generalist)
  df %>% 
    write_tsv("data/output/Aves/number_summary.tsv")
}