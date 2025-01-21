### Calculate species diet niche breadth
### Use EltonTraits diet categories

library(tidyverse)
library(purrr)

# Raw data directory
data_dir <- "raw_data/"

# output directory
output_dir <- "derived_data/"

## EltonTraits 1.0 diet data with AOUs
bird_traits <- read.csv(paste0(data_dir, "breeding_bird_diet_mig_forage_traits.csv"), stringsAsFactors = F)

## Species list 
species_list <- read.csv(paste0(data_dir, "species_list_mbbs.csv"), stringsAsFactors = F)

## Filter to species of interest
bird_diets <- bird_traits %>%
  filter(aou %in% species_list$aou)

## Total possible diet categories
max_diet <- 10

## Shannon evenness of diet categories
diet_niche_breadth <- bird_diets %>%
  select(genus, species, english_common_name, aou, binomial, Diet.Inv:Diet.PlantO) %>%
  pivot_longer(names_to = "diet_cat", values_to = "prop", Diet.Inv:Diet.PlantO) %>%
  group_by(genus, species, english_common_name, aou, binomial) %>%
  summarize(shannonE_diet = -sum( (prop/100) * log(prop/100) / log(max_diet), na.rm = T))

write.csv(diet_niche_breadth, paste0(output_dir, "diet_niche_breadth_mbbs.csv"), row.names = F)
