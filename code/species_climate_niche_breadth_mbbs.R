## Climate niche breadth
## Hypervolume of species breeding range based on BIOCLIM variables

library(tidyverse)
library(purrr)
require(raster)
require(maps)
library(sf)
library(hypervolume)
library(rasterDT)
#library(geodata)

## BioArk directory
info <- sessionInfo()
bioark <- ifelse(grepl("apple", info$platform), "/Volumes", "\\\\BioArk")

## species range maps directory
range_dir <- paste0(bioark, "/hurlbertlab/GIS/birds/All/All/")
#for ijbg mbbs purposes
range_dir <- "Z:/GIS/birds/All/All/"

## output directory
output_dir <- "derived_data/"

## Species list 
species_list <- read.csv("raw_data/species_list_mbbs.csv", stringsAsFactors = F)
#testing, cut down
species_list <- species_list %>% filter(!english_common_name %in% c("Brown-headed Nuthatch", "Yellow-throated Warbler", "Northern Cardinal", "Yellow-breasted Chat","House Finch"))

# Match BBS taxonomy with breeding range polygon taxonomy
fix_mismatch <- read.csv("derived_data/fix_breedingrange_genus_mismatch.csv", stringsAsFactors = F) %>%
  dplyr::select(-file, -spp_name) %>%
  filter(!is.na(old_genus)) %>%
  mutate(new_binomial = paste(old_genus, species))

bbs_spp <- species_list %>%
  mutate(binomial = paste(genus, species, sep = " ")) %>%
  filter(!grepl("unid.", english_common_name), !grepl("hybrid", english_common_name)) %>%
  mutate_at(c("binomial"), ~case_when(grepl("Colaptes auratus", .) ~ "Colaptes auratus",
                                      grepl("Junco hyemalis", .) ~ "Junco hyemalis",
                                      grepl("Setophaga coronata", .) ~ "Dendroica coronata",
                                      TRUE ~ .)) %>%
  left_join(fix_mismatch) %>%
  mutate(matched_name = ifelse(!is.na(old_genus), new_binomial, binomial),
         matched_filename = gsub(" ", "_", matched_name)) %>%
  filter(!is.na(matched_name))

range_files <- data.frame(file = list.files(range_dir)) %>%
  filter(grepl(".shp", file)) %>%
  mutate(spp_name = word(file, 1, 2, sep = "_"),
         file_binomial = gsub("_", " ", spp_name)) 

spp_list <- range_files %>%
  right_join(bbs_spp, by = c("file_binomial" = "matched_name")) %>%
  #add new columns for comparison
  mutate(climate_vol_1.4 = NA,
         climate_vol_2.1 = NA)

# WorldClim
#wc1.4_10m_bio <- list.files("raw_data/wc1.4_10m_bio/", pattern = ".bil$", full.names = TRUE)
wc2.1_10m_bio <- list.files("raw_data/wc2.1_10m_bio/", pattern = ".tif$", full.names = TRUE)
#wc2.1_10m_bio_min5 <- wc2.1_10m_bio[1:2]

climatelayerslist <- list(wc2.1_10m_bio)

for(wc in 1:length(climatelayerslist)) {
  
  climatelayers <- raster::stack(climatelayerslist[[wc]])
  
  #already trimmed so the datasource only has layers 5, 10, and 18.
  climatelayers_ss <- climatelayers
  
  nlayers(climatelayers_ss) #yep, prints 3 as expected.
  
  #Worldclim crs
  bio_crs <- st_crs(climatelayers_ss)
  
  # z-transform climate layers to make axes comparable
  for (i in 1:nlayers(climatelayers_ss)) {
    climatelayers_ss[[i]] <- (climatelayers_ss[[i]] - cellStats(climatelayers_ss[[i]], 'mean')) / cellStats(climatelayers_ss[[i]], 'sd') 
  }
  print("z-tranformed")
  
  #crop extent
  climatelayers_ss_cropped = crop(climatelayers_ss, extent(-150,-50,15,60))
  
  #Okay, now let's bring in the bird data.
  for(a in 31:nrow(spp_list)){ 
    
    species <- spp_list$file[a]
    print(species)
    
    climate_vol <- climate_hypervolume(species)
    
    #assign climate vol to correct column of df
    if(wc == 1){
      spp_list$climate_vol_1.4[a] <- climate_vol
    } else if (wc == 2){
      spp_list$climate_vol_2.1[a] <- climate_vol
    }
    
    #and actually! let's save ALONG THE WAY
    write.csv(spp_list,
              paste0(output_dir, "2.1", spp_list$spp_name[a], ".csv"),
              row.names = F)
    
  } # end bird data loop
}


### Climate hypervolume function
### Input: species file path for range shapefile
### Output: climate hypervolume

climate_hypervolume <- function(species) {
  ## For each species: read in breeding range shapefile
  
  br <- read_sf(paste0(range_dir, species))
  
  # Use extant (PRESENCE 1-3) breeding and resident ranges (SEASONAL 1-2) only
  # Re-project to land cover CRS
  breeding_range <- br %>%
    filter(PRESENCE %in% c(1:3), SEASONAL %in% c(1:2)) %>%
    st_transform(bio_crs)
  
  ## Rasterize shapefile
  # fasterizeDT()
  
  null_rast <- raster(extent(br), res = res(climatelayers_ss_cropped), crs = bio_crs)
  
  br_rast <- fasterizeDT(breeding_range, null_rast)
  
  ## Range raster to df of coordinates
  br_coords <- rasterToPoints(br_rast)
  
  br_df <- data.frame(lon = br_coords[, 1], lat = br_coords[, 2])
  
  ## Extract clim vars
  br_clim <- extract(climatelayers_ss_cropped, br_df)
  
  br_nona <- na.omit(br_clim)
  
  if(nrow(br_nona) == 0) {
    return(NA) # Breeding range not in North America
  } else {
    ## Calculate hypervolume
    br_hyper <- hypervolume_gaussian(br_nona)
    vol <- get_volume(br_hyper)
    
    print(paste(Sys.time(), species, "complete"))
    
    return(vol)
  }
  

}


