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

#climatelayers <- getData('worldclim', var='bio', res=10, path=tempdir())
#climatelayers_ss = climatelayers[[c(5, 10, 18)]]
#WorldClim doesn't work and can't change the URL (which doesn't exist now) that it's pulling from. Shift to using geodata download, but geodata downloads only worldclim 2_1, and doesn't support previous versions. Which is causing problems with matching data.
#OH. The paper Grace links uses worldclim 2.1 anyway. huh, soo. don't need to use 1.4 at all? Okay. so. let's just try switching the res to 1k/30s and see if that doesn't fix this mismatch. Ignore using 1.4
#climatelayers <- geodata::worldclim_global(var = "bio", res = 0.5, path = tempdir())
#welp, no worries, changing the URL didn't fix it. That's okay. Instead, we'll just load it in. Downloaded from: url = "https://geodata.ucdavis.edu/climate/worldclim/1_4/grid/cur/bio_10m_bil.zip"
#1.4 download, load in
  #climatelayers <- list.files("raw_data/bio_10m_bil/", pattern = ".bil$", full.names = TRUE)
  climatelayers <- list.files("raw_data/wc2.1_30s/", pattern = ".tif$", full.names = TRUE)
  climatelayers <- stack(climatelayers)
climatelayers_ss = climatelayers[[c(5, 10, 18)]] #warm quarter max, mean, and precipitation
#so, previous download used version 1.4 (cannot open URL 'https://biogeo.ucdavis.edu/data/climate/worldclim/1_4/grid/cur/bio_10m_bil.zip') and geodata supports only 2.1 at the moment. We will just re-calc all these scores for all the birds and take a look at how they match up with prev.
#New download link is https://geodata.ucdavis.edu/climate/worldclim/1_4/grid/cur/bio_10m_bil.zip
#Okay so Northern Cardinal doesn't match. With Grace's original work it was caluclated at 1.099151035 and with the new version of worldclim it's calculated at 0.374591968
#soooo, okay now it's downloading 1.4 but it's still NOT getting a match. Now cardinal hypervolume is being calculated as 4.649983191.
#Is the problem the layers I'm trimming? Are those not the right layers? (5,10,18) - check that they layers haven't moved and confirm from Grace's paper what layer designations she INTENDED to use. Nope, these are the right ones.
#"To measure climate niche breadth, we calculated a hypervolume for each species using warmest quarter temperature and precipitation variables at 1-km resolution from BioClim"
#Ah! At the 1-km resolution! I think I've accidentally downloaded from the 10 minute.
#Okay so, not exactlyyyyy 1km, but I think 30 seconds is what she must have meant.
#But the specs here are wrong, because it is trying to download the 10m version from Worldclim. This info abt how to do this comes from the hypervolume package


#At the LEAST, I need to use the 1.4 to confirm I can replicate the Cardinal calculation. From there, if they match, we can see about if I then switch it to the new worldclim that includes an extra 30 years of climate data or not. 
#LOL the geodata server is down for maintence and expected back by Dec 21. pft.
#AH. this may ALSO be why the above link is not working, because the server is down.
#Yes :) Once the server came back up the path worked perfectly.

# WorldClim crs
bio_crs <- st_crs(climatelayers)

# z-transform climate layers to make axes comparable
for (i in 1:nlayers(climatelayers_ss)) {
  climatelayers_ss[[i]] <- (climatelayers_ss[[i]] - cellStats(climatelayers_ss[[i]], 'mean')) / cellStats(climatelayers_ss[[i]], 'sd') 
}
#seems like this worked okay, like climatelayers_ss[[1]] looks like it's still got plenty of information in it.

#this doesn't work, and not sure what's going wrong, but theoretical terra implementation. If something goes wrong below and numbers don't turn out as expected..come back to this area to check
#for (i in 1:nlyr(climatelayers_ss)) {
#  climatelayers_ss[[i]] <- (climatelayers_ss[[i]] - terra::global(climatelayers_ss[[i]], 'mean', na.rm = TRUE)) / terra::global(climatelayers_ss[[i]], 'sd', na.rm = TRUE) 
#}

climatelayers_ss_cropped = crop(climatelayers_ss, extent(-150,-50,15,60)) #this extent is most of North America (minus some parts of Canada), most of South America (except the tip) and some of Africa and Europe. 

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

## Calculate climate hypervolume based on each species range map

#this is done so nicely, but is causing problems! let's skip dplyr and do a simple loop.
spp_list$climate_vol <- 999

#for(i in 10:nrow(spp_list)){
#cardinal testing
for(i in 4){
  
  species <- spp_list$file[i]
  print(species)
  
  climate_vol <- climate_hypervolume(species)
  
  spp_list$climate_vol[i] <- climate_vol
  
  #and actually! let's save ALONG THE WAY
  write.csv(spp_list,
            paste0(output_dir, "climate_niche_breadth_mbbs_bioclim2.1_30s_upto_", spp_list$spp_name[i], ".csv"),
            row.names = F)
  
}





spp_hypervol <- spp_list %>%
  filter(!is.na(file)) %>%
  filter(file != "Alauda_arvensis_22717415.shp", file != "Phylloscopus_borealis_22715316.shp") %>%
  mutate(climate_vol = map_dbl(file, ~climate_hypervolume(.))) %>%
  dplyr::select(file, spp_name, file_binomial, aou, english_common_name, family, genus, species, binomial, old_genus, 
                new_binomial, matched_filename, species_code, climate_vol)
  
write.csv(spp_hypervol,
            paste0(output_dir, "climate_niche_breadth_mbbs_ACFL.csv"),
            row.names = F)


#okay well, it's stalling out. We should cut the species list down to one species n see if it's actually broke or just taking a long time to run bc ""lots"" (not THAT many lots) of species


