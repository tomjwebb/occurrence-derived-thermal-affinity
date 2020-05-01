## R functions to implement the temperature matching of marine species occurrence records

get_temp_summ_by_sp <- function(sp_id, save_all_recs = TRUE, check_match = TRUE, use_defaults = TRUE){
  # function to get OBIS records for a given species and match to IAP and Bio-Oracle temperatures
  # if save_all_recs == TRUE, this full matched dataset will be saved before summarising
  # function returns a summary of temperature affinity of the species
  
  # First, check if this species has previously been matched, if check_match is TRUE
  # (This is an option in case the user wants to update all records)
  
  if(check_match == TRUE){
    out_path <- ifelse(use_defaults, paste0(file.path(getwd()), "/t_matched_obis_recs"),
                       paste0(bespoke_path, "/t_matched_obis_recs"))
    # set species_done flag to FALSE
    species_done <- FALSE
    
    # check if the directory exists:
    if(dir.exists(file.path(out_path))){
      # get list of aphia ids from the files in this directory:
      done_files <- list.files(path = out_path)
      done_aphias <- done_files %>%
        word(sep = "_") %>%
        str_replace(pattern = "aphia", replacement = "") %>%
        as.integer()
      
      # check if species sp_id is in this list
      species_done <- sp_id %in% done_aphias
    }
    
    # get the file if it exists
    if(species_done == TRUE){
      which_file <- which(done_aphias == sp_id)
      
      # what if more than one file matches?
      if(length(which_file) > 1){
        # find the most recent one
        which_file <- which_file[order(done_files[which_file], decreasing = TRUE) == 1]
      }
      
      sp_temp <- read_csv(file = paste(out_path, done_files[which_file], sep = "/"),
                          col_types = cols(
                            decimalLongitude = col_double(),
                            decimalLatitude = col_double(),
                            depth = col_double(),
                            eventDate = col_datetime(format = ""),
                            scientificName = col_character(),
                            aphiaID = col_integer(),
                            year = col_integer(),
                            month = col_character(),
                            depth0 = col_double(),
                            grid_depth = col_integer(),
                            grid_lon = col_integer(),
                            grid_lat = col_integer(),
                            iap_t = col_double(),
                            iap_sst = col_double(),
                            iap_sbt = col_double(),
                            iap_grid_bottom_depth = col_integer(),
                            bo_sst = col_double(),
                            bo_sbt = col_double()
                          )) %>%
        t_summary() %>%
        mutate(species_id = sp_id) %>%
        dplyr::select(species_id, everything())
    } else {
      sp_temp <- get_obis_recs(species_id = sp_id) %>%
        get_iap_gridded_t() %>%
        get_bio_oracle_t() %>%
        save_full_recs(save_recs = save_all_recs) %>%
        t_summary() %>%
        mutate(species_id = sp_id) %>%
        dplyr::select(species_id, everything())
    }
  } else {
    sp_temp <- get_obis_recs(species_id = sp_id) %>%
      get_iap_gridded_t() %>%
      get_bio_oracle_t() %>%
      save_full_recs(save_recs = save_all_recs) %>%
      t_summary() %>%
      mutate(species_id = sp_id) %>%
      dplyr::select(species_id, everything())
  }
  
  sp_temp
}

## Getting OBIS records
get_obis_recs <- function(species_id, missing_check = FALSE,
                          fields = c("decimalLongitude", "decimalLatitude",
                                     "minimumDepthInMeters", "maximumDepthInMeters", "depth",
                                     "eventDate", "year", "month",
                                     "scientificName", "aphiaID")
){
  # Fuction to get OBIS records for a given species_id, which must be a recognised WoRMS Aphia ID
  
  # NB OBIS returns records from all taxa gathered under the same valid Aphia ID; the aphia ID returned is that of the taxon as recorded, not necessarily the valid ID, so in order that the final dataset is correctly named we add back in the 'correct' ID here as valid_AphiaID
  
  if(missing_check == TRUE){
    # catch invalid / unrecognised AphiaIDs here - but recommend doing this prior to calling these functions
    if(length(checklist(taxonid = species_id)) > 1){
      # get OBIS records for a given species ID, add year and month, set negative and missing depth to 0
      obis_recs <- robis::occurrence(
        taxonid = species_id, fields = fields) %>%
        as_tibble()
      if(!"year" %in% names(obis_recs)){
        obis_recs <- obis_recs %>% 
          mutate(year = "NA")}
      if(!"month" %in% names(obis_recs)){
        obis_recs <- obis_recs %>% 
          mutate(month = NA)}
      obis_recs <- mutate(obis_recs,
                          depth = as.numeric(depth),
                          year = formatC(year),
                          month = formatC(as.numeric(month), width = 2, flag = "0"),
                          depth0 = case_when(
                            is.na(depth) ~ 0,
                            depth < 0 ~ 0,
                            TRUE ~ depth),
                          valid_AphiaID = species_id)
    } else {
      # at present just returns an empty tibble, which causes problems with other functions further down the pipeline, hence recommend checking AphiaIDs prior to calling
      obis_recs <- tibble()
    }
  } else {
    obis_recs <- robis::occurrence(
      taxonid = species_id, fields = fields) %>%
      as_tibble()
    if(!"year" %in% names(obis_recs)){
      obis_recs <- obis_recs %>% 
        mutate(year = "NA")}
    if(!"month" %in% names(obis_recs)){
      obis_recs <- obis_recs %>% 
        mutate(month = NA)}
    obis_recs <- mutate(
      obis_recs,
      depth = as.numeric(depth),
      year = formatC(year),
      month = formatC(as.numeric(month), width = 2, flag = "0"),
      depth0 = case_when(is.na(depth) ~ 0,
                         depth < 0 ~ 0,
                         TRUE ~ depth),
      valid_AphiaID = species_id
    )
  }
  # return the OBIS records
  obis_recs
}

## Matching OBIS Records to Bio-ORACLE Temperatures
get_bio_oracle_t <- function(obis_recs, layercodes = c("BO_sstmean", "BO2_tempmean_bdmean"),
                             use_defaults = TRUE){
  # Function to match a set of OBIS occurrence recrods to mean SST and SBT from Bio-ORACLE
  # Set path for where these two temperature datasets will be stored
  bo_path <- ifelse(use_defaults,
                    paste0(file.path(getwd()), "/biooracle"),
                    paste0(bespoke_path, "/biooracle"))
  if(!dir.exists(file.path(bo_path))){
    dir.create(path = file.path(bo_path),
               recursive = FALSE, showWarnings = FALSE)}
  # load the layers
  bo_t_dat <- load_layers(layercodes,
                          equalarea = TRUE, datadir = "biooracle")
  # Turn the OBIS occurrence locations into spatial points
  points <- SpatialPoints(
    obis_recs[,c("decimalLongitude", "decimalLatitude")],
    lonlatproj)
  # Reproject (could avoid this by setting equalarea = FALSE)
  points <- spTransform(points, equalareaproj)
  # Extract the temperatures for each point
  bo_sst <- raster::extract(bo_t_dat[[1]], points)
  bo_sbt <- raster::extract(bo_t_dat[[2]], points)
  bo_temp <- as_tibble(cbind(bo_sst, bo_sbt))
  
  # add these temperatures back to the OBIS records and return
  bind_cols(obis_recs, bo_temp)
  
}

## Matching OBIS Records to IAP Gridded Temperature
get_iap_gridded_t <- function(obis_recs, use_defaults = TRUE, bespoke_path = NULL){
  # Function to match a set of OBIS occurrence records to IAP gridded global temperature by month, year, lat, lon, and depth. Returns SST, SBT, and temperature at depth.
  # NB: This will take a long time to run first time, as it has to download the climate data. Subsequent runs will be much quicker, assuming the path remains the same
  # Process:
  # 1. get grid cell for occurrence dat
  # 2. group by month, year, deal with occs missing date
  # 3. get relevant IAP data file given month, year
  # 4. match occurrences to temperature by date, lat, lon, depth
  # 5. run over a all months / years
  
  #################################################
  # Function to run the temperature matching by latitude, longitude, time, and depth over all dates
  get_iap_temp_by_lat_lon_time_depth <- function(
    occurrence_dat = obis_recs, fpath = iap_path){
    
    # group the data by year and month, and get temperature for each record
    matched_temps <- occurrence_dat %>%
      group_by(year, month) %>%
      do(iap_t = get_t_by_grid(., fpath = fpath))
    matched_temps <- unnest(dplyr::select(matched_temps, iap_t))
    
    matched_temps
  }
  
  #################################################
  
  get_gridcell <- function(dll_df,
                           depth_vals = c(1, 5, seq(10, 100, by = 10), seq(120, 200, by = 20),
                                          seq(250, 900, by = 50), seq(1000, 1800, by = 100), 2000),
                           lon_vals = 1:360,lat_vals = seq(-89.5, 89.5)){
    
    # Function to get grid cell index from continuous lat/lon/depth given specified lat/lon/depth grid
    # Defaults here are for a global 1 degree grid, longitude in degrees East, latitude in degrees North, with standard depth bands following the World Ocean Atlas
    # Returns a grid index (depth, lon, lat) for one or a vector of locations
    # Accepts a 3-column tibble with cols depth, lon, and lat (in that order, col names not important)
    
    # find index for depth, lon and lat
    grid_depth <- findInterval(pull(dll_df[, 1]), depth_vals, all.inside = TRUE)
    grid_lon <- findInterval(pull(dll_df[, 2]), lon_vals, all.inside = TRUE)
    grid_lat <- findInterval(pull(dll_df[, 3]), lat_vals, all.inside = TRUE)
    # add to data frame
    grid_ids <- bind_cols(dll_df,
                          grid_depth = grid_depth, grid_lon = grid_lon, grid_lat = grid_lat)
    # return
    grid_ids
  }
  
  #################################################
  
  get_iap_by_month_year <- function(path = iap_path, month = "01", year = "1990") {
    # Function to get the IAP gridded temperature data for a given month and year and save it to a specified location
    # build up the file path to search
    fpath <- file.path(path, paste0("CZ16_1_2000m_Temp_year_", year, "_month_", month, ".nc"))
    # check for existence of required data; get it if it does not exist
    if (!file.exists(fpath)) {
      dir.create(dirname(fpath), recursive = TRUE, showWarnings = FALSE)
      dat_url <- paste0(
        "ftp://ds1.iap.ac.cn/ftp/cheng/CZ16_v3_IAP_Temperature_gridded_1month_netcdf/Monthly/",
        basename(fpath))
      download.file(dat_url, destfile = fpath)
    }
    
    # Open file as nc
    iap_t <- nc_open(filename = fpath)
    # get the sea temperature data for this month
    sea_temp <- ncvar_get(iap_t, varid = "temp")
    # close the file connection
    nc_close(iap_t)
    # return the temperature data
    sea_temp
  }
  
  #################################################
  
  get_t_by_grid <- function(loc_dat, fpath = iap_path){
    
    # Function to get temperature by grid cell for a given month and year
    
    # get month and year to query - assumes data has month, year variables and these are the same for all occurrences
    # (NB datasets with multiple months / years are dealt with in subsequent function)
    loc_m <- loc_dat$month[1]
    loc_y <- loc_dat$year[1]
    st_mat <- get_iap_by_month_year(month = loc_m, year = loc_y, path = fpath)
    st_loc <- st_mat[
      as.matrix(dplyr::select(loc_dat, grid_depth, grid_lon, grid_lat))]
    # sea surface temperature (temperature at depth layer 1)
    sst_loc <- st_mat[
      as.matrix(cbind(1, dplyr::select(loc_dat, grid_lon, grid_lat)))]
    # sea bottom temperature (deepest non-NA temp - will be sea bed in shallow water, t at 2000m in deep)
    # function to find the deepest non-NA temp
    min_na <- function(x){
      if(sum(is.na(x)) == 0){
        min_na <- length(x) + 1
      } else {
        min_na <- min(which(is.na(x)))
      }
      min_na[min_na > 1] <- min_na[min_na > 1] - 1
      min_na
    }
    # apply this over the full st_mat matrix
    sb_depth <- apply(st_mat, c(2,3), min_na)
    # get sea bottom depth for each observation
    sb_obs <- sb_depth[
      as.matrix(cbind(dplyr::select(loc_dat, grid_lon, grid_lat)))]
    # finally, get the temperature for this combination of depth, lon and lat
    sbt_loc <- st_mat[
      as.matrix(cbind(sb_obs, dplyr::select(loc_dat, grid_lon, grid_lat)))]
    # add the three temperature values back to the locations data, also add sb_depth to cross ref with survey depth
    loc_dat <- bind_cols(loc_dat,
                         iap_t = st_loc, iap_sst = sst_loc, iap_sbt = sbt_loc,
                         iap_grid_bottom_depth = sb_obs)
    loc_dat
  }
  
  # for later binding - get a unique row id for each record
  obis_recs <- mutate(rowid_to_column(obis_recs, "recordID"),
                      lon_deg_east = case_when(
                        decimalLongitude >= 0 ~ decimalLongitude,
                        decimalLongitude < 0 ~ 360 + decimalLongitude))
  # get grid cell for each depth0, lon and lat
  iap_occs <- obis_recs %>%
    dplyr::select(depth0, lon_deg_east, decimalLatitude) %>%
    get_gridcell()
  
  iap_occs <- bind_cols(dplyr::select(obis_recs, recordID, year, month), iap_occs)
  
  # filter out records with no date, or date out of range, or year but no month
  iap_occs <- iap_occs %>%
    mutate(yr = as.numeric(year)) %>%
    filter(!is.na(yr) & yr >= 1940 & yr < 2018 & month != "NA" & !is.na(month))
  
  if(nrow(iap_occs) > 0){
    # set path
    iap_path <- ifelse(use_defaults, paste0(file.path(getwd()), "/iap_gridded"),
                       paste0(bespoke_path, "/iap_gridded"))
    iap_temps <- get_iap_temp_by_lat_lon_time_depth(
      occurrence_dat = iap_occs, fpath = iap_path)
    iap_temps <- left_join(obis_recs,
                           dplyr::select(iap_temps, recordID, grid_depth:grid_lat,
                                         iap_t:iap_grid_bottom_depth), by = c("recordID"))
  } else {
    iap_temps <- mutate(obis_recs,
                        grid_depth = NA, grid_lon = NA, grid_lat = NA,
                        iap_t = NA, iap_sst = NA, iap_sbt = NA, iap_grid_bottom_depth = NA)
  }
  
  # Tidy up data
  iap_temps <- dplyr::select(iap_temps, -recordID, -lon_deg_east)
  
  # Return
  iap_temps
}

## Saving Temperature-Matched OBIS Records to File
save_full_recs <- function(rec_df, save_recs = TRUE, use_defaults = TRUE, bespoke_path = NULL){
  # if save_recs == TRUE, save the full set of obis records + IAP + BO temperature for a species

  if(save_recs == TRUE){
    out_path <- ifelse(use_defaults, paste0(file.path(getwd()),
      "/t_matched_obis_recs"),
    paste0(bespoke_path, "/t_matched_obis_recs"))
  
  if(!dir.exists(file.path(out_path))){
    dir.create(path = file.path(out_path), recursive = FALSE, showWarnings = FALSE)}
  
  # paste together the filename
  sp_filename <- paste0("aphia", rec_df$valid_AphiaID[1],
    "_obis_iap_bo_", Sys.Date(), ".csv")
  
  # write the file
  write_csv(x = rec_df, path = file.path(paste(out_path, sp_filename, sep = "/")))
  }

  # Return the (unchanged) data to pass to next function
  rec_df
}

## Summarising the Temperature Affinity of a Species
t_summary <- function(t_matched_dat){
  # Function to get a range of temperature summary stats from a matched obis-IAP-bio-oracle data frame
  counts <- summarise(t_matched_dat, n_obis_rec = n())
  missings <- miss_var_summary(dplyr::select(
    t_matched_dat, iap_t:iap_sbt, bo_sst:bo_sbt))
  missings_df <- t(missings[, "n_miss"])
  colnames(missings_df) <- paste0(pull(missings, variable), "_NA")
  missings_df <- as_tibble(missings_df)
  # define separate functions for 5% and 95% quantiles
  q5 <- function(x, na.rm = TRUE){stats::quantile(x, 0.05, na.rm = TRUE)}
  q95 <- function(x, na.rm = TRUE){stats::quantile(x, 0.95, na.rm = TRUE)}

  # get a range of summary stats over all variables in the dataset
  t_stats <- summarise_at(t_matched_dat,
    vars(iap_t:iap_sbt, bo_sst:bo_sbt),
    tibble::lst(mean, min, max, median, sd, mad, q5, q95), na.rm = TRUE)
  
  # Tidy up and return the species-level summary
  t_summ <- bind_cols(counts, missings_df, t_stats)

  t_summ
}