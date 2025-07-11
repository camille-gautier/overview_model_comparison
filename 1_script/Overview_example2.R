#! ---------------------------------------------------------------------------------------
#!
#! Description       : Create graphs and maps to compare various model for different variables
#!
#! Authors           : Camille Gautier <camille.gautier@ucalgary.ca>
#!
#! Creation date     : 2025-04-09 09:54:28
#!
#! Modifications     : 2025-05-15 16:18:13 Script reviewing + FUSE model integration (Cyril Thébault)
#!
#! Comments          :
#!
#! ---------------------------------------------------------------------------------------

library(ncdf4)
library(ncdf4.helpers)
library(ggplot2)
library(readxl)
library(dplyr)
library(CFtime)
library(lubridate)
library(sf)
library(viridisLite)
library(tidyr)
library(dplyr)
library(ggrepel)
library(ggpubr)
library(reshape2)
library(gridExtra)
library(RColorBrewer)

options(scipen = 100, digits = 4)

####---------------------------------------------
####
#### CONTROLE FILE
####
####---------------------------------------------

setwd("/Users/camillegautier/Documents/Visualization/github")

#Catchment name/number for the title of the sheet
catchment_name="Hay_river"

#Paths
#SUMMA outputs path, daily data with the variables:
#scalarSWE, scalarTotalET, mLayerVolFracWat and scalarRainfall in netcdf
summa_path="0_data/Example_2/SUMMA/output/run1_day.nc"

#Mizuroute output path, daily data with the variable:
#IRFroutedRunoff in netcdf
mizuroute_path="0_data/Example_2/SUMMA/output/"

#summa_zLocalAttributes.nc file
summa_attributes="0_data/Example_2/SUMMA/input/summa_zLocalAttributes.nc"

#MESH outputs path, daily data with the variables:
#RFF_D.csv, Q0_D.csv, ET_D.csv, SWE_D.csv, THLQSOL_D_IG1.csv, THLQSOL_D_IG3.csv and THLQSOL_D_IG3.csv
mesh_path="0_data/Example_2/MESH/output/"

#MESH forcing file path
mesh_forcing="0_data/Example_2/MESH/input/MESH_forcing.nc"

#MESH model path that includes MESH_drainage_database.nc
mesh_attributes="0_data/Example_2/MESH/input/MESH_drainage_database.nc"

#FUSE
fuse_path=""

#Observed streamflow 
Observed="0_data/Example_2/Observed/07OB001_DISCHARGE.csv"
#obs="/Users/camillegautier/Documents/Hay_river/new_discretization/DISCHARGE_DATA.csv"

#spatial discretisation to consider: "lumped", "subbasin"
#     - SUMMA: "lumped" or "subbasin"
#     - MESH: "lumped" or "subbasin"
#     - FUSE: "lumped"

spa = "subbasin"

#Target ID: focus on a specific subbasin or hru (can be found in the basin shapefile or attributes), if NA the outlet is chosen by default
#It must be consistent for all the models
#Only used if spa != "lumped"
target_ID=82029954

#TRUE if you want to aggregate all the subbasins upstream the target subbasin or FALSE if you only want the specific target subbasin
#Only used if spa == "subbasin"
upstream= TRUE


#Shapefiles
basin_shp="0_data/Example_2/Shapefiles/hay_basin_ag.shp"
river_shp="0_data/Example_2/Shapefiles/hay_river_ag.shp"

#Output_path
output_path="2_output/Example_2/"

####---------------------------------------------
####
#### FUNCTIONS
####
####---------------------------------------------

# find basin outlet
get_outlet <- function(basin_info){
  
  # Extract unique subbasins with no NextDownID
  outlet <- sapply(basin_info[names(basin_info) != "FUSE"], function(df) {
    unique(df$subbasin[is.na(df$NextDownID)|df$NextDownID==0|df$NextDownID==-9999])
  })
  
  if (length(unique(outlet)) == 1) {
    outlet <- unique(outlet)
    return(outlet)
  } else {
    stop("Mismatch: Models do not share the same subbasin with NextDownID = 0 or NA.")
  }
  
}

# compile information about ID (HRU, subbasin, lumped), area and NextDownID (only at the subbasin level)
get_info_bassin <- function(attribute_list, shp_list){
  
  basin_info=list()
  
  if (attribute_list$SUMMA!=""){
    
    # attributes
    summa_attributes =  attribute_list$SUMMA
    
    #open attributes file
    nc_attributes=nc_open(summa_attributes)
    hru <- ncvar_get(nc_attributes, "hruId")
    subbasin <- ncvar_get(nc_attributes, "hru2gruId")
    area <- ncvar_get(nc_attributes, "HRUarea")
    nc_close(nc_attributes)
    
    # shp
    summa_shp =  shp_list$SUMMA
    
    #open shp
    basinshp <- read_sf(summa_shp$basin)
    rivershp <- read_sf(summa_shp$river)
    
    if ("NextDownID" %in% names(basinshp)==TRUE){
      NextDownID <- basinshp$NextDownID[match(subbasin, basinshp$COMID)]
    }
    if ("NextDownID" %in% names(rivershp)==TRUE){
      NextDownID <- rivershp$NextDownID[match(subbasin, rivershp$COMID)]
    }
    
    info_df = data.frame(hru = hru, subbasin = subbasin, lumped = NA, area = area, NextDownID = NextDownID)
    
    basin_info$SUMMA <- info_df
    
  }
  
  if (attribute_list$MESH!=""){
    
    # attributes
    mesh_attributes =  attribute_list$MESH
    
    #open attributes file
    nc_attributes=nc_open(mesh_attributes)
    subbasin <- ncvar_get(nc_attributes, "subbasin")
    area <- ncvar_get(nc_attributes, "GridArea")
    NextDownID_Index <- ncvar_get(nc_attributes, "Next")
    NextDownID_Index[NextDownID_Index == 0] <- NA
    NextDownID <- subbasin[NextDownID_Index]
    nc_close(nc_attributes)
    
    info_df = data.frame(hru = subbasin, subbasin = subbasin, lumped = NA, area = area, NextDownID = NextDownID)
    
    basin_info$MESH <- info_df 
    
  }
  
  if (attribute_list$FUSE!=""){
    
    basinshp <- read_sf(summa_shp$basin)
    
    area = sum(basinshp$HRU_area)
    
    info_df = data.frame(hru = NA, subbasin = NA, lumped = NA, area = area, NextDownID = NA)
    
    basin_info$FUSE <- info_df 
    
  }
  
  outlet_id = get_outlet(basin_info)
  
  basin_info <- lapply(basin_info, function(df) {
    df$lumped <- outlet_id
    return(df)
  })
  
  return(basin_info)
}

# read models outputs and format them in a same way
read_results_models <- function(model_path_list, basin_info) {
  
  output_dfs <- lapply(names(model_path_list), function(x) list())
  names(output_dfs) <- names(model_path_list)
  
  #Read mizuroute
  #Streamflow routed
  if (model_path_list$SUMMA[["mizuRoute"]]!=""){
    
    mizu_file = list.files(model_path_list$SUMMA[["mizuRoute"]], full.names = TRUE, pattern = "mizuRoute")
    
    #open nc file
    nc_sim=nc_open(mizu_file)
    
    # Get raw time values
    time_raw <- ncvar_get(nc_sim, "time")
    
    # Get origin from units
    time_units <- ncatt_get(nc_sim, "time", "units")$value
    
    # Extract origin date string
    origin_str <- sub("seconds since ", "", time_units)
    origin <- as.POSIXct(origin_str, tz = "UTC")
    
    # Create actual POSIXct time values
    time_vals <- origin + time_raw
    
    
    dates <- as.Date(time_vals, format="%Y-%m-%d")
    
    #get subbasin
    subbasin <- ncvar_get(nc_sim, "basinID")
    
    # get streamflow
    qsim <- ncvar_get(nc_sim, "IRFroutedRunoff")
    
    #close nc file
    nc_close(nc_sim)
    
    Q=data.frame(qsim)
    
    if (ncol(Q)>1) {
      Q=t(Q)
    } else {
      Q=qsim
    }
    
    # create streamflow dataframe
    Q=data.frame(Q)
    colnames(Q) = subbasin
    rownames(Q) = dates
    
    # Find the minimum date (start of simulation)
    start_date <- min(as.Date(rownames(Q)))
    
    # Keep only data after the first year
    Q <- Q[as.Date(rownames(Q)) >= (start_date + 365), ,drop = FALSE]
    
    output_dfs$SUMMA$Qrouted <- Q
  }
  
  if (model_path_list$SUMMA[["SUMMA"]]!=""){
    #open outputs
    #Read SUMMA
    nc_sim=nc_open(model_path_list$SUMMA[["SUMMA"]])
    time_raw <- ncvar_get(nc_sim, "time")
    # Get raw time values
    # Get origin from units
    time_units <- ncatt_get(nc_sim, "time", "units")$value
    # Extract origin date string
    origin_str <- sub("seconds since ", "", time_units)
    origin <- as.POSIXct(origin_str, tz = "UTC")
    dates <- as.Date(time_vals, format="%Y-%m-%d")
    # Create actual POSIXct time values
    time_vals <- origin + time_raw
    #get subbasin (gru) & hru
    subbasin <- ncvar_get(nc_sim, "gruId")
    hru <- ncvar_get(nc_sim, "hruId")
    # Find the minimum date (start of simulation)
    start_date <- min(as.Date(dates))
    #Q
    if (any(startsWith(names(nc_sim$var), "scalarTotalRunoff"))){
      var_name <- grep("^scalarTotalRunoff", names(nc_sim$var), value = TRUE)[1]
      Q=data.frame(ncvar_get(nc_sim, var_name))
      
      if (ncol(Q)>1) {
        Q=t(Q)
      }
      colnames(Q) = hru
      Q <- as.data.frame(Q*86400*1000)
      rownames(Q) = dates
      # Keep only data after the first year
      Q <- Q[as.Date(rownames(Q)) >= (start_date + 365),  ,drop = FALSE]
      output_dfs$SUMMA$Q<- Q
    }
    #SWE
    if (any(startsWith(names(nc_sim$var), "scalarSWE"))){
      var_name <- grep("^scalarSWE", names(nc_sim$var), value = TRUE)[1]
      SWE=data.frame(ncvar_get(nc_sim, var_name))
      if (ncol(SWE)>1) {
        SWE=t(SWE)
      }
      
      SWE <- as.data.frame(SWE)
      colnames(SWE) = hru
      rownames(SWE) = dates
      # Keep only data after the first year
      SWE <- SWE[as.Date(rownames(SWE)) >= (start_date + 365),  ,drop = FALSE]
      output_dfs$SUMMA$SWE<- SWE
      #ET
    }
    if (any(startsWith(names(nc_sim$var), "scalarTotalET"))){
      var_name <- grep("^scalarTotalET", names(nc_sim$var), value = TRUE)[1]
      ET=data.frame(ncvar_get(nc_sim, var_name))
      if (ncol(ET)>1) {
        ET=t(ET)
      }
      ET <- as.data.frame(-1*3600*24*ET)
      colnames(ET) = hru
      rownames(ET) = dates
      # Keep only data after the first year
      ET <- ET[as.Date(rownames(ET)) >= (start_date + 365), ,drop = FALSE]
      output_dfs$SUMMA$ET <- ET
      #Soil Moisture
    }
    if (any(startsWith(names(nc_sim$var), "mLayerVolFracWat"))){
      var_name <- grep("^mLayerVolFracWat", names(nc_sim$var), value = TRUE)[1]
      SOIL = ncvar_get(nc_sim, var_name)
      dim=nc_sim$dim$gru$vals
      # Create a list for valid soil layers
      soil_list <- list()
      if (dim(dim)==1){
        num_layers <- dim(SOIL)[1]
        SOIL=t(SOIL)
        # Loop through each soil layer
        for (i in 1:num_layers) {
          layer_data <- SOIL[,i]
          if (!any(as.numeric(layer_data) == -9999)) {
            # Transpose and convert to data frame
            SOIL_I=as.data.frame(data.frame(layer_data))
            colnames(SOIL_I) <- hru
            # Keep only data after the first year
            rownames(SOIL_I) <- dates
            SOIL_I <- SOIL_I[as.Date(rownames(SOIL_I)) >= (start_date + 365),,drop = FALSE ]
            output_dfs$SUMMA[[paste0("SOIL_", i)]] <- SOIL_I
          }
        }
      } else {
        num_layers <- dim(SOIL)[2]
        # Loop through each soil layer
        for (i in 1:num_layers) {
          layer_data <- SOIL[, i, ]
          if (!any(as.numeric(layer_data) == -9999)) {
            # Transpose and convert to data frame
            SOIL_I=as.data.frame(t(data.frame(layer_data)))
            SOIL_I <- as.data.frame(SOIL_I)
            colnames(SOIL_I) <- hru
            rownames(SOIL_I) <- dates
            # Keep only data after the first year
            output_dfs$SUMMA[[paste0("SOIL_", i)]] <- SOIL_I
            SOIL_I <- SOIL_I[as.Date(rownames(SOIL_I)) >= (start_date + 365),,drop = FALSE ]
          }
        }
      }
    }
    if (any(startsWith(names(nc_sim$var), "pptrate"))){
      var_name <- grep("^pptrate", names(nc_sim$var), value = TRUE)[1]
      P=data.frame(ncvar_get(nc_sim, var_name))
      if (ncol(P)>1) {
        P=t(P)
      }
      P <- as.data.frame(3600*24*P)
      colnames(P) = hru
      rownames(P) = dates
      # Keep only data after the first year
      P <- P[as.Date(rownames(P)) >= (start_date + 365), ,drop = FALSE]
      output_dfs$SUMMA$P <- P
    }
    nc_close(nc_sim)
    
  }
  #Read MESH
  if (model_path_list$MESH!=""){
    
    # List the files
    mesh_files= list.files(model_path_list$MESH,pattern=".csv",full.names=TRUE)
    
    # Go through the different variables
    for (file in mesh_files) {
      
      # Define the name of each variable 
      name <- sub("\\.csv$", "", basename(file))
      if (grepl(pattern = "QO", x = name)){
        name_var="Qrouted"
      } else if (grepl(pattern = "RFF", x = name)){
        name_var="Q"
      } else if (grepl(pattern = "ET", x = name)){
        name_var="ET"
      } else if (grepl(pattern = "SNO", x = name)){
        name_var="SWE"
      } else if (grepl(pattern = "THLQSOL", x = name)){
        layerNumber = as.numeric(sub(".*IG", "", name))
        name_var=paste0("SOIL_", layerNumber)
      } else if (grepl(pattern = "PREC", x = name)){
        name_var="P"
      }
      
      # read the csv
      data=read.csv(file,sep=",",header = FALSE)
      
      # remove empty columns
      data <- data[, colSums(!is.na(data)) > 0]
      
      # format the dataframe
      dates = as.Date(data[,"V1"])
      data[, "V1"] <- NULL
      
      rownames(data) <- dates
      colnames(data) <- basin_info$MESH$subbasin
      
      # Information of instant Q (generated by each subbasin without routing) not available here
      # if (name_var=="Qrouted"){
      #   
      #   output_dfs$MESH$Q = sweep(data, 2, basin_info$MESH$area, FUN = function(q, a) (q * 86400 / a) * 1000)
      #   
      # }
      
      
      output_dfs$MESH[[name_var]]<- data
      # assign(name_var, data)
      
    }
    #Reading mesh forcing file to extract P
    if (file.exists(file =mesh_forcing)==TRUE){
      nc_forcing=nc_open(mesh_forcing)
      P_net=ncvar_get(nc_forcing,"RDRS_v2.1_A_PR0_SFC")
      subbasin=data.frame(c(1:nrow(ncvar_get(nc_forcing,"subbasin"))),ncvar_get(nc_forcing,"subbasin"))
      colnames(subbasin) <- c("rank","COMID")
      P_out=as.data.frame(P_net)*3600
      # Get raw time values
      time_raw <- ncvar_get(nc_forcing, "time")
      
      # Get origin from units
      time_units <- ncatt_get(nc_forcing, "time", "units")$value
      
      # Extract origin date string
      origin_str <- sub("hours since ", "", time_units)
      origin <- as.POSIXct(origin_str, tz = "UTC")
      
      # Create actual POSIXct time values
      time_vals <- origin + time_raw*3600
      
      colnames(P_out) = basin_info$MESH$subbasin
      P_out$day= as.Date(time_vals, format="%Y-%m-%d")
      # Keep only dates present in other MESH variables (e.g., Qrouted)
      
      P <- P_out[P_out$day %in% dates, , drop = FALSE]
      
      P_ag = aggregate(. ~ day, data = P, FUN = sum)
      rownames(P_ag) = P_ag$day
      P_ag$day = NULL
      
      output_dfs$MESH$P<- P_ag
    }
    
  }  
  
  #Read FUSE
  if (model_path_list$FUSE!=""){
    
    # outputs
    #fuse_file = list.files(model_path_list$FUSE, full.names = TRUE)
    
    #open outputs
    nc_sim=nc_open(model_path_list$FUSE)
    
    # Get raw time values
    time_raw <- ncvar_get(nc_sim, "time")
    
    # Get origin from units
    time_units <- ncatt_get(nc_sim, "time", "units")$value
    
    # Extract origin date string
    origin_str <- sub("days since ", "", time_units)
    origin <- as.Date(origin_str, tz = "UTC")
    
    # Create actual POSIXct time values
    time_vals <- origin + time_raw
    
    
    dates <- as.Date(time_vals, format="%Y-%m-%d")
    # Find the minimum date (start of simulation)
    start_date <- min(dates)
    
    #get subbasin (gru) & hru
    subbasin <- NA
    hru <- target_ID
    
    #Q
    if ("qsurf" %in% names(nc_sim$var)==TRUE){
      Q <- as.data.frame(ncvar_get(nc_sim, "qsurf"))
      colnames(Q) = hru
      rownames(Q) = dates
      # Keep only data after the first year
      Q <- Q[as.Date(rownames(Q)) >= (start_date + 365), , drop = FALSE]
      
      output_dfs$FUSE$Q<- Q
    }
    
    #Qrouted
    if ("q_routed" %in% names(nc_sim$var)==TRUE){
      Qrouted <- as.data.frame(data.frame(ncvar_get(nc_sim, "q_routed")*1000))
      colnames(Qrouted) = hru
      rownames(Qrouted) = dates
      # Keep only data after the first year
      Qrouted <- Qrouted[as.Date(rownames(Qrouted)) >= (start_date + 365),  , drop = FALSE]
      output_dfs$FUSE$Qrouted<- Qrouted
    }
    
    #SWE
    if ("swe_tot" %in% names(nc_sim$var)==TRUE){
      SWE <- data.frame(ncvar_get(nc_sim, "swe_tot"))
      colnames(SWE) = hru
      rownames(SWE) = dates
      # Keep only data after the first year
      SWE <- SWE[as.Date(rownames(SWE)) >= (start_date + 365), , drop = FALSE ]
      
      output_dfs$FUSE$SWE<- SWE
    }
    #ET
    # Get all variable names starting with "evap"
    evap_vars <- names(nc_sim$var)[startsWith(names(nc_sim$var), "evap")]
    
    if(length(evap_vars) > 0) {
      
      # Initialize ET sum as NULL
      ET_sum <- NULL
      
      for(varname in evap_vars) {
        # Extract data for this evap variable
        evap_data <- ncvar_get(nc_sim, varname)
        
        # If this is the first evap variable, start ET_sum with its data
        if(is.null(ET_sum)) {
          ET_sum <- evap_data
        } else {
          # Sum element-wise
          ET_sum <- ET_sum + evap_data
        }
      }
      
      # Convert summed vector/matrix to data.frame (adjust if 1 HRU or multiple HRUs)
      ET_df <- as.data.frame(ET_sum)
      
      colnames(ET_df) <- hru
      rownames(ET_df) <- dates
      
      # Keep only data after the first year, preserving dataframe structure
      ET_df <- ET_df[as.Date(rownames(ET_df)) >= (start_date + 365), , drop = FALSE]
      
      # Store output
      output_dfs$FUSE$ET <- ET_df
    }
    
    if ("ppt" %in% names(nc_sim$var)==TRUE){
      
      P <- as.data.frame(data.frame(ncvar_get(nc_sim, "ppt")))
      colnames(P) = hru
      rownames(P) = dates
      # Keep only data after the first year
      P <-P[as.Date(rownames(P)) >= (start_date + 365),  , drop = FALSE]
      
      output_dfs$FUSE$P <- P
      
    }
    nc_close(nc_sim)
    
  }
  
  
  return(output_dfs)
}

# find which IDs are upstream the target subbasin
FindUpstream <- function(IDs, NextDownID, focusID, upstream_comids = c(focusID)) {
  
  upstream_rows <- IDs[!(is.na(NextDownID)) & NextDownID == focusID]
  
  new_upstream <- upstream_rows[!upstream_rows %in% upstream_comids]
  
  upstream_comids <- c(upstream_comids, new_upstream)
  
  for (up_comid in new_upstream) {
    upstream_comids <- FindUpstream(IDs, NextDownID, up_comid, upstream_comids)
  }
  
  return(as.character(upstream_comids))
}

#Preparing data depending on:
#     - target ID (focus on a specific subbasin or hru, if NA use the outlet) 
#     - upstream condition (if TRUE: take every subbasin upstream, if FALSE: use only the target subbasin)
#     - spatial discretisation (lumped, subbasin or hru)
# Maybe should be split in 2 function (one for spatial consistency, the other for aggregation if upstream)
prepare_data <- function(output_dfs, basin_info, target_ID, upstream , spa ){
  
  if(spa != "subbasin" & upstream){
    message("Upstream option available only at subbasin scale, it has been set to FALSE")
    upstream = FALSE
  }
  
  # Use outlet if NA
  if (is.na(target_ID)==TRUE){
    target_ID = get_outlet(basin_info)
  }
  
  target_ID = as.character(target_ID)
  
  # Prepare data at the same spatial scale
  output_spa <- list()
  for(model in names(basin_info)){
    
    output_spa[[model]] <- list()
    
    # Loop over all variables except Qrouted
    variables_to_process <- setdiff(names(output_dfs[[model]]), "Qrouted")
    
    for (var_name in variables_to_process) {
      
      var_df <- output_dfs[[model]][[var_name]]
      
      # check the spatial discretisation of the result matrix
      group_levels <- c("hru", "subbasin", "lumped")
      group_counts <- sapply(group_levels, function(g) length(unique(basin_info[[model]][[g]])))
      if (!(ncol(var_df) %in% group_counts)) {
        stop("Cannot match the number of columns to any known spatial level")
      }
      matrix_spa <- group_levels[which(group_counts == ncol(var_df))]
      matrix_spa <- matrix_spa[length(matrix_spa)]
      
      if(which(group_levels == spa) < which(group_levels == matrix_spa)){
        output_spa[[model]][[var_name]] <- NA
        next
      }
      
      if (ncol(var_df)>1){
        # get ids and area
        area  <- setNames(basin_info[[model]]$area, basin_info[[model]][[spa]])
        
        #prepare the correspondance between the origin spatial discretisation (matrix_spa) and the target (spa)
        grouped_cols <- split(seq_along(names(area)), names(area))
        
        # Aggregate the outputs at the target spatial discretisation (spa) level
        weighted_var_df <- sweep(var_df, 2, area, FUN = "*")
        agg_df <- as.data.frame(sapply(grouped_cols, function(cols) {
          total_area <- sum(area[cols])
          rowSums(weighted_var_df[, cols, drop = FALSE]) / total_area
        }))
        agg_df <- agg_df[, order(as.numeric(colnames(agg_df))), drop = FALSE]
        output_spa[[model]][[var_name]] <- agg_df
      }
      else {
        output_spa[[model]][[var_name]] <- var_df
      }
      
    }
  }
  
  # Aggregate upstream subbasins
  if (upstream==TRUE){
    
    # get upstream ID from target_ID
    catchment_ids <- FindUpstream(IDs = basin_info$MESH$subbasin, NextDownID = basin_info$MESH$NextDownID, focusID = target_ID) ## /!\ To change: Should not depend on a model
    
    # get area for relevant catchment_ids
    relevant_hru_info <- basin_info[[1]][basin_info[[1]][[spa]] %in% catchment_ids, ]  ## /!\ To change: Should use a model that has information at spa level
    
    area_by_group <- relevant_hru_info %>%
      group_by(.data[[spa]]) %>%
      summarise(total_area = sum(area), .groups = "drop")
    
    # prepare weights for weighted average
    weights <- area_by_group$total_area / sum(area_by_group$total_area)
    names(weights) <- as.character(area_by_group[[spa]])
    
    # Loop over each model (e.g., SUMMA, MESH)
    output_format <- lapply(output_spa, function(model) {
      
      # For each model, aggregate all variables (Q, ET, SWE, etc.)
      lapply(model, function(df) {
        
        # Subset to relevant catchments
        df_subset <- df[, catchment_ids, drop = FALSE]
        
        # Align weights with column names
        matching_weights <- weights[colnames(df_subset)]
        
        # Apply weights using sweep (multiply columns by weights)
        weighted_df <- sweep(df_subset, 2, matching_weights, `*`)
        
        # Sum across columns to get weighted mean per time step
        vec <- rowSums(weighted_df)
        as.data.frame(setNames(list(vec), target_ID), row.names = names(vec), check.names = FALSE)
        
      })
    })
    
  } else {
    # return target_id
    output_format <- output_spa
  }
  
  # Include Qrouted
  # Maybe include it also in output_spa
  for (model_name in names(output_format)) {
    if ("Qrouted" %in% names(output_dfs[[model_name]])) {
      if (model_name=="FUSE"){
        output_format[[model_name]]$Qrouted <- output_dfs[[model_name]]$Qrouted*basin_info[[model]]$area/(1000000*24*3600)
      }
      else {
        output_format[[model_name]]$Qrouted <- output_dfs[[model_name]]$Qrouted[, target_ID, drop = FALSE]
      }
    }
  }
  
  
  return(list(output_spa = output_spa, output_format = output_format))
}

# Compute the water budget
compute_water_budget <- function(output_format) {
  vars <- c("Q", "ET", "P")
  
  water_budget <- lapply(names(output_format), function(model_name) {
    model_data <- output_format[[model_name]]
    
    # Skip model if not all required variables exist
    if (!all(vars %in% names(model_data))) return(NULL)
    
    target_ID <- colnames(model_data$Q)[1]  
    
    # Extract and add date column
    df_list <- lapply(vars, function(v) {
      df <- model_data[[v]]
      df$date <- as.Date(rownames(model_data[[v]]))
      colnames(df)[1] <- v
      return(df)
    })
    
    # Merge Q, ET, P by date
    merged <- Reduce(function(x, y) full_join(x, y, by = "date"), df_list)
    
    # Add year
    merged$year <- year(merged$date)
    
    # Count non-NA values per year for each variable
    yearly_counts <- merged %>%
      group_by(year) %>%
      summarise(
        P_n = sum(!is.na(P)),
        ET_n = sum(!is.na(ET)),
        Q_n = sum(!is.na(Q)),
        .groups = "drop"
      )
    
    # Filter years where all 3 variables have >= 365 values
    valid_years <- yearly_counts %>%
      filter(P_n >= 365, ET_n >= 365, Q_n >= 365) %>%
      pull(year)
    
    # Keep only data from valid years
    merged <- merged %>% filter(year %in% valid_years)
    
    # Annual sums
    annual <- merged %>%
      group_by(year) %>%
      summarise(
        P = sum(P, na.rm = TRUE),
        ET = sum(ET, na.rm = TRUE),
        Q = sum(Q, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        P_ET = P / ET,
        Q_P = Q / P
      )
    
    # Format each variable as 1-column data.frame with years as rownames
    budget_list <- lapply(c("P", "ET", "Q", "P_ET", "Q_P"), function(var) {
      df <- data.frame(annual[[var]])
      rownames(df) <- annual$year
      colnames(df) <- target_ID
      return(df)
    })
    names(budget_list) <- c("P", "ET", "Q", "P_ET", "Q_P")
    
    return(budget_list)
  })
  
  # Set names and remove NULLs
  names(water_budget) <- names(output_format)
  water_budget <- Filter(Negate(is.null), water_budget)
  
  return(water_budget)
}

####---------------------------------------------
####
#### MAPS
####
####---------------------------------------------

map_hru <- function(basin_shp, river_shp, target_ID, upstream){
  river=read_sf(river_shp)
  basin=read_sf(basin_shp)
  
  if (is.na(target_ID)==TRUE){
    target_ID = get_outlet(basin_info)
  }
  
  if (upstream==TRUE){
    catchment_ids <- FindUpstream(basin$COMID, basin$NextDownID, target_ID)
  } else {
    catchment_ids <- target_ID
  }
  
  # Avoid to have a fully red map when all subbasins are selected
  if(length(catchment_ids) == nrow(basin)){
    catchment_ids = NULL
  }
  
  if (spa=="lumped"){
    ggplot() + 
      geom_sf(data = basin, fill = "white", color = "gray72",size=0.1) +  # Optional: coloring the basin
      geom_sf(data = river, color="#1C86EE", size = 0.2) + 
      theme_bw() +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            legend.position = "right",
            legend.box = "vertical",
            panel.grid = element_blank())+
      ggtitle("Catchment map")+
      theme(plot.title = element_text(size = 10, hjust = 0.5)) 
  } else {
    ggplot() + 
      geom_sf(data = basin, fill = "white", color = "gray72",size=0.1) +  # Optional: coloring the basin
      geom_sf(data = basin[basin$COMID %in% catchment_ids,],fill= "firebrick1", color="firebrick4", size = 0.5) + 
      geom_sf(data = river, color="#1C86EE", size = 0.2) + 
      theme_bw() +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            legend.position = "right",
            legend.box = "vertical",
            panel.grid = element_blank())+
      ggtitle("Catchment map")+
      theme(plot.title = element_text(size = 10, hjust = 0.5)) 
  }
  
}


####---------------------------------------------
####
#### GRAPHS
####
####---------------------------------------------

monthly_regime_boxplot <- function(mylist, var_name) {
  color <- switch(var_name,
                  Q = "blue",
                  ET = "darkgreen",
                  SWE = "black")
  title <- switch(var_name,
                  Q = "Streamflow regime",
                  ET = "Evapotranspiration regime",
                  SWE = "SWE regime")
  y_lab <- switch(var_name,
                  Q = "Streamflow [mm/month]",
                  ET = "ET [mm/month]",
                  SWE = "SWE [mm/month]")
  
  combined_data <- bind_rows(
    lapply(names(mylist), function(model_name) {
      
      if (!var_name %in% names(mylist[[model_name]])) return(NULL)  # skip if variable missing
      
      df <- mylist[[model_name]][[var_name]]
      data.frame(
        date = as.Date(rownames(df)),
        value = df[,1],     # Maybe a better way to do that... here column 1 hard coded.
        model = model_name
      )
    })
  )
  
  var_data <- combined_data %>%
    mutate(YearMonth = format(date, "%Y-%m")) %>%
    group_by(YearMonth, model) %>%  # Group by both YearMonth and model
    summarize(discharge = sum(value, na.rm = TRUE), .groups = "drop") %>%  # Ensure model stays
    mutate(Date = as.Date(paste0(YearMonth, "-01"), format = "%Y-%m-%d"),  # Create Date for plotting
           Month = factor(format(Date, "%b"), levels = month.abb))  # Extract Month as short names
  
  # Colour line dictionary
  model_names <- names(output_format) 
  n_models <- length(model_names)
  colors <- suppressWarnings(brewer.pal(min(n_models, 8), "Set2"))
  if (n_models > 8) {
    colors <- colorRampPalette(brewer.pal(8, "Set2"))(n_models)
  }
  model_colors <- setNames(colors, model_names)
  
  # Short month names
  gg <- ggplot(var_data, aes(x = Month, y = discharge, color=model, fill=model)) +
    geom_boxplot( alpha = 0.3, outliers = FALSE) +
    scale_y_continuous(limits = c(0, NA)) +
    scale_color_manual(values = model_colors) +
    scale_fill_manual(values = model_colors)+
    labs(
      title =title,
      x = "Month",
      y = y_lab
    ) +
    theme_bw() +
    theme(
      axis.title.x = element_text(color = color),  # Color of x-axis title
      axis.title.y = element_text(color = color),  # Color of y-axis title
      axis.text.x = element_text(color = color),    # Color of x-axis text
      axis.text.y = element_text(color = color, angle = 90, hjust = 0.5),    # Color of y-axis text
      plot.title = element_text(color = color),     # Color of plot title
      panel.border = element_rect(color = color),
      legend.position="none")
  
  gg
  
}

annual_barplot <- function(mylist, var_name) {
  # Determine variable-specific labels and color
  color <- switch(var_name,
                  Q = "blue",
                  ET = "darkgreen",
                  SWE = "black")
  title <- switch(var_name,
                  Q = "Annual streamflow",
                  ET = "Annual evapotranspiration",
                  SWE = "Annual SWE")
  y_lab <- switch(var_name,
                  Q = "Streamflow [mm/year]",
                  ET = "ET [mm/year]",
                  SWE = "SWE [mm/year]")
  
  # Combine all models into one data frame
  combined_data <- bind_rows(
    lapply(names(mylist), function(model_name) {
      
      if (!var_name %in% names(mylist[[model_name]])) return(NULL)  # skip if variable missing
      
      df <- mylist[[model_name]][[var_name]]
      data.frame(
        date = as.Date(rownames(df)),
        value = df[,1],     # Maybe a better way to do that... here column 1 hard coded.
        model = model_name
      )
    })
  )
  
  
  # Calculate hydrological year and annual sum
  combined_annual <- combined_data %>%
    mutate(
      Year = as.numeric(format(date, "%Y")),
      Month = as.numeric(format(date, "%m")),
      Hydro_Year = ifelse(Month >= 10, Year + 1, Year)
    ) %>%
    group_by(model, Hydro_Year)%>%
    summarise(Annual_Value = sum(value, na.rm = TRUE),
              Na_Flag = any(is.na(value)), .groups = "drop")
  
  # Colour line dictionary
  model_names <- names(output_format) 
  n_models <- length(model_names)
  colors <- suppressWarnings(brewer.pal(min(n_models, 8), "Set2"))
  if (n_models > 8) {
    colors <- colorRampPalette(brewer.pal(8, "Set2"))(n_models)
  }
  model_colors <- setNames(colors, model_names)
  
  gg <- ggplot(combined_annual, aes(x = factor(Hydro_Year), y = Annual_Value, fill = model)) +
    geom_bar(stat = "identity", position = position_dodge(preserve = "single")) +
    scale_fill_manual(values = model_colors)+
    labs(
      title = title,
      y = y_lab,
      x = "Hydrological year [Oct–Sep]"
    ) +
    scale_x_discrete(breaks = seq(min(combined_annual$Hydro_Year),max(combined_annual$Hydro_Year),by=5)) +  # Set x-axis limits
    theme_bw() +
    theme(
      axis.title.x = element_text(color = color),  # Color of x-axis title
      axis.title.y = element_text(color = color),  # Color of y-axis title
      axis.text.x = element_text(color = color),    # Color of x-axis text
      axis.text.y = element_text(color = color, angle = 90, hjust = 0.5),    # Color of y-axis text
      plot.title = element_text(color = color),     # Color of plot title  # Color of legend title
      panel.border = element_rect(color = color),
      legend.position="none" # Reduce legend key size
    )
  
  gg
}

cumulative_frequency_curve <- function(mylist, var_name){
  color <- switch(var_name,
                  Q = "blue",
                  ET = "darkgreen",
                  SWE = "black")
  title <- switch(var_name,
                  Q = "CDF streamflow",
                  ET = "CDF evapotranspiration",
                  SWE = "CDF SWE")
  y_lab <- switch(var_name,
                  Q = "Streamflow [mm/day] (log)",
                  ET = "ET [mm/day]",
                  SWE = "SWE [mm/day]")
  
  combined_data <- bind_rows(
    lapply(names(mylist), function(model_name) {
      
      if (!var_name %in% names(mylist[[model_name]])) return(NULL)  # skip if variable missing
      
      df <- mylist[[model_name]][[var_name]]
      data.frame(
        date = as.Date(rownames(df)),
        value = df[,1],     # Maybe a better way to do that... here column 1 hard coded.
        model = model_name
      )
    })
  )
  
  quantiles <- bind_rows(
    lapply(names(mylist), function(model_name) {
      probs = c(0, 0.25, 0.5, 0.75, 1)
      vec <- quantile(mylist[[model_name]][[var_name]][,1], probs = probs, na.rm = TRUE)
      df <- data.frame(
        quantile = probs,
        value = vec,
        model = model_name
      )
      return(df)
    })
  )
  
  
  combined_data <- combined_data %>%
    group_by(model) %>%
    arrange(value) %>%
    mutate(cum_freq = seq_along(value) / length(value))
  
  # Colour line dictionary
  model_names <- names(output_format) 
  n_models <- length(model_names)
  colors <- suppressWarnings(brewer.pal(min(n_models, 8), "Set2"))
  if (n_models > 8) {
    colors <- colorRampPalette(brewer.pal(8, "Set2"))(n_models)
  }
  model_colors <- setNames(colors, model_names)
  
  gg <- ggplot(combined_data, aes(x = cum_freq, y = value, color = model)) +
    geom_line() +
    scale_color_manual(values = model_colors) +
    geom_point(data = quantiles, aes(x = quantile, y = value, color = model), size = 3) +
    geom_text_repel(
      data = quantiles,
      aes(x = quantile, y = value, label = round(value, 2), color = model),
      size = 3
    )+
    labs(
      title = title,
      x = "Cumulative Frequency",
      y = y_lab,
      color = NULL
    ) +
    theme_bw() +
    theme(
      axis.title.x = element_text(color = color),
      axis.title.y = element_text(color = color),
      axis.text.x = element_text(color = color),
      axis.text.y = element_text(color = color, angle = 90, hjust = 0.5),
      plot.title = element_text(color = color),
      legend.position="none",
      panel.border = element_rect(color = color))
  
  
  if(var_name == "Q" | var_name == "Qrouted"){
    gg <- gg + scale_y_log10()
  }
  
  gg
}

cumulative_curve <- function(mylist, var_name){
  color <- switch(var_name,
                  Q = "blue",
                  ET = "darkgreen",
                  SWE = "black")
  title <- switch(var_name,
                  Q = "Cumulative streamflow",
                  ET = "Cumulative evapotranspiration",
                  SWE = "Cumulative SWE")
  y_lab <- switch(var_name,
                  Q = "Streamflow [mm]",
                  ET = "ET [mm]",
                  SWE = "SWE [mm]")
  
  
  combined_data <- bind_rows(
    lapply(names(mylist), function(model_name) {
      
      if (!var_name %in% names(mylist[[model_name]])) return(NULL)  # skip if variable missing
      
      df <- mylist[[model_name]][[var_name]]
      data.frame(
        date = as.Date(rownames(df)),
        value = df[,1],     # Maybe a better way to do that... here column 1 hard coded.
        model = model_name
      )
    })
  )
  # Calculate maximal annual values
  # Prepare data with cumulative sum
  combined_data_cum <- combined_data %>%
    arrange(date) %>%
    group_by(model) %>%
    mutate(Cumulative = cumsum(replace_na(value, 0)))
  
  # Colour line dictionary
  model_names <- names(output_format) 
  n_models <- length(model_names)
  colors <- suppressWarnings(brewer.pal(min(n_models, 8), "Set2"))
  if (n_models > 8) {
    colors <- colorRampPalette(brewer.pal(8, "Set2"))(n_models)
  }
  model_colors <- setNames(colors, model_names)
  
  
  # Calculate quantiles and min/max values
  gg <- ggplot(combined_data_cum, aes(x = date, y = Cumulative,color=model)) +
    geom_line() + 
    scale_color_manual(values = model_colors) +
    labs(
      title = title,
      x = "Date",
      y = y_lab
    ) +
    theme_bw()+
    theme(axis.title.x = element_text(color = color),  # Color of x-axis title
          axis.title.y = element_text(color = color),  # Color of y-axis title
          axis.text.x = element_text(color = color),    # Color of x-axis text
          axis.text.y = element_text(color = color, angle = 90, hjust = 0.5),    # Color of y-axis text
          plot.title = element_text(color = color),     # Color of plot title
          legend.position="none",   # Color of legend title
          panel.border = element_rect(color = color))
  
  gg
}

soil_fraction <- function(mylist){
  
  combined_data <- bind_rows(
    lapply(names(mylist), function(model_name) {
      model_data <- mylist[[model_name]]
      
      # Extract all "SOIL_*" variables
      soil_vars <- names(model_data)[startsWith(names(model_data), "SOIL_")]
      
      if (length(soil_vars) == 0) return(NULL)  # Skip if no SOIL_* variables
      
      # Combine into a data.frame (one column per SOIL layer)
      soil_data <- do.call(cbind, lapply(soil_vars, function(var) {
        df <- model_data[[var]]
        # Ensure it is a 1-column dataframe or matrix
        setNames(df, var)
      }))
      
      # Add date and model
      data.frame(
        date = as.Date(rownames(soil_data)),
        soil_data,
        model = model_name,
        check.names = FALSE
      )
    })
  )
  
  
  
  # Reshape to long format
  long_data <- combined_data %>%
    pivot_longer(cols = starts_with("SOIL_"), 
                 names_to = "layer", 
                 values_to = "value")
  
  # Colour line dictionary
  model_names <- names(output_format) 
  n_models <- length(model_names)
  colors <- suppressWarnings(brewer.pal(min(n_models, 8), "Set2"))
  if (n_models > 8) {
    colors <- colorRampPalette(brewer.pal(8, "Set2"))(n_models)
  }
  model_colors <- setNames(colors, model_names)
  
  ggplot(long_data, aes(x = date, y = value, color = model, linetype = layer)) +
    geom_line() +
    scale_color_manual(values = model_colors) +
    scale_y_continuous(limits = c(0, NA)) +
    labs(
      title = "Soil Water Fraction",
      x = "Date",
      y = "Water Fraction",
      linetype = "Soil Layer"
    ) + guides(color = "none")+
    theme(legend.position = "inside")+
    theme_bw()
  
}

#Water budget

water_budget_simple <- function(mylist, Observed) {
  
  water_budget <- compute_water_budget(mylist)
  
  
  combined_data <-   bind_rows(
    lapply(names(water_budget), function(model_name) {
      model_vars <- water_budget[[model_name]]
      
      years <- as.integer(rownames(model_vars[[1]]))
      
      # Combine all variables into a single data.frame
      df <- data.frame(
        year = years,
        lapply(model_vars, function(df_var) df_var[[1]])
      )
      
      df$model <- model_name
      return(df)
    }),
    .id = NULL
  )
  
  combined_data$diff=combined_data$P-combined_data$Q-combined_data$ET
  
  combined_data$year = as.character(combined_data$year)
  
  mean_diff <- combined_data %>%
    group_by(model) %>%
    summarise(mean_diff = mean(diff, na.rm = TRUE), .groups = "drop")
  mean_diff$mean_diff <- as.numeric(mean_diff$mean_diff)
  
  # Colour line dictionary
  model_names <- names(output_format) 
  n_models <- length(model_names)
  colors <- suppressWarnings(brewer.pal(min(n_models, 8), "Set2"))
  if (n_models > 8) {
    colors <- colorRampPalette(brewer.pal(8, "Set2"))(n_models)
  }
  model_colors <- setNames(colors, model_names)
  
  ggplot(data = combined_data) +
    geom_point(aes(x = year, y = diff, color = model)) +
    scale_color_manual(values = model_colors) +
    labs(title = "Water budget",
         x = "Year", y = "P - Q - ET [mm/year]") +
    scale_x_discrete(breaks = seq(min(as.numeric(as.character(combined_data$year))), max(as.numeric(as.character(combined_data$year))), by = 5))+
    geom_text_repel(data = mean_diff,
                    aes(x = min(combined_data$year),
                        y = mean_diff,
                        label = paste("mean =", round(mean_diff, 2), "mm/year"),
                        color = model),
                    size = 3, hjust = 0, direction = "y", nudge_x = 0) +
    geom_abline(slope = 0, linetype = "dashed", color = "red") +
    coord_cartesian(ylim = c(min(combined_data$diff, 0),max(0, combined_data$diff)))+
    theme_bw()+
    theme(legend.position="none")
  
}

qqplot <- function(mylist,Observed) {
  combined_data <- bind_rows(
    lapply(names(mylist), function(model_name) {
      
      if (!"Qrouted" %in% names(mylist[[model_name]])) return(NULL)  # skip if variable missing
      
      df <- mylist[[model_name]]$Qrouted
      data.frame(
        date = as.Date(rownames(df)),
        value = df[,1],     # Maybe a better way to do that... here column 1 hard coded.
        model = model_name
      )
    })
  )
  
  # Colour line dictionary
  model_names <- names(output_format) 
  n_models <- length(model_names)
  colors <- suppressWarnings(brewer.pal(min(n_models, 8), "Set2"))
  if (n_models > 8) {
    colors <- colorRampPalette(brewer.pal(8, "Set2"))(n_models)
  }
  model_colors <- setNames(colors, model_names)
  
  na_strings <- c("-9999","","NA","#VALUE!")
  
  obs=data.frame(read.csv(Observed,sep=",",na.strings = na_strings))
  
  # Find column name that matches "date" (case-insensitive)
  date_col <- names(obs)[tolower(names(obs)) == "date"]
  # Convert DATE column to Date class
  obs$date <- as.Date(obs[[date_col]]) 
  
  merge = merge(combined_data, obs, by="date")
  
  colnames(merge) <- c("date","value","model","Date","Observed")
  
  # Extract each model's data frame and name the column accordingly
  model_data <- lapply(names(mylist), function(model_name) {
    if (!"P" %in% names(mylist[[model_name]])) return(NULL)
    
    df <- mylist[[model_name]]$P
    data.frame(
      date = as.Date(rownames(df)),
      value = df[,1]
    ) %>% 
      rename(!!model_name := value)  # Rename 'value' to model name
  }) 
  # Merge all model columns by date
  P_wide <- Reduce(function(x, y) merge(x, y, by = "date", all = TRUE), model_data)
  all.equal(P_wide$SUMMA, P_wide$MESH)
  
  ggplot() +
    geom_point(data=merge,aes(x = Observed, y = value, color = model), show.legend = FALSE ) +
    scale_color_manual(values = model_colors) +
    labs(title = "Streamflow QQ Plot",
         x = "Observed Streamflow [m3/s]", y = "Simulated streamflow [m3/s]") +
    coord_cartesian(ylim = c(min(combined_data$value, 0),max(0, combined_data$value)))+
    scale_x_continuous()+
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    theme_bw() +
    theme(
      axis.title.x = element_text(color = "blue"),  # Color of x-axis title
      axis.title.y = element_text(color = "blue"),  # Color of y-axis title
      axis.text.x = element_text(color = "blue"),    # Color of x-axis text
      axis.text.y = element_text(color = "blue", angle = 90, hjust = 0.5),    # Color of y-axis text
      plot.title = element_text(color = "blue"),     # Color of plot title  # Color of legend title
      panel.border = element_rect(color = "blue")
    )
  
}

#Timeserie streamflow
timeserie_streamflow <- function(mylist, Observed){
  
  combined_data <- bind_rows(
    lapply(names(mylist), function(model_name) {
      
      if (!"Qrouted" %in% names(mylist[[model_name]])) return(NULL)  # skip if variable missing
      
      df <- mylist[[model_name]]$Qrouted
      data.frame(
        date = as.Date(rownames(df)),
        value = df[,1],     # Maybe a better way to do that... here column 1 hard coded.
        model = model_name
      )
    })
  )
  # Colour line dictionary
  model_names <- names(output_format) 
  n_models <- length(model_names)
  colors <- suppressWarnings(brewer.pal(min(n_models, 8), "Set2"))
  if (n_models > 8) {
    colors <- colorRampPalette(brewer.pal(8, "Set2"))(n_models)
  }
  model_colors <- setNames(colors, model_names)
  
  na_strings <- c("-9999","","NA","#VALUE!")
  
  obs=data.frame(read.csv(Observed,sep=",",na.strings = na_strings))
  
  # Find column name that matches "date" (case-insensitive)
  date_col <- names(obs)[tolower(names(obs)) == "date"]
  # Convert DATE column to Date class
  obs$date <- as.Date(obs[[date_col]]) 
  
  merge = merge(combined_data, obs, by="date")
  
  colnames(merge) <- c("date","value","model","Date","Observed")
  
  # Extract each model's data frame and name the column accordingly
  model_data <- lapply(names(mylist), function(model_name) {
    if (!"P" %in% names(mylist[[model_name]])) return(NULL)
    
    df <- mylist[[model_name]]$P
    data.frame(
      date = as.Date(rownames(df)),
      value = df[,1]
    ) %>% 
      rename(!!model_name := value)  # Rename 'value' to model name
  }) 
  # Merge all model columns by date
  P_wide <- Reduce(function(x, y) merge(x, y, by = "date", all = TRUE), model_data)
  all.equal(P_wide$SUMMA, P_wide$MESH)
  
  ggplot() +
    geom_line(data=merge,aes(x = date, y = value, color = model), show.legend = FALSE ) +
    geom_line(data=merge,aes(x = date, y = Observed, linetype = "Observed"), color = "black", size = 0.5) +
    #+geom_bar(stat=)
  # Manual linetype legend for "Observed"
  scale_linetype_manual(
    name = "",  # Legend title (empty)
    values = c("Observed" = "dashed")
  ) +
    # Optional: keep color mapping for visual appearance only
    scale_color_manual(values = model_colors) +
    labs(title = "Streamflow timeseries",
         x = "Date", y = "Streamflow [m3/s]") +
    coord_cartesian(ylim = c(min(combined_data$value, 0),max(0, combined_data$value)))+
    theme_bw() +
    theme(
      axis.title.x = element_text(color = "blue"),  # Color of x-axis title
      axis.title.y = element_text(color = "blue"),  # Color of y-axis title
      axis.text.x = element_text(color = "blue"),    # Color of x-axis text
      axis.text.y = element_text(color = "blue", angle = 90, hjust = 0.5),    # Color of y-axis text
      plot.title = element_text(color = "blue"),     # Color of plot title  # Color of legend title
      panel.border = element_rect(color = "blue"),
      legend.position = "inside",
      legend.position.inside = c(0.9,0.9)
    )+ theme(legend.background=element_blank())
}


#Metric functions
KGE <- function(sim, obs) {
  
  suis_NA <- is.na(obs) | is.na(sim)
  
  qsim = sim[!suis_NA]
  qobs = obs[!suis_NA]
  
  sigma_Obs <- sd(qobs)
  sigma_Sim <- sd(qsim)
  alpha <- sigma_Sim/sigma_Obs
  ere <- cor(qobs, 
             qsim, method = "pearson")
  beta <- sum(qsim)/sum(qobs)
  return(list(KGE = 1-sqrt((ere-1)^2+(beta-1)^2+(alpha-1)^2),
              alpha = alpha,
              beta = beta,
              ere = ere))
} 

NSE <- function(sim, obs) {
  
  suis_NA <- is.na(obs) | is.na(sim)
  
  qsim = sim[!suis_NA]
  qobs = obs[!suis_NA]
  
  numerator = sum((qobs - qsim)^2)
  denominator = sum((qobs - mean(qobs))^2)
  
  nse = 1-(numerator/denominator)
  return(nse)
} 

Pbias <- function(sim, obs) {
  
  suis_NA <- is.na(obs) | is.na(sim)
  
  qsim = sim[!suis_NA]
  qobs = obs[!suis_NA]
  
  numerator = abs(sum(qobs-qsim))*100
  denominator = sum(qobs)
  
  pbias = numerator/denominator
  return(pbias)
} 

#Calculate metrics
calculate_metrics <- function(mylist, Observed){
  
  metrics=list()
  na_strings <- c("-9999","","NA","#VALUE!")
  
  obs=data.frame(read.csv(Observed,sep=",",na.strings = na_strings))
  
  # Find column name that matches "date" (case-insensitive)
  date_col <- names(obs)[tolower(names(obs)) == "date"]
  # Convert DATE column to Date class
  rownames(obs)= as.Date(obs[[date_col]]) 
  
  for (model in names(mylist)) {
    
    merged <- merge(mylist[[model]]$Qrouted, obs, by = "row.names")
    colnames(merged) = c("date","sim","DATE","obs")
    metrics[[model]]$KGE=KGE(merged$sim,merged$obs)$KGE
    metrics[[model]]$NSE=NSE(merged$sim,merged$obs)
    metrics[[model]]$PBIAS=Pbias(merged$sim,merged$obs)
  }
  
  metrics_text <- paste(
    unlist(lapply(names(metrics), function(model) {
      metric_lines <- sapply(names(metrics[[model]]), function(metric_name) {
        value <- round(metrics[[model]][[metric_name]], 2)
        paste0("  ", metric_name, " = ", value, ifelse(metric_name == "PBIAS", "%", ""))
      })
      c(paste0(model, ":"), metric_lines, "")  # Add model name and a blank line
    })),
    collapse = "\n"
  )
  
  ggplot() +
    annotate("text", x = 0, y = 1, label = paste0("Metrics: \n",metrics_text), size = 4) +
    theme_void()
}

####---------------------------------------------
####
#### SUMMARY SHEETS
####
####---------------------------------------------

summary_sheets_graphs <- function(catchment_name, output_dfs, output_format, target_ID, upstream, basin_shp, river_shp, output_path, Observed){
  
  if (is.na(target_ID)==TRUE){
    target_ID = get_outlet(basin_info)
  }
  
  if ((target_ID==get_outlet(basin_info) & upstream) | spa=="lumped")  {
    # Create an empty plot layout
    layout <- matrix(c(
      c(1, 5, 9, 13, 13),
      c(2, 6, 10, 13, 13),
      c(3, 7, 11, 14, 14),
      c(4, 8, 12, 15, 16),
      c(17, 18, 18, 18, 19))
      , nrow = 5, byrow = TRUE)
    
    
    CDF_Q <- cumulative_frequency_curve (output_format,"Q") 
    CUM_Q <- cumulative_curve (output_format,"Q")
    BAR_Q <- annual_barplot (output_format,"Q")
    REG_Q <- monthly_regime_boxplot (output_format,"Q")
    CDF_ET <- cumulative_frequency_curve (output_format,"ET")
    CUM_ET <- cumulative_curve (output_format,"ET")
    BAR_ET <- annual_barplot (output_format,"ET")
    REG_ET <- monthly_regime_boxplot (output_format,"ET")
    CDF_SWE <- cumulative_frequency_curve(output_format,"SWE")
    CUM_SWE <- cumulative_curve (output_format,"SWE")
    BAR_SWE <- annual_barplot (output_format,"SWE")
    REG_SWE <- monthly_regime_boxplot (output_format,"SWE")
    WB_SIM <- water_budget_simple(output_format) + theme(legend.position = "none")
    SOIL <- soil_fraction(output_format)
    MAP <- map_hru(basin_shp, river_shp, target_ID, upstream)
    QQ <- qqplot(output_format, Observed)
    Q_timeseries <- timeserie_streamflow(output_format, Observed)
    metrics <- calculate_metrics(output_format, Observed)
    
    mylegend <- get_legend(
      BAR_ET + 
        labs(fill = paste(catchment_name,"\n\nHRU: ",target_ID, "\n\nUpstream: ", upstream, "\n\n","Models")) +             
        theme(legend.position = "right",legend.text=element_text(size=12),legend.title=element_text(size=16)))
    
    gg <- arrangeGrob(
      CDF_Q, CUM_Q, BAR_Q, REG_Q, 
      CDF_ET, CUM_ET, BAR_ET, REG_ET, 
      CDF_SWE, CUM_SWE, BAR_SWE, REG_SWE, 
      MAP, SOIL, WB_SIM, mylegend,
      QQ, Q_timeseries, metrics,
      layout_matrix = layout
    )
  }
  else {
    # Create an empty plot layout
    layout <- matrix(c(
      c(1, 5, 9, 13, 13),
      c(2, 6, 10, 13, 13),
      c(3, 7, 11, 14, 14),
      c(4, 8, 12, 15, 16)), 
      nrow = 4, byrow = TRUE)
    
    
    CDF_Q <- cumulative_frequency_curve (output_format,"Q") 
    CUM_Q <- cumulative_curve (output_format,"Q")
    BAR_Q <- annual_barplot (output_format,"Q")
    REG_Q <- monthly_regime_boxplot (output_format,"Q")
    CDF_ET <- cumulative_frequency_curve (output_format,"ET")
    CUM_ET <- cumulative_curve (output_format,"ET")
    BAR_ET <- annual_barplot (output_format,"ET")
    REG_ET <- monthly_regime_boxplot (output_format,"ET")
    CDF_SWE <- cumulative_frequency_curve(output_format,"SWE")
    CUM_SWE <- cumulative_curve (output_format,"SWE")
    BAR_SWE <- annual_barplot (output_format,"SWE")
    REG_SWE <- monthly_regime_boxplot (output_format,"SWE")
    WB_SIM <- water_budget_simple(output_format) + theme(legend.position = "none")
    SOIL <- soil_fraction(output_format)
    MAP <- map_hru(basin_shp, river_shp, target_ID, upstream)
    
    mylegend <- get_legend(
      BAR_ET + 
        labs(fill = paste(catchment_name,"\n\nHRU: ",target_ID, "\n\nUpstream: ", upstream, "\n\n","Models")) +             
        theme(legend.position = "right",legend.text=element_text(size=12),legend.title=element_text(size=16)))
    
    gg <- arrangeGrob(
      CDF_Q, CUM_Q, BAR_Q, REG_Q, 
      CDF_ET, CUM_ET, BAR_ET, REG_ET, 
      CDF_SWE, CUM_SWE, BAR_SWE, REG_SWE, 
      MAP, SOIL, WB_SIM, mylegend,
      layout_matrix = layout
    )
  }
  
  ggsave(filename = paste0(output_path,catchment_name,"_Graphs.png"),
         plot = gg,
         height = 16, width = 22)
}


####---------------------------------------------
####
#### MAIN
####
####---------------------------------------------
model_path_list <- list(SUMMA  = c(SUMMA = summa_path, mizuRoute = mizuroute_path), 
                        MESH = mesh_path,
                        FUSE = fuse_path)

attribute_list <- list(SUMMA  = summa_attributes,
                       MESH = mesh_attributes,
                       FUSE = fuse_path)

shp_list <- list(SUMMA  = list(basin = basin_shp,
                               river = river_shp),
                 MESH = list(basin = basin_shp,
                             river = river_shp))

basin_info <- get_info_bassin(attribute_list, shp_list)
output_dfs <- read_results_models(model_path_list, basin_info)
outputs <- prepare_data(output_dfs, basin_info, target_ID, upstream, spa)
output_format <- outputs$output_format


summary_sheets_graphs(catchment_name,output_dfs, output_format, target_ID, upstream, basin_shp, river_shp, output_path, Observed)
