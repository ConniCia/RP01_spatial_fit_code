# Aim:
#   clean data so that it's ready to be read in during simulations
# Steps:
#   - pass country, type of spatial discretisation, scales
#     N.B. second element is a list of the same length as the first argument
#   - looped over passed arguments

wrapper_cleanRawData <- function (countries, grid_scales,
                                  scales_DEFAULT,
                                  grid_changes = "",
                                  grid_changes_Namibia = "") {
  
  for (i in seq_along(countries)) {
    
    country <- countries[i]
    
    # default
    cleanRawData(country, grid_scales, scales_DEFAULT, "DEFAULT")
    
    # changes to be applied to all countries
    for (grid_change in grid_changes) {
      cleanRawData(country, grid_scales, scales_DEFAULT, grid_change)
    }
    
    # changes to be applied only to Namibia
    if (country == "Namibia") {
      for (grid_change in grid_changes_Namibia) {
        cleanRawData(country, grid_scales, scales_DEFAULT, grid_change)
      }
    }
    
  }
  
}


cleanRawData <- function (country, grid_scales, scales_DEFAULT, grid_change) {
  
  
  tic(paste(country, grid_change, "preparation"))
  
  # grid_scale limit for sensitivity analyses ---------------------------------#
  if (grid_change != "DEFAULT")
    grid_scales <- intersect(grid_scales, scales_DEFAULT)
  
  # output folder -------------------------------------------------------------#
  dir_out <- file.path(dir_cln, country)
  if (grid_change != "")
    dir_out <- paste_(dir_out, grid_change)
  if(!dir.exists(dir_out))
    dir.create(dir_out)
  
  
  # what admins to merge (will only be applied to Namibia) --------------------#
  
  ## default to merge: admins that fall through at least one of the grid scales
  merge.me   <- c("330207", "330409", "330510", "331001", "331006", "331303")
  merge.with <- c("330206", "330402", "330509", "330702", "331005", "331302")
  delete.me <- NULL
  keep.me <- NULL
  
  if (grepl("hole", grid_change)) {
    
    ## to merge: admins that are surrounded by other admins like the hole of a bagel
    merge.me   <- c(merge.me,   c("330104", "330501", "330502", "330503", "330504", "330505", "330506", "330507", "330508", "331110"))
    merge.with <- c(merge.with, c("330103", "330509", "330509", "330509", "330509", "330509", "330509", "330509", "330509", "331102"))
    
  }
  
  if (grepl("flow", grid_change)) {
    
    ## to merge: codes of admins with no travel data -> merge to closest neighbour (best road connection, looking at road maps)
    merge.me   <- c(merge.me,   c("330304", "331002", "331004"))
    merge.with <- c(merge.with, c("330306", "331008", "331007"))
    
  }
  
  if (grepl("noNE", grid_change)) {
    
    ## to delete: codes of protruding admins in the north-east
    delete.me <- c("330405", "331301", "331302", "331303", "331304", "331305", "331306")
    
  }
  
  if (grepl("onlydenseN", grid_change)) {
    
    ## to keep: codes of high population density admins
    keep.me <- c("3304", "3307", "3309", "3310", "3311")
    keep.me <- unlist(lapply(keep.me, function(x) sprintf(paste0(x, "%02d"), 1:12)))
    
  }
  
  if (grepl("random", grid_change)) {
    
    ## This is done later in the code!
    
  }
  
  
  #----------------------------------------------------------------------------#
  # read input files and clean, list admin codes that have travel data         #
  #----------------------------------------------------------------------------#
  
  #-- names file --------------------------------------------------------------#
  
  # read in units file - (Kenya = admin level 1, Namibia = admin level 2)
  # only keep columns relating to admin code and name
  dropped <- if ( country == "Kenya" ) 2 else c(2,3)
  INPUT_names <- fread(paste0(dir_raw, country, "_names.txt"),
                       drop = dropped, # drop unused columns
                       col.names  = c("code", "name"),
                       colClasses = "character")
  # remove special characters from names
  INPUT_names[, name := gsub("[_/]", " ", name)]
  
  
  #-- population file ---------------------------------------------------------#
  
  # read in population per ~km^2 (=: pixel)
  # only keep columns relating to admin code, longitude, latitude, ambient population
  dropped <- if ( country == "Kenya" ) 4 else NULL
  col.classes <- if ( country == "Kenya" ) c("numeric", "numeric", "integer", "character", "character") else c("numeric", "numeric", "integer", "character")
  INPUT_pop <- fread(paste0(dir_raw, country, "_pop.txt"),
                     header = FALSE,
                     drop = dropped, # drop unused columns
                     col.names = c("lon", "lat", "pop", "code"),
                     colClasses = col.classes)
  
  
  #-- flow data file ----------------------------------------------------------#
  
  # read in flow data between admins
  if ( country == "Kenya" ) {
    INPUT_flow <- fread(paste0(dir_raw, country, "_flow_data.csv"),
                        select = 1:3, # only use specified columns
                        col.names  = c("origin", "destination", "amount"),
                        colClasses = c("character", "character", "integer", rep("character", 7)))
    # remove special characters from names
    INPUT_flow[, `:=`(origin      = gsub("[.]", " ", origin),
                      destination = gsub("[.]", " ", destination))]
    # replace names by codes
    for (col in c("origin", "destination")) {
      INPUT_flow <- merge(INPUT_flow, INPUT_names, by.x=col, by.y="name")
      setnames(INPUT_flow, "code", paste0("code_", col))
    }
    INPUT_flow[, c("origin", "destination") := NULL]
    # transform to origin-destination matrix
    INPUT_flow <- dcast.data.table(INPUT_flow, code_origin ~ code_destination, value.var = "amount", fill=0, drop = FALSE)
    INPUT_flow <- as.matrix(INPUT_flow, rownames="code_origin")
    
    # quality check
    if( !all(INPUT_flow >= 0) )
      stop("Something wrong with the flow data!!!")
  } else {
    INPUT_flow <- read.table(paste0(dir_raw, country, "_flow_data.csv"),
                             sep=",",
                             header = T,
                             stringsAsFactors=FALSE,
                             row.names=1) # first column of df actually represents rownames
    # remove special characters from names
    colnames(INPUT_flow) <- gsub("[.]", " ", colnames(INPUT_flow))
    # replace names by codes
    rownames(INPUT_flow) <- INPUT_names$code[match(rownames(INPUT_flow), INPUT_names$name)]
    colnames(INPUT_flow) <- INPUT_names$code[match(colnames(INPUT_flow), INPUT_names$name)]
    # transpose matrix (so that rows=origins, columns=destinations)
    INPUT_flow <- t(as.matrix(INPUT_flow))
    # quality check: sum of each row has to be between 0 and 1
    if( !all( INPUT_flow >= 0 & rowSums(INPUT_flow) < 1.0001 ) )
      stop("Something wrong with the flow data!!!")
  }
  
  # sort by code
  INPUT_flow <- INPUT_flow[order(rownames(INPUT_flow)), order(colnames(INPUT_flow))]
  # quality check
  if ( !identical(rownames(INPUT_flow), colnames(INPUT_flow)) )
    stop("Row and column names of the flow data do not correspond!!!")
  # remove rows and columns with 0 data
  remove <- vector("integer")
  for (i in nrow(INPUT_flow):1) {
    if (sum(INPUT_flow[i, ]) < 0.0001) {
      remove <- c(remove, i)
    }
  }
  if (length(remove) > 0)
    INPUT_flow <- INPUT_flow[-remove, -remove]
  
  
  # N.B. I don't need to work on INPUT_names anymore, it's only purpose was to sub admin names for admin codes
  rm(INPUT_names)
  
  
  #----------------------------------------------------------------------------#
  # merge admins                                                               #
  #----------------------------------------------------------------------------#
  
  if ( country == "Namibia" ) {
    
    # perform each merge separately
    for (i in 1:length(merge.me)) {
      # population
      if (merge.me[i] %in% INPUT_pop[, unique(code)]) {
        # change admin code
        INPUT_pop[code == merge.me[i], code := merge.with[i]]
      }
      # flow data
      if (merge.me[i] %in% rownames(INPUT_flow)) {
        # sum correct rows
        INPUT_flow[merge.with[i], ] <- INPUT_flow[merge.me[i], ] + INPUT_flow[merge.with[i], ]
        # sum correct columns
        INPUT_flow[, merge.with[i]] <- INPUT_flow[, merge.me[i]] + INPUT_flow[, merge.with[i]]
      }
    }
    # flow data: delete rows and columns
    INPUT_flow <- INPUT_flow[ !rownames(INPUT_flow) %in% merge.me, !colnames(INPUT_flow) %in% merge.me ]
    # sort by admin code
    INPUT_flow <- INPUT_flow[order(rownames(INPUT_flow)), order(colnames(INPUT_flow))]
    
  }
  
  
  #----------------------------------------------------------------------------#
  #-- print pixel population of all admins (not just the ones used in simulations) #
  #----------------------------------------------------------------------------#
  
  # with empty pixels, so I know how to colour them/where to draw the borders when plotting!
  if (grid_change == "DEFAULT")
    fwrite(INPUT_pop, file = file.path(dir_cln, paste0("my_total_", country, "_pixels.tsv")))
  
  
  #----------------------------------------------------------------------------#
  # remove admins that have no or 0 travel data or 0 ambient population        #
  #----------------------------------------------------------------------------#
  
  # delete rows referring to pixels with 0 ambient population
  INPUT_pop <- INPUT_pop[pop > 0]
  
  # list admin codes
  codes_INPUT_pop <- INPUT_pop[, unique(code)]
  
  # list admin codes
  codes_INPUT_flow <- rownames(INPUT_flow)
  
  # compute admin codes to KEEP: the ones that have both travel data AND population
  codes_INPUT <- intersect(codes_INPUT_pop, codes_INPUT_flow)
  # population: delete rows
  INPUT_pop <- INPUT_pop[code %in% codes_INPUT]
  # flow data: delete rows and columns
  INPUT_flow <- INPUT_flow[ rownames(INPUT_flow) %in% codes_INPUT, colnames(INPUT_flow) %in% codes_INPUT]
  
  
  #----------------------------------------------------------------------------#
  # filter admins according to the chosen sensitivity analysis                 #
  #----------------------------------------------------------------------------#
  
  # delete admins
  if (length(delete.me) > 0) {
    # population: delete rows
    INPUT_pop <- INPUT_pop[! code %in% delete.me]
    # flow data: delete rows and columns
    INPUT_flow <- INPUT_flow[ ! rownames(INPUT_flow) %in% delete.me, ! colnames(INPUT_flow) %in% delete.me]
  }
  
  # sample random admins
  if (grepl("random", grid_change)) {
    
    num <- str_extract_all(grid_change, "\\d+")[[1]]
    nr_samples <- num[1]
    set.seed(num[2])
    keep.me <- sample(unique(INPUT_pop$code), nr_samples, replace = FALSE)
    
  }
  
  # restrict to specific admins to keep
  if (length(keep.me) > 0) {
    # population: delete rows
    INPUT_pop <- INPUT_pop[code %in% keep.me]
    # flow data: delete rows and columns
    INPUT_flow <- INPUT_flow[rownames(INPUT_flow) %in% keep.me, colnames(INPUT_flow) %in% keep.me]
  }
  
  
  #----------------------------------------------------------------------------#
  #-- print pixel population of admins I use in C simulation ------------------#
  #----------------------------------------------------------------------------#
  
  # with empty pixels, so I know how to colour them/where to draw the borders when plotting!
  if (grid_change == "DEFAULT") {
    fwrite(INPUT_pop, file = file.path(dir_cln, paste0("my_C_", country, "_pixels.tsv")))
    fwrite(INPUT_pop[, .(unique_admins = sort(unique(code)))], 
           file = file.path(dir_cln, paste0("my_C_", country, "_admins.tsv")))
  }
  
  toc()
  
  
  #----------------------------------------------------------------------------#
  # Compute input files for C simulation (admin level)                         #
  #----------------------------------------------------------------------------#
  tic(paste(country, grid_change, "admin"))
  
  
  #-- population output -------------------------------------------------------#
  
  # compute population and population-weighted centroid by admin code
  dt_pop <- INPUT_pop[, .(apop = sum(pop),
                          clon = sum(pop * lon) / sum(pop),
                          clat = sum(pop * lat) / sum(pop)),
                      by = .(code)]
  # sort by admin code
  dt_pop <- setorder(dt_pop, code)
  # write "admin__patch_population"
  fwrite(dt_pop[, .(apop)],
         file=file.path(dir_out, "admin__patch_population.tsv"),
         col.names=FALSE)
  # write "admin__admin_population"
  fwrite(dt_pop[, .(apop)],
         file=file.path(dir_out, "admin__admin_population.tsv"),
         col.names=FALSE)
  # write "admin__patches_per_admin"
  fwrite(dt_pop[, .N, by=code][,.(N)],
         file=file.path(dir_out, "admin__patches_per_admin.tsv"),
         col.names=FALSE)
  # write "admin__metadata"
  metadata <- data.table(names =c("nr_admins",                   "nr_patches"),
                         values=c(dt_pop[,length(unique(code))], dt_pop[,.N]))
  fwrite(metadata,
         file=file.path(dir_out, "admin__metadata.tsv"),
         sep="\t",
         col.names=FALSE)
  
  
  #-- distance output ---------------------------------------------------------#
  
  # distance in m (shortest distance assuming earth is perfect sphere)
  distance <- distm(cbind(dt_pop[,clon], dt_pop[,clat]),
                    fun = distHaversine)
  # convert *upper diagonal* into integer vector
  distance <- as.integer(distance[t(upper.tri(distance))])
  # write "admin_patch_distance"
  fwrite(as.data.table(distance),
         file=file.path(dir_out, "admin__patch_distance.tsv"),
         col.names=FALSE)
  
  
  #-- flow data output --------------------------------------------------------#
  # make Namibian flow data integer by multiplying each row (=origin) by the corresponding population size
  mat_flow <- INPUT_flow
  if (country == "Namibia") {
    mat_flow <- INPUT_flow * dt_pop[, apop] + 0.5
    class(mat_flow) <- "integer"
  }
  # write "admin__admin_flow_data"
  fwrite(as.data.table(mat_flow),
         file=file.path(dir_out, "admin__admin_flow_data.tsv"),
         sep="\t",
         col.names=FALSE)
  
  toc()
  
  #----------------------------------------------------------------------------#
  # Compute input files for C simulation (loop over different grid scales)     #
  #----------------------------------------------------------------------------#
  
  for (grid_scale in grid_scales) {
    
    tic(paste(country, grid_change, grid_scale, "km"))
    #-- admin population output -----------------------------------------------#
    
    # compute population by admin code
    dt_admin_pop <- INPUT_pop[, .(apop = sum(pop)),
                              by = .(code)]
    # sort by admin code
    dt_admin_pop <- setorder(dt_admin_pop, code)
    # write "<grid_scale>_admin_population"
    fwrite(dt_admin_pop[, .(apop)],
           file=file.path(dir_out, paste0(grid_scale, "__admin_population.tsv")),
           col.names=FALSE)
    
    #-- compute grid ----------------------------------------------------------#
    
    # compute *projected* lower left corner of the grid
    g_origin <- INPUT_pop[ , c(lon = min(lon), lat = min(lat)) ]
    # grid_scale in grid metrics
    g_scale <- as.integer(grid_scale) / 120
    
    # shift grid by half the grid side in *one direction*
    if (country == "Kenya" & grid_scale == "20") {
      if (grepl("lon", grid_change))
        g_origin[1] <- g_origin[1] - g_scale/2
      if (grepl("lat", grid_change))
        g_origin[2] <- g_origin[2] - g_scale/2
    }
    
    # for each pixel, compute the coordinates of the corresponding grid square (lower left corner)
    INPUT_pop[, `:=`(glon = floor( (lon - g_origin[1]) / g_scale ) * g_scale + g_origin[1],
                     glat = floor( (lat - g_origin[2]) / g_scale ) * g_scale + g_origin[2])]
    
    #-- population output -----------------------------------------------------#
    
    # compute majority admin, population and population-weighted centroid by patch
    dt_patch_pop <- INPUT_pop[, .(gcode  = my.mode(code),
                                  gpop   = sum(pop),
                                  clon   = sum(pop * lon) / sum(pop),
                                  clat   = sum(pop * lat) / sum(pop)),
                              by = .(glon, glat)]
    # sort by grid admin code
    dt_patch_pop <- setorder(dt_patch_pop, gcode)
    # write "<grid_scale>_patch_population"
    fwrite(dt_patch_pop[, .(gpop)],
           file=file.path(dir_out, paste0(grid_scale, "__patch_population.tsv")),
           col.names=FALSE)
    # write "<grid_scale>_patches_per_admin"
    fwrite(dt_patch_pop[, .N, by=gcode][,.(N)],
           file=file.path(dir_out, paste0(grid_scale, "__patches_per_admin.tsv")),
           col.names=FALSE)
    # write "<grid_scale>_metadata"
    metadata <- data.table(names =c("nr_admins",                         "nr_patches"),
                           values=c(dt_patch_pop[,length(unique(gcode))], dt_patch_pop[,.N]))
    fwrite(metadata,
           file=file.path(dir_out, paste0(grid_scale, "__metadata.tsv")),
           sep="\t",
           col.names=FALSE)
    
    
    #-- distance output -------------------------------------------------------#
    
    # distance in m (shortest distance assuming earth is perfect sphere)
    distance <- distm(cbind(dt_patch_pop[,clon], dt_patch_pop[,clat]),
                      fun = distHaversine)
    # convert *upper diagonal* into integer vector
    distance <- as.integer(distance[lower.tri(distance)])
    # write "<grid_scale>_patch_distance"
    fwrite(as.data.table(distance),
           file = file.path(dir_out, paste0(grid_scale, "__patch_distance.tsv")),
           col.names=FALSE)
    
    
    #-- flow data output ------------------------------------------------------#
    # make Namibian flow data integer by multiplying each row (=origin) by the corresponding population size
    mat_patch_flow <- INPUT_flow
    if (country == "Namibia") {
      mat_patch_flow <- INPUT_flow * dt_admin_pop[, apop] + 0.5
      class(mat_patch_flow) <- "integer"
    }
    
    # write "<grid_scale>__admin_flow_data"
    fwrite(as.data.table(mat_patch_flow),
           file=file.path(dir_out, paste0(grid_scale, "__admin_flow_data.tsv")),
           sep="\t",
           col.names=FALSE)
    
    toc()
  } # end of loop over "grid_scales"
  
  
}
