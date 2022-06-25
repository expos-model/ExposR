# Copyright (C) President and Fellows of Harvard College 2021

# This program is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public
#   License along with this program.  If not, see
#   <http://www.gnu.org/licenses/>.

# EXPOS uses a digital elevation model (dem) to estimate exposed and
# protected areas for a given wind direction and inflection angle. The
# resulting topograhic exposure maps can be combined with output from 
# the Hurrecon model to estimate hurricane wind damage across a region.
# EXPOS contains two main functions:

# 1. The expos_model function estmates topopgrahic exposure for a specified
# wind direction and inflection angle. The input file is assumed to be a
# raster of elevation values in GeoTiff format with missing values represented
# by zero. Cells may be rectangular. If a geographic coordinate system is used,
# horizontal and vertical units are assumed to be degrees and meters, 
# respectively; otherwise horizontal and vertical units must be the same. 
# Columns are assumed to be closely aligned with true North (0 degrees); 
# if not, the map orientation (in degrees) must be specified. The name of 
# the input file is assumed to be "dem.tif".

# The output file is a raster file in GeoTiff format with the following
# values: 0 = missing data, 1 = protected, 2 = exposed. Output files
# are named "expos-xxx-yy.tif" where xxx is the wind direction and yy
# is the inflection angle.

# 2. The expos_damage function estimates regional hurricane damage as a 
# function of topographic exposure to peak wind direction at each raster cell.
# If a cell is protected, the enhanced Fujita scale (EF) rating is reduced
# by a specified amount. This function requires a hurricane GeoTiff file 
# created by Hurrecon, eight exposure files created by Expos (N, NE, E, etc), 
# and a reprojection file in csv format that contains lat/long coordinates 
# for the lower left and upper right corners of the digital elevation model.

# The output file is a raster file in GeoTiff format with the following values:
# 0 = missing, 1 = no damage, 2 = F0 damage, 3 = F1 damage, 4 = F2 damage, 
# 5 = F3 damage, 6 = F4 damage, 7 = F5 damage. Output files are named 
# "hhhh-damage-yy-z.tif" where hhhh is the hurricane ID, yy is the inflection 
# angle, and z is the reduction in EF rating for protected areas.

# Emery R. Boose
# June 2022

# R version 4.2.0

# Required packages:
#  raster

### INTERNAL FUNCTIONS ####################################

# create environment to store path

exp_env <- new.env(parent=emptyenv())

#' get_path returns the path for the current set of model runs.
#' If not set, a message is displayed.
#' @return current path
#' @noRd

get_path <- function() {
    # display message if not set
    if (!exists("exp_path", envir=exp_env)) {
        stop("Path not set. Please use expos_set_path.", call. = FALSE)

    # otherwise return current path
    } else {
        exp_path <- exp_env[["exp_path"]]
        invisible(exp_path)
    }
}

#' check_file_exists displays an error message and stops execution if
#' the specified file does not exist.
#' @param file_name name of file
#' @return no return value
#' @noRd

check_file_exists <- function(file_name) {
    if (file.exists(file_name) == FALSE) {
        stop("File not found: ", file_name)
    }
}

#' get_subdirectory returns the subdirectory for a given filename
#' and stops execution if the filename is not supported.
#' @param filename name of file
#' @return subdirectory
#' @noRd

get_subdirectory <- function(filename) {
    if (filename == "dem") {
        subdir <- "/dem/"
    } else if (grepl("expos", filename)) {
        subdir <- "/exposure/"
    } else if (grepl("damage", filename)) {
        subdir <- "/damage/"
    } else {
        stop("File name not supported")
    }

    return(subdir)
}

#' get_row_order returns TRUE if the row order remains unchanged 
#' (quadrants I-II) or FALSE if the row order is reversed
#' (quadrants III-IV).
#' @param wind_direction wind direction (degrees)
#' @return row order (TRUE or FALSE)
#' @noRd

get_row_order <- function(wind_direction) {
    if (wind_direction > 90 && wind_direction < 270) {
        row_order <- FALSE

    } else {
        row_order <- TRUE
    }

    return(row_order)
}

#' get_col_order returns TRUE if the column order remains unchanged
#' (quadrants I, IV) or FALSE if the column order is reversed
#' (quadrants II-III).
#' @param wind_direction wind direction (degrees)
#' @return column order (TRUE or FALSE)
#' @noRd

get_col_order <- function(wind_direction) {
    if (wind_direction > 180 && wind_direction < 360) {
        col_order <- TRUE

    } else {
        col_order <- FALSE
    }

    return(col_order)
}

#' get_transposed_wind_direction returns the transposed wind direction
#' after the actual wind direction is shifted to quadrant II.
#' @param wdir wind direction (degrees)
#' @return transposed wind direction (degrees)
#' @noRd

get_transposed_wind_direction <- function(wdir) {
    # quadrant I
    if (wdir >= 0 && wdir <= 90) {
        t_dir <- 360 - wdir

    # quadrant IV
    } else if (wdir > 90 && wdir < 180) {
        t_dir <- wdir + 180

    # quadrant III
    } else if (wdir >= 180 && wdir < 270) {
        t_dir <- 540 - wdir;

    # quadrant II
    } else if (wdir >= 270 && wdir <= 360) {
        t_dir  <- wdir;
    }

    return(t_dir)
}

#' west_north_west creates and saves a raster of exposure values for
#' transposed wind directions between 270 degrees and the cell diagonal 
#' (WNW). The transposed matrix of elevation values is processed in column
#' major order.
#' @param wind_direction wind direction (degrees)
#' @param inflection_angle inflection angle (degrees)
#' @param t_dir transposed wind direction (degrees)
#' @param lat_long whether coordinate system is latitude/longitude (degrees)
#' @return raster of modeled exposure values
#' @noRd

west_north_west <- function(wind_direction, inflection_angle, t_dir, lat_long) {
    # get path
    exp_path <- get_path()
    
    # convert 1 degree of latitude to meters
    deg2meters <- 111195
 
    # read dem file in GeoTiff format
    dem_file <- paste(exp_path, "/dem/dem.tif", sep="")
    check_file_exists(dem_file)
    dem_r <- raster::raster(dem_file)
  
    # get number of rows & columns
    nrows <- dim(dem_r)[1]
    ncols <- dim(dem_r)[2]

    # get extent
    xmn <- raster::extent(dem_r)[1]
    xmx <- raster::extent(dem_r)[2]
    ymn <- raster::extent(dem_r)[3]
    ymx <- raster::extent(dem_r)[4]
  
    # calculate cell dimensions
    cell_x <- (xmx-xmn)/ncols
    cell_y <- (ymx-ymn)/nrows

    # adjust if lat/long
    if (lat_long == TRUE) {
        lat_mid <- (xmx-xmn)/2
        cell_x <- cell_x*deg2meters*cos(lat_mid*pi/180)
        cell_y <- cell_y*deg2meters
    }

    # set exposure values
    pro_value <- 1
    exp_value <- 2
    
    # get row & column order
    row_order <- get_row_order(wind_direction)
    col_order <- get_col_order(wind_direction)

    # flip raster as needed
    if (row_order == TRUE && col_order == TRUE) {
        rr <- dem_r

    } else if (row_order == FALSE && col_order == TRUE) {
        rr <- raster::flip(dem_r, direction="y")

    } else if (row_order == TRUE && col_order == FALSE) {
        rr <- raster::flip(dem_r, direction="x")

    } else if (row_order == FALSE && col_order == FALSE) {
        xx  <- raster::flip(dem_r, direction="y")
        rr <- raster::flip(xx, direction="x")
    }

    # create dem matrix
    dem_m <- raster::as.matrix(rr)

    # create exposure array
    expos_m <- matrix(0, nrows, ncols)

    # create vectors for intermediate values
    p_shift <- vector(mode="numeric", length=ncols)  # shift value for each column
    h_pos   <- vector(mode="numeric", length=nrows)  # height of land column number
    h_elev  <- vector(mode="numeric", length=nrows)  # height of land elevation
  
    # get tangent of inflection angle
    tan_inf <- tan(inflection_angle*pi/180)

    # get adjustment for transposed wind direction
    adj <- tan((t_dir - 270)*pi/180)
  
    # get shift value for each column
    for (j in 1:ncols) {
        p_shift[j] <- round(adj*(j-1)*cell_x/cell_y)
    }

    # calculate exposure values
    for (j in 1:ncols) {    
        # display every 10th col number
        if (j %% 10 == 0) {
            cat("\rcol", j)
        }

        # first column exposed by default
        if (j == 1) {
            for (i in 1:nrows) {
                h_pos[i] <- 1
                h_elev[i] <- dem_m[i, 1]
                expos_m[i, 1] <- exp_value
            }
  
        } else {
            # shift for current column (0 or 1)
            shift <- p_shift[j] - p_shift[j-1]
 
            # shift by one row
            if (shift == 1) {
                for (i in (nrows-1):1) {
                    h_pos[i+1] <- h_pos[i]
                    h_elev[i+1] <- h_elev[i]
                }

                h_pos[1] <- j
                h_elev[1] <- dem_m[1, j]
                expos_m[i, j] <- exp_value
            }

            for (i in (shift+1):nrows) {
                # exposed (higher elevation)
                if (dem_m[i, j] >= h_elev[i]) {
                    h_pos[i] <- j
                    h_elev[i] <- dem_m[i, j]
                    expos_m[i, j] <- exp_value
        
                } else {
                    x_dist <- (p_shift[j] - p_shift[h_pos[i]])*cell_x
                    y_dist <- (j - h_pos[i])*cell_y
                    xy_dist <- sqrt(x_dist^2 + y_dist^2)
                    z_dist <- xy_dist * tan_inf
          
                    # exposed (beyond wind shadow)
                    if (dem_m[i, j] >= h_elev[i] - z_dist) {
                        h_pos[i] <- j
                        h_elev[i] <- dem_m[i, j]
                        expos_m[i, j] <- exp_value
        # 
                    # protected (in wind shadow)
                    } else {
                        expos_m[i, j] <- pro_value
                    }
                }
            }
        }
    }
    
    # set missing values to zero
    mask <- dem_m
    mask[mask != 0] <- 1
    zz <- mask * expos_m

    # create raster
    rr <- raster::raster(nrows=nrows, ncols=ncols, xmn=xmn, xmx=xmx, 
        ymn=ymn, ymx=ymx, vals=zz)

    # flip raster as needed
    if (row_order == TRUE && col_order == TRUE) {
        expos_r <- rr

    } else if (row_order == FALSE && col_order == TRUE) {
        expos_r <- raster::flip(rr, direction="y")

    } else if (row_order == TRUE && col_order == FALSE) {
        expos_r <- raster::flip(rr, direction="x")

    } else if (row_order == FALSE && col_order == FALSE) {
        xx  <- raster::flip(rr, direction="y")
        expos_r <- raster::flip(xx, direction="x")
    }

    # copy coordinate reference system from dem
    raster::crs(expos_r) <- raster::crs(dem_r)

    # return modeled values as raster
    invisible(expos_r)
}

#' north_north_west creates and saves a raster of exposure values for
#' transposed wind directions between the cell diagonal and 360 degrees 
#' (NNW). The transposed matrix of elevation values is processed in row
#' major order.
#' @param wind_direction wind direction (degrees)
#' @param inflection_angle inflection angle (degrees)
#' @param t_dir transposed wind direction (degrees)
#' @param lat_long whether coordinate system is latitude/longitude (degrees)
#' @return raster of modeled exposure values
#' @noRd

north_north_west <- function(wind_direction, inflection_angle, t_dir, lat_long) {
    # get path
    exp_path <- get_path()
    
    # convert 1 degree of latitude to meters
    deg2meters <- 111195

    # read dem file in GeoTiff format
    dem_file <- paste(exp_path, "/dem/dem.tif", sep="")
    check_file_exists(dem_file)
    dem_r <- raster::raster(dem_file)
  
    # get number of rows & columns
    nrows <- dim(dem_r)[1]
    ncols <- dim(dem_r)[2]

    # get extent
    xmn <- raster::extent(dem_r)[1]
    xmx <- raster::extent(dem_r)[2]
    ymn <- raster::extent(dem_r)[3]
    ymx <- raster::extent(dem_r)[4]
  
    # calculate cell dimensions
    cell_x <- (xmx-xmn)/ncols
    cell_y <- (ymx-ymn)/nrows

    # adjust if lat/long
    if (lat_long == TRUE) {
        lat_mid <- (xmx-xmn)/2
        cell_x <- cell_x*deg2meters*cos(lat_mid*pi/180)
        cell_y <- cell_y*deg2meters
    }

    # set exposure values
    pro_value <- 1
    exp_value <- 2

    # get row & column order
    row_order <- get_row_order(wind_direction)
    col_order <- get_col_order(wind_direction)

    # flip raster as needed
    if (row_order == TRUE && col_order == TRUE) {
        rr <- dem_r

    } else if (row_order == FALSE && col_order == TRUE) {
        rr <- raster::flip(dem_r, direction="y")

    } else if (row_order == TRUE && col_order == FALSE) {
        rr <- raster::flip(dem_r, direction="x")

    } else if (row_order == FALSE && col_order == FALSE) {
        xx  <- raster::flip(dem_r, direction="y")
        rr <- raster::flip(xx, direction="x")
    }

    # create dem matrix
    dem_m <- raster::as.matrix(rr)
  
    # create exposure array
    expos_m <- matrix(0, nrows, ncols)

    # create vectors for intermediate values
    p_shift <- vector(mode="numeric", length=nrows)  # shift value for each row
    h_pos   <- vector(mode="numeric", length=ncols)  # height of land row number
    h_elev  <- vector(mode="numeric", length=ncols)  # height of land elevation
  
    # get tangent of inflection angle
    tan_inf <- tan(inflection_angle*pi/180)

    # get adjustment for transposed wind direction
    adj <- tan((360 - t_dir)*pi/180)
  
    # get shift value for each row
    for (i in 1:nrows) {
        p_shift[i] <- round(adj*(i-1)*cell_y/cell_x)
    }

    # calculate exposure values
    for (i in 1:nrows) {
        # display every 10th row number
        if (i %% 10 == 0) {
            cat("\rrow", i)
        }

        # first row is exposed by default
        if (i == 1) {
            for (j in 1:ncols) {
                h_pos[j] <- 1
                h_elev[j] <- dem_m[1, j]
                expos_m[1, j] <- exp_value
            }
        } else {
            # shift for current row (0 or 1)
            shift <- p_shift[i] - p_shift[i-1];

            # shift by one column
            if (shift == 1) {
                for (j in (ncols-1):1) {
                    h_pos[j+1] <- h_pos[j]
                    h_elev[j+1] <- h_elev[j]
                }
      
                h_pos[1] <- i
                h_elev[1] <- dem_m[i, 1]
                expos_m[i, j] <- exp_value
            }

            for (j in (shift+1):ncols) {
                # exposed (higher elevation)
                if (dem_m[i, j] >= h_elev[j]) {
                    h_pos[j] <- i
                    h_elev[j] <- dem_m[i, j]
                    expos_m[i, j] <- exp_value
        
                } else {
                    x_dist <- (p_shift[i] - p_shift[h_pos[j]])*cell_x
                    y_dist <- (i - h_pos[j])*cell_y
                    xy_dist <- sqrt(x_dist^2 + y_dist^2)
                    z_dist <- xy_dist * tan_inf
          
                    # exposed (beyond wind shadow)
                    if (dem_m[i, j] >= h_elev[j] - z_dist) {
                        h_pos[j] <- i
                        h_elev[j] <- dem_m[i, j]
                        expos_m[i, j] <- exp_value
          
                    # protected (in wind shadow)
                    } else {
                        expos_m[i, j] <- pro_value
                    }
                }
            }
        }
    }
  
    # set missing values to zero
    mask <- dem_m
    mask[mask != 0] <- 1
    zz <- mask * expos_m

    # create raster
    rr <- raster::raster(nrows=nrows, ncols=ncols, xmn=xmn, xmx=xmx, 
        ymn=ymn, ymx=ymx, vals=zz)

    # flip raster as needed
    if (row_order == TRUE && col_order == TRUE) {
        expos_r <- rr

    } else if (row_order == FALSE && col_order == TRUE) {
        expos_r <- raster::flip(rr, direction="y")

    } else if (row_order == TRUE && col_order == FALSE) {
        expos_r <- raster::flip(rr, direction="x")

    } else if (row_order == FALSE && col_order == FALSE) {
        xx  <- raster::flip(rr, direction="y")
        expos_r <- raster::flip(xx, direction="x")
    }

    # copy coordinate reference system from dem
    raster::crs(expos_r) <- raster::crs(dem_r)

    # return modeled values as raster
    invisible(expos_r)
}


### UTILITY FUNCTIONS #####################################

#' @title
#' Utility Functions
#' @description
#' expos_set_path sets the path for the current set of model runs.
#' @param exp_path path for current model runs
#' @param console whether to display messages in console
#' @return no return value
#' @export
#' @rdname utility

expos_set_path <- function(exp_path, console=TRUE) {
    if (exp_path == "") {
        stop("Need to enter a path")

    } else if (dir.exists(exp_path) == FALSE) {
        stop("Path does not exist")
    }

    exp_env[["exp_path"]] <- exp_path

    if (console == TRUE) {
        cat("Path set to", exp_path, "\n")
    }
}

#' @description
#' expos_get_path returns the current path for a set of model runs.
#' @param console whether to display messages in console
#' @return current path
#' @export
#' @rdname utility

expos_get_path <- function(console=TRUE) {
    if (exists("exp_path", envir=exp_env)) {
        exp_path <- exp_env[["exp_path"]]

        if (console == TRUE) {
            cat(exp_path, "\n")
        }

        invisible(exp_path)

    } else {
        if (console == TRUE) {
            cat("Path not set\n")
        }

        invisible(NULL)
    }        
}


### MODELING FUNCTIONS ####################################

#' @title
#' Modeling Functions
#' @description
#' expos_model uses a raster file of elevation values, a specified wind
#' direction, and a specified inflection angle to create a raster file
#' of wind exposure values (0 = missing data, 1 = protected, 2 = exposed).
#' If a geographic coordinate system is used, horizontal and vertical units 
#' are assumed to be degrees and meters, respectively; otherwise horizontal 
#' and vertical units must be the same. Columns are assumed to be closely 
#' aligned with true North (0 degrees); if not, the map orientation must 
#' be specified. The name of the input file is assumed to be "dem.tif".

#' @param wind_direction wind direction (degrees)
#' @param inflection_angle inflection angle (degrees)
#' @param lat_long whether coordinate system is latitude/longitude (degrees)
#' @param orient map orientation (degrees)
#' @param save whether to save results to file
#' @param console whether to display messages in console
#' @return raster of modeled exposure values
#' @export
#' @rdname modeling

expos_model <- function(wind_direction, inflection_angle, lat_long=FALSE, orient=0,
    save=TRUE, console=TRUE) {
    
    # get path
    exp_path <- get_path()

    # announcement
    if (console == TRUE) {
        cat("... Modeling exposure ...\n")
    }

    # convert 1 degree of latitude to meters
    deg2meters <- 111195

    # check wind direction
    if (wind_direction < 0 || wind_direction > 360) {
        stop("Please supply wind direction in range 0-360 degrees")
    }

    # check inflection angle
    if (inflection_angle < 0 || inflection_angle > 90) {
        stop("Please supply inflection angle in range 0-90 degrees")
    }

    # read dem file in GeoTiff format
    dem_path <- paste(exp_path, "/dem/dem.tif", sep="")
    check_file_exists(dem_path)
    dem_r <- raster::raster(dem_path)
 
    # get number of rows & columns
    nrows <- dim(dem_r)[1]
    ncols <- dim(dem_r)[2]

    # get extent
    xmn <- raster::extent(dem_r)[1]
    xmx <- raster::extent(dem_r)[2]
    ymn <- raster::extent(dem_r)[3]
    ymx <- raster::extent(dem_r)[4]
  
    # calculate cell dimensions
    cell_x <- (xmx-xmn)/ncols
    cell_y <- (ymx-ymn)/nrows

    # adjust if lat/long
    if (lat_long == TRUE) {
        lat_mid <- (xmx-xmn)/2
        cell_x <- cell_x*deg2meters*cos(lat_mid*pi/180)
        cell_y <- cell_y*deg2meters
    }

    # get angle of cell diagonal
    cell_diagonal <- 360 - 180*atan(cell_x/cell_y)/pi;

    # adjust wind direction for map orientation
    wind_direction <- wind_direction + orient

    if (wind_direction > 360) {
        wind_direction <- wind_direction - 360
    }

    # get transposed wind direction
    t_dir <- get_transposed_wind_direction(wind_direction)
  
    # create exposure map
    if (t_dir < cell_diagonal) {
        expos_r <- west_north_west(wind_direction, inflection_angle, t_dir, lat_long)

    } else {
        expos_r <- north_north_west(wind_direction, inflection_angle, t_dir, lat_long)
    }

    # output
    if (save == TRUE) {
        # save modeled values in a Geotiff file
        expos_file = paste(exp_path, "/exposure/expos-", formatC(wind_direction, width=3, flag="0"), "-", 
            formatC(inflection_angle, width=2, flag="0"), ".tif", sep="")

        rgdal::setCPLConfigOption("GDAL_PAM_ENABLED", "FALSE")
        raster::writeRaster(expos_r, expos_file, overwrite=TRUE)
    
        if (console == TRUE) {
            cat("\nSaving to", expos_file, "\n")
        }
    }

    # return modeled values as raster
    invisible(expos_r)
}

#' @description
#' expos_damage uses output from Hurrecon and Expos to create a raster
#' of hurricane wind damage where topograhic exposure at each location
#' is determined by peak wind direction. If a location is protected, 
#' the enhanced Fujita scale rating is reduced by a specified amount.
#' This function requires a hurricane tif file created by Hurrecon, 
#' eight exposure files created by Expos (N, NE, E, etc), and a reprojection
#' file in csv format that contains lat long coordinates for the lower left 
#' and upper right corners of the digital elevation model.
#' @param hurricane hurricane name (as it appears in tif file)
#' @param inflection_angle inflection angle (degrees)
#' @param protect how much to reduce damage in protected areas (Fujita 
#' scale ratings)
#' @param save whether to save results to file
#' @param console whether to display messages in console
#' @return raster of modeled wind damage values
#' @export
#' @rdname modeling

expos_damage <- function(hurricane, inflection_angle, protect, save=TRUE, 
    console=TRUE) {
    
    # get path
    exp_path <- get_path()

    # announcement
    if (console == TRUE) {
        cat("... Modeling damage ...\n")
    }

    # check protect value
    if (protect < 0 || protect > 6) {
        stop("Please supply protect in range 0-6")
    }

    # read exposure files
    ee_r <- list()

    for (i in 1:8) {
        wind_direction <- (i-1)*45

        expos_file <- paste(exp_path, "/exposure/expos-", formatC(wind_direction, width=3, flag="0"), "-", 
            formatC(inflection_angle, width=2, flag="0"), ".tif", sep="")
    
        ee_r[[i]] <- raster::raster(expos_file)
    }

    # read dem file
    dem_file <- paste(exp_path, "/dem/dem.tif", sep="")
    dem_r <- raster::raster(dem_file)

    dem_rows <- dim(dem_r)[1]
    dem_cols <- dim(dem_r)[2]

    dem_xmn <- raster::extent(dem_r)[1]
    dem_xmx <- raster::extent(dem_r)[2]
    dem_ymn <- raster::extent(dem_r)[3]
    dem_ymx <- raster::extent(dem_r)[4]

    # read hurrecon file
    hur_file <- paste(exp_path, "/damage/", hurricane, ".tif", sep="")
    ff_r <- raster::raster(hur_file, 2)  # enhanced Fujita scale
    cc_r <- raster::raster(hur_file, 4)  # cardinal wind direction (1-8)

    hur_rows <- dim(ff_r)[1]
    hur_cols <- dim(ff_r)[2]

    hur_xmn <- raster::extent(ff_r)[1]
    hur_xmx <- raster::extent(ff_r)[2]
    hur_ymn <- raster::extent(ff_r)[3]
    hur_ymx <- raster::extent(ff_r)[4]

    # read reproject file
    reproject_file <- paste(exp_path, "/damage/reproject.csv", sep="")
    rp <- read.csv(reproject_file, header=TRUE)

    lat_0 <- rp$lat_0
    lon_0 <- rp$lon_0
    lat_1 <- rp$lat_1
    lon_1 <- rp$lon_1

    # convert rasters to matrices
    ee_m <- list()

    for (i in 1:8) {
        ee_m[[i]] <- raster::as.matrix(ee_r[[i]])
    }

    dem_m <- raster::as.matrix(dem_r)
    ff_m  <- raster::as.matrix(ff_r)
    cc_m  <- raster::as.matrix(cc_r)

    # create damage matrix (0 = missing, 1 = no damage)
    dam_m <- dem_m
    dam_m[dam_m != 0] <- 1

    # calculate damage values
    for (i in 1:dem_rows) {
        if (console == TRUE) {
            # display every 10th row number
            if (i %% 10 == 0) {
                cat("\rrow", i)
            }
        }

        for (j in 1:dem_cols) {
            # skip missing values in dem
            if (dem_m[i, j] != 0) {
                # get lat long coordinates
                hur_x <- lon_0 + (lon_1 - lon_0)*(j - 0.5)/dem_cols
                hur_y <- lat_1 - (lat_1 - lat_0)*(i - 0.5)/dem_rows

                # get row & col in hurricane file (note: row 1 = top of raster)
                hur_row <- ceiling(hur_rows - hur_rows*(hur_y - hur_ymn)/(hur_ymx - hur_ymn))
                hur_col <- ceiling(hur_cols*(hur_x - hur_xmn)/(hur_xmx - hur_xmn))

                # # check if in Hurrecon output
                if (hur_row >= 1 && hur_row <= hur_rows && hur_col >= 1 && hur_col <= hur_cols) {
                    # get peak wind direction (1-8)
                    pdir <- cc_m[hur_row, hur_col]
                   
                    if (pdir == 0) {
                        # no damage if no peak wind direction
                        dam_m[i, j] <- 1

                    } else {
                        # get topographic exposure
                        exposure <- ee_m[[pdir]][i, j]
                        # get fujita scale damage (0-7)
                        damage <- ff_m[hur_row, hur_col]

                        # protected
                        if (exposure == 1) {
                            # reduce damage
                            pro_damage <- damage - protect
                            
                            if (pro_damage < 1) {
                                pro_damage <- 1
                            }   

                            dam_m[i, j] <- pro_damage

                        # exposed
                        } else {
                            dam_m[i, j] <- damage
                        }
                    }

                # otherwise set to missing
                } else {
                   dam_m[i, j] <- 0
                }
            }
        }
    }

    # create raster of modeled results
    dam_r <- raster::raster(nrows=dem_rows, ncols=dem_cols, xmn=dem_xmn, xmx=dem_xmx, 
        ymn=dem_ymn, ymx=dem_ymx, vals=dam_m)

    # copy coordinate reference system from dem
    raster::crs(dam_r) <- raster::crs(dem_r)

    if (save == TRUE) {
        # save modeled results in GeoTiff file
        dam_file <- paste(exp_path, "/damage/", hurricane, "-damage-", 
            formatC(inflection_angle, width=2, flag="0"), "-", protect, ".tif", sep="")
        
        rgdal::setCPLConfigOption("GDAL_PAM_ENABLED", "FALSE")
        raster::writeRaster(dam_r, dam_file, overwrite=TRUE)

        if (console == TRUE) {
            cat("\nSaving to", dam_file, "\n")
        }
    }

    # return modeled results as raster
    invisible(dam_r)
}


### SUMMARIZING FUNCTIONS #################################

#' @title
#' Summarizing Functions
#' @description
#' expos_summarize displays summary information for a specified raster
#' file, including the number of rows and columns, spatial extent, cell
#' height and width, and minimum and maximum value.
#' @param filename name of input raster file
#' @param console whether to display results in console
#' @return a string containing summary information
#' @export
#' @rdname summarizing

expos_summarize <- function(filename, console=TRUE) {
    # get path
    exp_path <- get_path()

    # announcement
    if (console == TRUE) {
        cat("... Summarizing raster ...\n")
    }
    
    # get subdirectory
    subdir <- get_subdirectory(filename)

    # read file in GeoTiff format
    file_path <- paste(exp_path, subdir, filename, ".tif", sep="")
    check_file_exists(file_path)
    rr <- raster::raster(file_path)

    # get number of rows & columns
    nrows <- dim(rr)[1]
    ncols <- dim(rr)[2]

    # get extent
    xmn <- raster::extent(rr)[1]
    xmx <- raster::extent(rr)[2]
    ymn <- raster::extent(rr)[3]
    ymx <- raster::extent(rr)[4]
  
    # calculate cell dimensions
    cell_x <- (xmx-xmn)/ncols
    cell_y <- (ymx-ymn)/nrows

    # get min & max values
    val_min <- raster::minValue(rr)
    val_max <- raster::maxValue(rr)

    # create display string
    st <- paste("Rows: ", nrows, "  Columns: ", ncols, "\n", sep="")
    st <- paste(st, "Northing: ", round(ymn), " to ", round(ymx), "\n", sep="")
    st <- paste(st, "Easting: ", round(xmn), " to ", round(xmx), "\n", sep="")
    st <- paste(st, "Cell height: ", round(cell_y), "\n", sep="")
    st <- paste(st, "Cell width: ", round(cell_x), "\n", sep="")
    st <- paste(st, "Values: ", round(val_min, 1), " to ", round(val_max, 1), "\n", sep="")
    
    # display results in console
    if (console == TRUE) {
        cat(st)
    }

    invisible(st)
}


### PLOTTING FUNCTIONS ####################################

#' @title
#' Plotting Functions
#' @description
#' expos_plot creates a plot of a specified raster file. Optional arguments
#' include plot title, horizontal units, vertical units, vector (boundary
#' files) and color palette.
#' @param filename name of input raster file
#' @param title plot title
#' @param h_units horizontal units
#' @param v_units vertical units
#' @param vector whether to display vectory boundary files
#' @param colormap color palette
#' @return no return value
#' @export
#' @rdname plotting

expos_plot <- function(filename, title="", h_units="meters", v_units="meters",
    vector=TRUE, colormap="default", console=TRUE) {
    
    # get path
    exp_path <- get_path()

    if (console == TRUE) {
        cat("... Plotting raster ...\n")
    }

    # get subdirectory
    subdir <- get_subdirectory(filename)

    # read file in GeoTiff format
    file_path <- paste(exp_path, subdir, filename, ".tif", sep="")
    check_file_exists(file_path)
    rr <- raster::raster(file_path)

    # get vector boundary file
    if (vector == TRUE) {
        boundaries_file <- paste(exp_path, "/vector/boundaries.shp", sep="")
        check_file_exists(boundaries_file)
        boundaries <- rgdal::readOGR(boundaries_file)
    }

    rr_min <- raster::minValue(rr)
    rr_max <- raster::maxValue(rr)

    # default palettes
    if (length(colormap) == 1) {
        if (colormap == "default") {
            # exposure map
            if (grepl("expos", filename)) {
                cmap <- c("white", "grey", "blue")

            # damage map
            } else if (grepl("damage", filename)) {
                all_cols <- c("white", "grey", "purple", "blue", "green", "yellow", "orange", "red")
                cmap <- c(all_cols[rr_min+1])

                if (rr_max > rr_min) {
                    for (i in (rr_min+2):(rr_max+1)) {
                        cmap <- append(cmap, all_cols[i])
                    }
                }

            # dem, etc
            } else {
                cmap <- rev(terrain.colors(255))
            }
        
        # raster palette
        } else if (colormap == "r_default") {
            cmap <- rev(terrain.colors(255))
        }

    # user-specified palette
    } else {
        cmap <- colormap
    }

    # create plot
    par(mar=c(5.1, 4.6, 4.1, 2.1))

    if (grepl("dem", filename)) {
        if (title == "") {
            title <- "Elevation"
        }
        v_units_str <- paste("  ", v_units, sep="")
        raster::plot(rr, main=title, xlab=h_units, ylab=h_units,
            legend.args=list(text=v_units_str, line=1), col=cmap)
        if (vector == TRUE) {
            raster::plot(boundaries, add=TRUE)
        }
    
    } else if (grepl("expos", filename)) {
        if (title == "") {
            x <- strsplit(filename, "-")[[1]]
            title <- paste("Exposure ", x[2], "-", x[3], sep="")
        }
        vals <- c(0, 1, 2)
        labs <- c("", "Pro", "Exp")
        arg <- list(at=vals, labels=labs)
        raster::plot(rr, main=title, xlab=h_units, ylab=h_units, axis.args=arg, 
            legend.args=list(text='  Exposure', line=1), col=cmap)
        if (vector == TRUE) {
            raster::plot(boundaries, add=TRUE)
        }

    } else if (grepl("damage", filename)) {
        if (title == "") {
            x <- strsplit(filename, "-")[[1]]
            title <- paste(x[1], "-", x[2], " Damage ", x[4], "-", x[5],sep="")
        }
        vals <- c(0, 1, 2, 3, 4, 5, 6, 7)
        labs <- c("", "None", "EF0", "EF1", "EF2", "EF3", "EF4", "EF5")
        raster::plot(rr, main=title, xlab=h_units, ylab=h_units,
            axis.args=list(at=vals, labels=labs), 
            legend.args=list(text='  EF Scale', line=1), col=cmap)
        if (vector == TRUE) {
            raster::plot(boundaries, add=TRUE)
        }
    
    } else {
        if (title == "") {
            title <- filename
        }
        raster::plot(rr, main=title, xlab=h_units, ylab=h_units, col=cmap)
        if (vector == TRUE) {
            raster::plot(boundaries, add=TRUE)
        }
    }
}

