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

###############################################################################

# The EXPOS model uses a digital elevation model (DEM) to estimate exposed
# and protected areas for a given hurricane wind direction and inflection angle.
# The resulting topograhic exposure maps can be combined with output from the 
# HURRECON model to estimate hurricane wind damage across a region.

### INTERNAL FUNCTIONS ####################################

# create environment to store path

exp_env <- new.env(parent=emptyenv())

#' get_exp_path returns the path for the current set of model runs.
#' If not set, an error message is displayed.
#' @return current path
#' @noRd

get_exp_path <- function() {
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

#' get_subdirectory returns the subdirectory for a given file name
#' and stops execution if the file name is not supported.
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

#' check_raster_value checks to see if the specified raster layer
#' contains the specified value.
#' @param raster name of raster
#' @param layer number of layer
#' @param value value to check
#' @return whether raster layer contains value (TRUE or FALSE)
#' @noRd

check_raster_value <- function(raster, layer, value) {
    vv <- terra::unique(raster[[layer]])
    
    for (i in 1:nrow(vv)) {
        if (vv[i, 1] == value) {
            return(TRUE)
        }
    }

    return(FALSE)
}

#' check_lat_long tries to determine if the raster coordinate system
#' is latitude/longitude in degrees.  It returns TRUE if X values are 
#' between -180 and 180 and Y values are between -90 and 90; otherwise
#' it returns FALSE.
#' @param xmin minimum X value
#' @param xmax maximum X value
#' @param ymin minimum Y value
#' @param ymax maximum Y value
#' @return whether coordinates are lat/long (TRUE or FALSE)
#' @noRd

check_lat_long <- function(xmin, xmax, ymin, ymax) {
    if (xmin >= -180 && xmax <= 180 && ymin >= -90 && ymax <= 90) {
        return(TRUE)
    } else {
        return(FALSE)
    }
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

#' west_north_west creates and returns a raster of exposure values for
#' transposed wind directions between 270 degrees and the cell diagonal 
#' (WNW). The transposed matrix of elevation values is processed in column
#' major order.
#' @param wind_direction wind direction (degrees)
#' @param inflection_angle inflection angle (degrees)
#' @param t_dir transposed wind direction (degrees)
#' @param lat_long whether coordinate system is latitude/longitude
#' @return raster of modeled exposure values
#' @noRd

west_north_west <- function(wind_direction, inflection_angle, t_dir, lat_long) {
    # get path
    exp_path <- get_exp_path()
    
    # convert 1 degree of latitude to meters
    deg2meters <- 111195
 
    # read DEM file in GeoTiff format
    dem_file <- paste(exp_path, "/dem/dem.tif", sep="")
    check_file_exists(dem_file)
    dem_r <- terra::rast(dem_file)
  
    # get number of rows & columns
    nrows <- terra::nrow(dem_r)
    ncols <- terra::ncol(dem_r)

    # get extent
    xmin <- terra::ext(dem_r)[1]
    xmax <- terra::ext(dem_r)[2]
    ymin <- terra::ext(dem_r)[3]
    ymax <- terra::ext(dem_r)[4]
  
    # calculate cell dimensions
    cell_x <- (xmax-xmin)/ncols
    cell_y <- (ymax-ymin)/nrows

    # adjust if lat/long
    if (lat_long == TRUE) {
        lat_mid <- ymin + (ymax-ymin)/2
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
        rr <- terra::flip(dem_r, direction="vertical")

    } else if (row_order == TRUE && col_order == FALSE) {
        rr <- terra::flip(dem_r, direction="horizontal")

    } else if (row_order == FALSE && col_order == FALSE) {
        xx  <- terra::flip(dem_r, direction="vertical")
        rr <- terra::flip(xx, direction="horizontal")
    }

    # create dem matrix
    dem_m <- terra::as.matrix(rr, wide=TRUE)

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
            message(paste("\rcol", j), appendLF=FALSE)
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
    rr <- terra::rast(nrows=nrows, ncols=ncols, xmin=xmin, xmax=xmax, 
        ymin=ymin, ymax=ymax, vals=zz)

    # flip raster as needed
    if (row_order == TRUE && col_order == TRUE) {
        expos_r <- rr

    } else if (row_order == FALSE && col_order == TRUE) {
        expos_r <- terra::flip(rr, direction="vertical")

    } else if (row_order == TRUE && col_order == FALSE) {
        expos_r <- terra::flip(rr, direction="horizontal")

    } else if (row_order == FALSE && col_order == FALSE) {
        xx  <- terra::flip(rr, direction="vertical")
        expos_r <- terra::flip(xx, direction="horizontal")
    }

    # copy coordinate reference system from dem
    terra::crs(expos_r) <- terra::crs(dem_r)

    # return modeled values as raster
    invisible(expos_r)
}

#' north_north_west creates and returns a raster of exposure values for
#' transposed wind directions between the cell diagonal and 360 degrees 
#' (NNW). The transposed matrix of elevation values is processed in row
#' major order.
#' @param wind_direction wind direction (degrees)
#' @param inflection_angle inflection angle (degrees)
#' @param t_dir transposed wind direction (degrees)
#' @param lat_long whether coordinate system is latitude/longitude
#' @return raster of modeled exposure values
#' @noRd

north_north_west <- function(wind_direction, inflection_angle, t_dir, lat_long) {
    # get path
    exp_path <- get_exp_path()
    
    # convert 1 degree of latitude to meters
    deg2meters <- 111195

    # read dem file in GeoTiff format
    dem_file <- paste(exp_path, "/dem/dem.tif", sep="")
    check_file_exists(dem_file)
    dem_r <- terra::rast(dem_file)
  
    # get number of rows & columns
    nrows <- terra::nrow(dem_r)
    ncols <- terra::ncol(dem_r)

    # get extent
    xmin <- terra::ext(dem_r)[1]
    xmax <- terra::ext(dem_r)[2]
    ymin <- terra::ext(dem_r)[3]
    ymax <- terra::ext(dem_r)[4]
  
    # calculate cell dimensions
    cell_x <- (xmax-xmin)/ncols
    cell_y <- (ymax-ymin)/nrows

    # adjust if lat/long
    if (lat_long == TRUE) {
        lat_mid <- ymin + (ymax-ymin)/2
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
        rr <- terra::flip(dem_r, direction="vertical")

    } else if (row_order == TRUE && col_order == FALSE) {
        rr <- terra::flip(dem_r, direction="horizontal")

    } else if (row_order == FALSE && col_order == FALSE) {
        xx  <- terra::flip(dem_r, direction="vertical")
        rr <- terra::flip(xx, direction="horizontal")
    }

    # create dem matrix
    dem_m <- terra::as.matrix(rr, wide=TRUE)
  
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
            message(paste("\rrow", i), appendLF=FALSE)
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
    rr <- terra::rast(nrows=nrows, ncols=ncols, xmin=xmin, xmax=xmax, 
        ymin=ymin, ymax=ymax, vals=zz)

    # flip raster as needed
    if (row_order == TRUE && col_order == TRUE) {
        expos_r <- rr

    } else if (row_order == FALSE && col_order == TRUE) {
        expos_r <- terra::flip(rr, direction="vertical")

    } else if (row_order == TRUE && col_order == FALSE) {
        expos_r <- terra::flip(rr, direction="horizontal")

    } else if (row_order == FALSE && col_order == FALSE) {
        xx  <- terra::flip(rr, direction="vertical")
        expos_r <- terra::flip(xx, direction="horizontal")
    }

    # copy coordinate reference system from dem
    terra::crs(expos_r) <- terra::crs(dem_r)

    # return modeled values as raster
    invisible(expos_r)
}


### UTILITY FUNCTIONS #####################################

#' @title
#' Utility Functions
#' @description
#' expos_set_path sets the path for the current set of model runs.
#' @param exp_path path for current model runs
#' @return no return value
#' @export
#' @rdname utility

expos_set_path <- function(exp_path) {
    if (exp_path == "") {
        stop("Need to enter a path")

    } else if (dir.exists(exp_path) == FALSE) {
        stop("Path does not exist")
    }

    exp_env[["exp_path"]] <- exp_path
        message(paste("Path set to", exp_path))
}

#' @description
#' expos_get_path returns the current path for a set of model runs.
#' @return current path
#' @export
#' @rdname utility

expos_get_path <- function() {
    if (exists("exp_path", envir=exp_env)) {
        exp_path <- exp_env[["exp_path"]]

        message(exp_path)
        invisible(exp_path)

    } else {
        message("Path not set")
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
#' The user can specify if coordinates are lat/long; otherwise lat/long 
#' is assumed if X values are between -180 and 180 and Y values are between
#' -90 and 90. If lat/long, horizontal and vertical units are assumed 
#' to be degrees and meters, respectively; otherwise horizontal and vertical 
#' units must be the same. Columns are assumed to be closely aligned with 
#' true North (0 degrees); if not, the map orientation (azimuth) must be 
#' specified in degrees. The name of the input file is assumed to be "dem.tif".

#' @param wind_direction wind direction (degrees)
#' @param inflection_angle inflection angle (degrees)
#' @param lat_long whether coordinate system is latitude/longitude
#' @param orient map orientation (degrees)
#' @param save whether to save results to file
#' @param exp_path path for current set of model runs
#' @return raster of modeled exposure values
#' @export
#' @examples
#' exp_path <- system.file("", package="ExposR", mustWork=TRUE)
#' expos_model(wind_direction=135, inflection_angle=6, save=FALSE, exp_path=exp_path)
#' @rdname modeling

expos_model <- function(wind_direction, inflection_angle, lat_long=NULL, orient=0,
    save=TRUE, exp_path=NULL) {
    
    # get path
    if (!is.null(exp_path)) {
        expos_set_path(exp_path)
    } else {
        exp_path <- get_exp_path()
    }

    # announcement
    message("... Modeling exposure ...")

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
    dem_file <- paste(exp_path, "/dem/dem.tif", sep="")
    check_file_exists(dem_file)
    dem_r <- terra::rast(dem_file)
 
    # get number of rows & columns
    nrows <- terra::nrow(dem_r)
    ncols <- terra::ncol(dem_r)

    # get extent
    xmin <- terra::ext(dem_r)[1]
    xmax <- terra::ext(dem_r)[2]
    ymin <- terra::ext(dem_r)[3]
    ymax <- terra::ext(dem_r)[4]
 
    # calculate cell dimensions
    cell_x <- (xmax-xmin)/ncols
    cell_y <- (ymax-ymin)/nrows

    # check for lat/long
    if (is.null(lat_long)) {
        lat_long <- check_lat_long(xmin, xmax, ymin, ymax)
    }

    # adjust if lat/long
    if (lat_long == TRUE) {
        lat_mid <- ymin + (ymax-ymin)/2
        cell_x <- cell_x*deg2meters*cos(lat_mid*pi/180)
        cell_y <- cell_y*deg2meters
    }

    # get angle of cell diagonal
    cell_diagonal <- 360 - 180*atan(cell_x/cell_y)/pi;

    # adjust wind direction for map orientation
    wind_direction <- wind_direction - orient

    if (wind_direction < 0) {
        wind_direction <- wind_direction + 360
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

        terra::writeRaster(expos_r, expos_file, overwrite=TRUE, 
            gdal=c("GDAL_PAM_ENABLED", "FALSE"))
    
        message(paste("\nSaving to", expos_file))
    }

    # return modeled values as raster
    invisible(expos_r)
}

#' @description
#' expos_damage uses output from the EXPOS and HURRECON models to create 
#' a raster of hurricane wind damage where topographic exposure at each 
#' location is determined by peak wind direction. If a location is protected, 
#' the enhanced Fujita scale rating from HURRECON is reduced by a specified 
#' amount. This function requires a hurricane file in GeoTiff format created 
#' by HURRECON, exposure files created by EXPOS for the eight cardinal wind 
#' directions (N, NE, E, etc), and a reprojection file in CSV format 
#' (reproject.csv) that contains lat/long coordinates in degrees for the 
#' lower left and upper right corners of the digital elevation model.
#' @param hurricane hurricane name (as it appears in tif file)
#' @param inflection_angle inflection angle (degrees)
#' @param protect how much to reduce damage in protected areas (number of 
#' Fujita scale ratings)
#' @param save whether to save results to file
#' @param exp_path path for current set of model runs
#' @return raster of modeled wind damage values
#' @export
#' @rdname modeling

expos_damage <- function(hurricane, inflection_angle, protect, save=TRUE, 
    exp_path=NULL) {
    
    # get path
    if (!is.null(exp_path)) {
        expos_set_path(exp_path)
    } else {
        exp_path <- get_exp_path()
    }

    # announcement
    message("... Modeling damage ...")

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
    
        ee_r[[i]] <- terra::rast(expos_file)
    }

    # read dem file
    dem_file <- paste(exp_path, "/dem/dem.tif", sep="")
    dem_r <- terra::rast(dem_file)

    dem_rows <- terra::nrow(dem_r)
    dem_cols <- terra::ncol(dem_r)

    dem_xmin <- terra::ext(dem_r)[1]
    dem_xmax <- terra::ext(dem_r)[2]
    dem_ymin <- terra::ext(dem_r)[3]
    dem_ymax <- terra::ext(dem_r)[4]

    # read hurrecon file
    hur_file <- paste(exp_path, "/damage/", hurricane, ".tif", sep="")
    hur_r <- terra::rast(hur_file)
    ff_r <- hur_r[[2]]  # enhanced Fujita scale
    cc_r <- hur_r[[4]]  # cardinal wind direction (1-8)

    hur_rows <- terra::nrow(ff_r)
    hur_cols <- terra::ncol(ff_r)

    hur_xmin <- terra::ext(ff_r)[1]
    hur_xmax <- terra::ext(ff_r)[2]
    hur_ymin <- terra::ext(ff_r)[3]
    hur_ymax <- terra::ext(ff_r)[4]

    # read reproject file
    reproject_file <- paste(exp_path, "/damage/reproject.csv", sep="")
    rp <- utils::read.csv(reproject_file, header=TRUE)

    lat_0 <- rp$lat_0
    lon_0 <- rp$lon_0
    lat_1 <- rp$lat_1
    lon_1 <- rp$lon_1

    # convert rasters to matrices
    ee_m <- list()

    for (i in 1:8) {
        ee_m[[i]] <- terra::as.matrix(ee_r[[i]], wide=TRUE)
    }

    dem_m <- terra::as.matrix(dem_r, wide=TRUE)
    ff_m  <- terra::as.matrix(ff_r, wide=TRUE)
    cc_m  <- terra::as.matrix(cc_r, wide=TRUE)

    # create damage matrix (0 = missing, 1 = no damage)
    dam_m <- dem_m
    dam_m[dam_m != 0] <- 1

    # calculate damage values
    for (i in 1:dem_rows) {
        # display every 10th row number
        if (i %% 10 == 0) {
            message(paste("\rrow", i), appendLF=FALSE)
        }

        for (j in 1:dem_cols) {
            # skip missing values in dem
            if (dem_m[i, j] != 0) {
                # get lat long coordinates
                hur_x <- lon_0 + (lon_1 - lon_0)*(j - 0.5)/dem_cols
                hur_y <- lat_1 - (lat_1 - lat_0)*(i - 0.5)/dem_rows

                # get row & col in hurricane file (note: row 1 = top of raster)
                hur_row <- ceiling(hur_rows - hur_rows*(hur_y - hur_ymin)/(hur_ymax - hur_ymin))
                hur_col <- ceiling(hur_cols*(hur_x - hur_xmin)/(hur_xmax - hur_xmin))

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
    dam_r <- terra::rast(nrows=dem_rows, ncols=dem_cols, xmin=dem_xmin, xmax=dem_xmax, 
        ymin=dem_ymin, ymax=dem_ymax, vals=dam_m)

    # copy coordinate reference system from dem
    terra::crs(dam_r) <- terra::crs(dem_r)

    if (save == TRUE) {
        # save modeled results in GeoTiff file
        dam_file <- paste(exp_path, "/damage/", hurricane, "-damage-", 
            formatC(inflection_angle, width=2, flag="0"), "-", protect, ".tif", sep="")
        
        terra::writeRaster(dam_r, dam_file, overwrite=TRUE, 
            gdal=c("GDAL_PAM_ENABLED", "FALSE"))

        message(paste("\nSaving to", dam_file))
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
#' height and width, and minimum and maximum value. The user can specify 
#' if coordinates are lat/long; otherwise lat/long is assumed if X values 
#' are between -180 and 180 and Y values are between -90 and 90.
#' @param filename name of input raster file
#' @param lat_long whether coordinate system is latitude/longitude
#' @param console whether to display results in console
#' @param exp_path path for current set of model runs
#' @return a string containing summary information
#' @export
#' @rdname summarizing

expos_summarize <- function(filename, lat_long=NULL, console=TRUE, exp_path=NULL) {
    # get path
    if (!is.null(exp_path)) {
        expos_set_path(exp_path)
    } else {
        exp_path <- get_exp_path()
    }

    # announcement
    message("... Summarizing raster ...")
    
    # convert 1 degree of latitude to meters
    deg2meters <- 111195

    # get subdirectory
    subdir <- get_subdirectory(filename)

    # read file in GeoTiff format
    file_path <- paste(exp_path, subdir, filename, ".tif", sep="")
    check_file_exists(file_path)
    rr <- terra::rast(file_path)

    # get number of rows & columns
    nrows <- terra::nrow(rr)
    ncols <- terra::ncol(rr)

    # get extent
    xmin <- terra::ext(rr)[1]
    xmax <- terra::ext(rr)[2]
    ymin <- terra::ext(rr)[3]
    ymax <- terra::ext(rr)[4]
  
    # calculate cell dimensions
    cell_x <- (xmax-xmin)/ncols
    cell_y <- (ymax-ymin)/nrows

    # check for lat/long
    if (is.null(lat_long)) {
        lat_long <- check_lat_long(xmin, xmax, ymin, ymax)
     }

    b_units <- ""
    c_units <- ""
    v_units <- ""

    # adjust if lat/long
    if (lat_long == TRUE) {
        lat_mid <- ymin + (ymax-ymin)/2
        cell_x <- cell_x*deg2meters*cos(lat_mid*pi/180)
        cell_y <- cell_y*deg2meters

        b_units <- " degrees"
        c_units <- " meters"
        if (filename == "dem") {
            v_units <- " meters"
        }
    }

    # get min & max values
    val_min <- terra::minmax(rr)[1]
    val_max <- terra::minmax(rr)[2]

    # create display string
    st <- paste("Rows: ", nrows, "  Columns: ", ncols, "\n", sep="")
    st <- paste(st, "Northing: ", round(ymin, 6), " to ", round(ymax, 6), b_units, "\n", sep="")
    st <- paste(st, "Easting: ", round(xmin, 6), " to ", round(xmax, 6), b_units, "\n", sep="")
    st <- paste(st, "Cell height: ", round(cell_y), c_units, "\n", sep="")
    st <- paste(st, "Cell width: ", round(cell_x), c_units, "\n", sep="")
    st <- paste(st, "Values: ", round(val_min), " to ", round(val_max), v_units, "\n", sep="")
    
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
#' expos_plot creates a plot of a raster file. The user can specify if
#' coordinates are lat/long; otherwise lat/long is assumed if X values 
#' are between -180 and 180 and Y values are between -90 and 90. Optional 
#' arguments include plot title, horizontal units, vertical units, vector 
#' boundary files, and color palette. 
#' @param filename name of input raster file
#' @param title plot title
#' @param lat_long whether coordinate system is latitude/longitude
#' @param h_units horizontal units
#' @param v_units vertical units
#' @param vector whether to display vectory boundary files
#' @param colormap color palette
#' @param exp_path path for current set of model runs
#' @return no return value
#' @export
#' @rdname plotting

expos_plot <- function(filename, title="", lat_long=NULL, h_units="meters", 
    v_units="meters", vector=TRUE, colormap="default", exp_path=NULL) {
    
    # get path
    if (!is.null(exp_path)) {
        expos_set_path(exp_path)
    } else {
        exp_path <- get_exp_path()
    }

    message("... Plotting raster ...")

    # get subdirectory
    subdir <- get_subdirectory(filename)

    # read file in GeoTiff format
    file_path <- paste(exp_path, subdir, filename, ".tif", sep="")
    check_file_exists(file_path)
    rr <- terra::rast(file_path)

    # get vector boundary file
    if (vector == TRUE) {
        boundaries_file <- paste(exp_path, "/vector/boundaries.shp", sep="")
        check_file_exists(boundaries_file)
        boundaries <- terra::vect(boundaries_file)
    }

    rr_min <- terra::minmax(rr)[1]
    rr_max <- terra::minmax(rr)[2]
    rr_unique <- terra::unique(rr)

    # check for lat/long
    if (is.null(lat_long)) {
        xmin <- terra::ext(rr)[1]
        xmax <- terra::ext(rr)[2]
        ymin <- terra::ext(rr)[3]
        ymax <- terra::ext(rr)[4]

        lat_long <- check_lat_long(xmin, xmax, ymin, ymax)
    }

    # adjust units
    if (lat_long == TRUE) {
        h_units = "degrees"
        v_units = "meters"
    }

    # get labels & colors
    if (grepl("expos", filename)) {
        exp_all_labs <- c("", "protected", "exposed")
        exp_all_cols <- c("white", "grey", "green")

        exp_labs <- c(exp_all_labs[rr_min+1])
        exp_cols <- c(exp_all_cols[rr_min+1])

        if (rr_max > rr_min) {
            for (i in (rr_min+1):(rr_max)) {
                if (check_raster_value(rr, 1, i)) {
                    exp_labs <- append(exp_labs, exp_all_labs[i+1])
                    exp_cols <- append(exp_cols, exp_all_cols[i+1])
                }
            }
        }

    } else if (grepl("damage", filename)) {
        dam_all_labs <- c("", "None", "EF0", "EF1", "EF2", "EF3", "EF4", "EF5")
        dam_all_cols <- c("white", "grey", "purple", "blue", "green", "yellow", "orange", "red")

        dam_labs <- c(dam_all_labs[rr_min+1])
        dam_cols <- c(dam_all_cols[rr_min+1])

        if (rr_max > rr_min) {
            for (i in (rr_min+1):(rr_max)) {
                if (check_raster_value(rr, 1, i)) {
                    dam_labs <- append(dam_labs, dam_all_labs[i+1])
                    dam_cols <- append(dam_cols, dam_all_cols[i+1])
                }
            }
        }
    }

    # default palettes
    if (length(colormap) == 1) {
        if (colormap == "default") {
            # exposure map
            if (grepl("expos", filename)) {
                cmap <- exp_cols

            # damage map
            } else if (grepl("damage", filename)) {
                cmap <- dam_cols

            # dem, etc
            } else {
                cmap <- rev(grDevices::terrain.colors(255))
            }
        
        # raster palette
        } else if (colormap == "r_default") {
            cmap <- rev(grDevices::terrain.colors(255))
        }

    # user-specified palette
    } else {
        cmap <- colormap
    }

    # create plot
    oldpar <- graphics::par(no.readonly=TRUE)
    on.exit(graphics::par(oldpar))

    graphics::par(mar=c(5.1, 4.6, 4.1, 2.1))

    if (grepl("dem", filename)) {
        if (title == "") {
            title <- "Elevation"
        }
        v_units_str <- paste("  ", v_units, sep="")
        terra::plot(rr, main=title, xlab=h_units, ylab=h_units, type="continuous",
            plg=list(title=v_units_str), col=cmap)
        if (vector == TRUE) {
            terra::plot(boundaries, add=TRUE)
        }
    
    } else if (grepl("expos", filename)) {
        if (title == "") {
            x <- strsplit(filename, "-")[[1]]
            title <- paste("Exposure ", x[2], "-", x[3], sep="")
        }
         terra::plot(rr, main=title, xlab=h_units, ylab=h_units, type="classes",
            plg=list(title='  Exposure'), levels=exp_labs, all_levels=TRUE, col=cmap)
        if (vector == TRUE) {
            terra::plot(boundaries, add=TRUE)
        }

    } else if (grepl("damage", filename)) {
        if (title == "") {
            x <- strsplit(filename, "-")[[1]]
            title <- paste(x[1], "-", x[2], " Damage ", x[4], "-", x[5],sep="")
        }
        terra::plot(rr, main=title, xlab=h_units, ylab=h_units, type="classes",
            plg=list(title='  EF Scale'), levels=dam_labs, all_levels=TRUE, col=cmap)
        if (vector == TRUE) {
            terra::plot(boundaries, add=TRUE)
        }
    
    } else {
        if (title == "") {
            title <- filename
        }
        terra::plot(rr, main=title, xlab=h_units, ylab=h_units, col=cmap)
        if (vector == TRUE) {
            terra::plot(boundaries, add=TRUE)
        }
    }
}

