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

# The EXPOS model uses a digital elevation model to estimate exposed and
# protected areas for a given wind direction and inflection angle.

# The input file is assumed to be a raster file in GeoTiff format with 
# missing values represented by zero.  Cells may be rectangular but 
# horizontal and vertical units must be the same. Columns are assumed
# to be closely aligned with true north (if not, wind direction values
# must be adjusted accordingly). The name of the input file is 
# assumed to be "dem.tif".

# The output file is a raster file in GeoTiff format with the following
# values: 0 = missing data, 1 = protected, 2 = exposed. Output files
# are named "expos-xxx-yy.tif" where xxx is the wind direction and yy
# is the inflection angle.

# Emery R. Boose
# November 2021

# R version 4.1.1

# Required packages:
#  raster

### INTERNAL FUNCTIONS ####################################


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
#' @param save whether to save results to file
#' @param console whether to display messages in console
#' @return raster of modeled exposure values
#' @noRd

west_north_west <- function(wind_direction, inflection_angle, t_dir, save, console) {
    # get current working directory
    cwd <- getwd()
 
    # read dem file in GeoTiff format
    dem_file <- paste(cwd, "/dem.tif", sep="")
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

    # output
    if (save == TRUE) {
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

        # save modeled values in a Geotiff file
        expos_file = paste(cwd, "/expos-", formatC(wind_direction, width=3, flag="0"), "-", 
            formatC(inflection_angle, width=2, flag="0"), ".tif", sep="")

        raster::writeRaster(expos_r, expos_file, overwrite=TRUE)
    
        if (console == TRUE) {
            cat("\nSaving to", expos_file, "\n")
        }
    }
  
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
#' @param save whether to save results to file
#' @param console whether to display messages in console
#' @return raster of modeled exposure values
#' @noRd

north_north_west <- function(wind_direction, inflection_angle, t_dir, save, console) {
    # get current working directory
    cwd <- getwd()
 
    # read dem file in GeoTiff format
    dem_file <- paste(cwd, "/dem.tif", sep="")
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

    # output
    if (save == TRUE) {
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

        # save modeled values in a Geotiff file
        expos_file = paste(cwd, "/expos-", formatC(wind_direction, width=3, flag="0"), "-", 
            formatC(inflection_angle, width=2, flag="0"), ".tif", sep="")

        raster::writeRaster(expos_r, expos_file, overwrite=TRUE)
    
        if (console == TRUE) {
            cat("\nSaving to", expos_file, "\n")
        }
    }
  
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

    setwd(exp_path)

    if (console == TRUE) {
        cat("Path set to", exp_path)
    }
}


### MODELING FUNCTIONS ####################################

#' @title
#' Modeling Functions
#' @description
#' expos_model uses a raster file of elevation values, a specified wind
#' direction, and a specified inflection angle to create a raster file
#' of wind exposure values (0 = missing data, 1 = protected, 2 = exposed).
#' @param wind_direction wind direction (degrees)
#' @param inflection_angle inflection angle (degrees)
#' @param save whether to save results to file
#' @param console whether to display messages in console
#' @return no return value
#' @export
#' @rdname modeling

expos_model <- function(wind_direction, inflection_angle, save=TRUE, console=TRUE) {
    # get current working directory
    cwd <- getwd()
 
    # check wind direction
    if (wind_direction < 0 || wind_direction > 360) {
        stop("Please supply wind direction in range 0-360 degrees")
    }

    # check inflection angle
    if (inflection_angle < 0 || inflection_angle > 90) {
        stop("Please supply inflection angle in range 0-90 degrees")
    }

    # read dem file in GeoTiff format
    dem_path <- paste(cwd, "/dem.tif", sep="")
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

    # get angle of cell diagonal
    cell_diagonal <- 360 - 180*atan(cell_x/cell_y)/pi;

    # get transposed wind direction
    t_dir <- get_transposed_wind_direction(wind_direction)
  
    # create exposure map
    if (t_dir < cell_diagonal) {
        west_north_west(wind_direction, inflection_angle, t_dir, save, console)

    } else {
        north_north_west(wind_direction, inflection_angle, t_dir, save, console)
    }
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
    # get current working directory
    cwd <- getwd()
 
    # read file in GeoTiff format
    file_path <- paste(cwd, "/", filename, ".tif", sep="")
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
    st <- paste(st, "Values: ", round(val_min), " to ", round(val_max), "\n", sep="")
    
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
#' expos_plot creates a plot of a specified raster file.
#' @param filename name of input raster file
#' @return no return value
#' @export
#' @rdname plotting

expos_plot <- function(filename) {
    # get current working directory
    cwd <- getwd()
 
    # read file in GeoTiff format
    file_path <- paste(cwd, "/", filename, ".tif", sep="")
    check_file_exists(file_path)
    rr <- raster::raster(file_path)
 
    # create plot
    par(mar=c(5.1, 4.6, 4.1, 2.1))

    raster::plot(rr, main=filename)
}

