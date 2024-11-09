library(ExposR)
library(testthat)

local_edition(3)

# get expected values
exp_path <- system.file("", package="ExposR", mustWork=TRUE)
exposure_expected <- paste0(exp_path, '/exposure/expos-090-06.tif')
damage_expected <- paste0(exp_path, '/damage/AL1938-06-damage-06-2.tif')

# copy expected values to R temporary directory
tdir <- tempdir()
dir.create(paste0(tdir, '/dem'))
dir.create(paste0(tdir, '/exposure'))
dir.create(paste0(tdir, '/damage'))

file.copy(paste0(exp_path, '/dem/dem.tif'), paste0(tdir, '/dem/dem.tif'))

file.copy(paste0(exp_path, '/exposure/expos-000-06.tif'), paste0(tdir, '/exposure/expos-000-06.tif'))
file.copy(paste0(exp_path, '/exposure/expos-045-06.tif'), paste0(tdir, '/exposure/expos-045-06.tif'))
file.copy(paste0(exp_path, '/exposure/expos-090-06.tif'), paste0(tdir, '/exposure/expos-090-06.tif'))
file.copy(paste0(exp_path, '/exposure/expos-135-06.tif'), paste0(tdir, '/exposure/expos-135-06.tif'))
file.copy(paste0(exp_path, '/exposure/expos-180-06.tif'), paste0(tdir, '/exposure/expos-180-06.tif'))
file.copy(paste0(exp_path, '/exposure/expos-225-06.tif'), paste0(tdir, '/exposure/expos-225-06.tif'))
file.copy(paste0(exp_path, '/exposure/expos-270-06.tif'), paste0(tdir, '/exposure/expos-270-06.tif'))
file.copy(paste0(exp_path, '/exposure/expos-315-06.tif'), paste0(tdir, '/exposure/expos-315-06.tif'))

file.copy(paste0(exp_path, '/damage/AL1938-06.tif'), paste0(tdir, '/damage/AL1938-06.tif'))
file.copy(paste0(exp_path, '/damage/reproject.csv'), paste0(tdir, '/damage/reproject.csv'))

# get new values
expos_model(wind_direction=90, inflection_angle=6, exp_path=tdir)
exposure_new <- paste0(tdir, '/exposure/expos-090-06.tif')

expos_damage(hurricane="AL1938-06", inflection_angle=6, protect=2, exp_path=tdir)
damage_new <- paste0(tdir, '/damage/AL1938-06-damage-06-2.tif')

# test expos_summary
test_that("expos_summarize", {
	expect_snapshot_value(expos_summarize(filename="dem", exp_path=exp_path), style="serialize", cran=FALSE)
})

# test expos_model
test_that("expos_model", {
	expect_snapshot_file(exposure_new, exposure_expected, cran=FALSE)
})

# test expos_damage
test_that("expos_damage", {
	expect_snapshot_file(damage_new, damage_expected, cran=FALSE)
})

