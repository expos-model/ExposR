library(ExposR)
library(testthat)

# Test expos_summary
test.expected <- system.file("dem", "dem_summary.expected", package="ExposR", mustWork=TRUE)
exp_path <- system.file("", package="ExposR", mustWork=TRUE)
expect_snapshot_value(expos_summarize(filename="dem", exp_path=exp_path), test.expected, 
	style="serialize", cran=FALSE)

# Test expos_model
test.expected <- system.file("exposure", "expos-090-06.tif", package="ExposR", mustWork=TRUE)
exp_path <- system.file("", package="ExposR", mustWork=TRUE)
expect_snapshot_value(expos_model(wind_direction=90, inflection_angle=6, save=FALSE, 
	exp_path=exp_path), test.expected, style="serialize", cran=FALSE)

# Test expos_damage
test.expected <- system.file("damage", "AL1938-06-damage-06-2.tif", package="ExposR", mustWork=TRUE)
exp_path <- system.file("", package="ExposR", mustWork=TRUE)
expect_snapshot_value(expos_damage(hurricane="AL1938-06", inflection_angle=6, protect=2, save=FALSE,
	exp_path=exp_path), test.expected, style="serialize", cran=FALSE)
