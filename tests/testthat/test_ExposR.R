library(ExposR)
library(testthat)

exp_path <- system.file("", package="ExposR", mustWork=TRUE)

# Test expos_summary
expect_snapshot_value(expos_summarize(filename="dem", exp_path=exp_path),
	style="serialize", cran=FALSE)

# Test expos_model
expect_snapshot_value(expos_model(wind_direction=90, inflection_angle=6, save=FALSE,
	exp_path=exp_path), style="serialize", cran=FALSE)

# Test expos_damage
expect_snapshot_value(expos_damage(hurricane="AL1938-06", inflection_angle=6, protect=2, save=FALSE,
	exp_path=exp_path), style="serialize", cran=FALSE)
