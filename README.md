## EXPOS Model

The EXPOS model uses a digital elevation model (dem) to estimate exposed and
protected areas for a given wind direction and inflection angle. The resulting 
topograhic exposure maps can be combined with output from the Hurrecon model 
to estimate hurricane wind damage across a region

EXPOS contains two main functions:

1. The <i>expos_model</i> function estmates topopgrahic exposure for a specified
wind direction and inflection angle. The input file is assumed to be a raster of 
elevation values in GeoTiff format with missing values represented by zero. Cells 
may be rectangular. If a geographic coordinate system is used, horizontal and 
vertical units are assumed to be degrees and meters, respectively; otherwise 
horizontal and vertical units must be the same. Columns are assumed to be closely 
aligned with true North (0 degrees); if not, the map orientation (in degrees) must 
be specified. The name of the input file is assumed to be "dem.tif".

The output file is a raster file in GeoTiff format with the following values: 
0 = missing data, 1 = protected, 2 = exposed. Output files are named 
"expos-xxx-yy.tif" where xxx is the wind direction and yy is the inflection angle.

In previous studies, a spatial resolution of 30 meters and an inflection angle of 
6 degrees were found to work well (see below). Note that increasing the inflection 
angle tends to decrease the size and number of protected areas while increasing 
the odds that these sites are protected.

2. The <i>expos_damage</i> function estimates regional hurricane damage where
topographic exposure at each location is determined by peak wind direction. If a 
location is protected, the enhanced Fujita scale (EF) rating is reduced by a
specified amount. This function requires a hurricane GeoTiff file created by 
Hurrecon, exposure files created by Expos for the eight cardinal wind directions 
(N, NE, E, etc), and a reprojection file in csv format that contains lat/long 
coordinates for the lower left and upper right corners of the digital elevation 
model.

The output file is a raster file in GeoTiff format with the following values:
0 = missing, 1 = no damage, 2 = F0 damage, 3 = F1 damage, 4 = F2 damage, 
5 = F3 damage, 6 = F4 damage, 7 = F5 damage. Output files are named 
"hhhh-damage-yy-z.tif" where hhhh is the hurricane ID, yy is the inflection 
angle, and z is the reduction in EF rating for protected areas.

EXPOS is available in both R (ExposR) and Python (ExposPython) versions. 
The model is an updated version of the original EXPOS model written in Borland 
Pascal for use with Idrisi (see below for details).

Please note: both versions are currently under development and subject to change.

## Getting Started

Here are the basic steps for using the model. Please see below for more details.

1. Download the R or Python version of the model from GitHub.
2. Create a directory for a particular study area with subdirectories as described 
below.
3. Copy the digital elevation file to the dem subdirectory and rename it "dem.tif".
4. Download geographic and political boundary shapefiles for the desired study area.
Copy these files to the vector subdirectory and rename so the first name of each file 
is "boundaries".
5. Use the <i>expos_model</i> function to create exposure maps for different wind 
directions and inflection angles.
6. Use the <i>expos_damage</i> function to create maps of enhanced Fujita scale wind 
damage for particular hurricanes.

## Details

All user functions begin with "expos". The wind direction and inflection angle must be
specified in degrees.

The user specifies a directory (<i>exp_dir</i>) for a given study area. Input and output
files are stored on the following subdirectories of this directory:


```{r}
exp_dir/dem
exp_dir/exposure
exp_dir/damage
exp_dir/vector
```

The dem subdirectory contains the input elevation file. The exposure subdirectory 
contains output files from the <i>expos_model</i> function. The damage subdirectory 
contains input files from the Hurrecon model and output files from the 
<i>expos_damage</i> function. Shapefiles that contain geographic and political 
boundaries for viewing results are stored on the vector subdirectory.

To run the model, run <i>expos.R</i> (R) or <i>expos.py</i> (Python). 

The R version may also be installed as an R package using the devtools package:

```{r}	
devtools::install_github("expos-model/ExposR")
```

This has the advantage of providing readily accessible help messages for each 
function.

## Model Functions

```{r}	
expos_set_path
expos_get_path
expos_model
expos_damage
expos_summarize
expos_plot
```

The <i>expos_set_path</i> function sets the path for the current set of model runs. 
The <i>expos_get_path</i> function returns the current path. Use <i>expos_set_path</i> 
before using other functions.

The <i>expos_model</i> function creates a raster file of topographic wind exposure 
as a function of wind direction and inflection angle.

The <i>expos_damage</i> function uses output from Hurrecon and Expos to create a 
raster file of wind damage where topograhic exposure at each location is determined 
by peak wind direction. If a location is protected, the enhanced Fujita scale 
rating is reduced by a specfied amount.

The <i>expos_summarize</i> function displays summary information for a specified 
raster file, including the number of rows and columns, spatial extent, cell height 
and width, and minimum and maximum value.

The <i>expos_plot</i> function creates a plot of a specified raster file.

## Examples


```{r}
expos_set_path("c:/expos/r/mass_30m")
expos_get_path()
expos_model(wind_direction=90, inflection_angle=6)
expos_damage("AL1938-06", inflection_angle=6, protect=2)
expos_summarize("dem")

expos_plot("dem")
expos_plot("expos-090-06")
expos_plot("AL1938-06-damage-06-2")
```

## History

The original EXPOS model was written in Borland Pascal and depended on Idrisi 
for spatial visualization. The model was used in published studies of the ecological 
impacts of historical hurricanes in New England and Puerto Rico:

* Boose, E. R., Chamberlin, K. E., Foster, D. R. 2001. Landscape and regional impacts 
of hurricanes in New England. Ecological Monographs 71: 27-48.

* Boose, E. R., Serrano, M. I., Foster, D. R. 2004. Landscape and regional impacts of 
hurricanes in Puerto Rico. Ecological Monographs 74: 335-352.

