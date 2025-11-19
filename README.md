# TempEst2
TempEst 2/SCHEMA is "TEMPerature ESTimation, version 2, using Seasonal Conditions Historical Estimation with Modeled daily Anomaly". It estimates stream water temperature at a point, for streams of any size, using a data-driven model based on satellite remote sensing data.

The [GitHub repository](https://github.com/mines-ciroh/tempest2) contains the model source code (`model.R`), Google Earth Engine data retrieval script (`eeretrieval.py`) and example data collection points (`datapts.py`), and R Notebooks demonstrating model usage (`demo.Rmd`) and a large set of model performance tests (`validation.Rmd`), used to make the performance claims reproducible.  In addition, the HydroShare data (see below) includes knitted output of both Notebooks (`demo.pdf`, `validation.pdf`) and the full training/testing dataset (`AllData.csv`, `ExtData.csv` including maximum temperatures).

In addition to the GEE data retrieval script, there is a utility script in Python, `points_above.py`, for generating lists of coordinates in the correct format, evenly-spaced above given target site(s) - this is useful for quickly generating study points over rivers of interest.

## Quick Guide

Dependencies in R: `tidyverse`, `fields`. The validation notebook also requires `terra`, `sf, `ggspatial`, and `hydroGOF`. Additional dependencies for Google Earth Engine data retrieval: Python with the GEE API configured.

Note that it is advisable to run at least 5 years and preferably 10, even if you are interested in a shorter period, so that TempEst 2 has adequate data to characterize the point of interest (see period-of-record tests in the paper).

### Using a Trained Model

A pre-trained model is available in the GitHub [Releases](https://github.com/mines-ciroh/TempEst2/releases) as `model.rda`.  `load("model.rda")` will load a function called `model`.  The argument to `model` is a data frame with columns:

* id (site ID; character)
* lon (longitude; decimal degrees East)
* lat (latitude; decimal degrees North)
* day (Julian day; integer)
* date (as an R Date)
* elevation (mean in a 500-m radius about point of interest; meters)
* water, shrubland, grassland, barren (land cover, abundance in a 500-m radius about point of interest; unitless, m2/m2)
* lst (land surface temperature, daily daytime mean in a 500-m radius about point of interest; Celsius)
* humidity (specific humidity, daily mean in a 500-m radius about point of interest; unitless, kg/kg)

Using this data frame as the argument to `model`, the function will return the same data frame with three or five new columns:

* temp.mod (modeled daily mean temperature; Celsius)
* temp.doy (modeled day-of-year mean temperature; Celsius)
* temp.anom (modeled temperature anomaly relative to day-of-year mean; Celsius)
* temp.plus (modeled maximum temperature relative to day-of-year mean; Celsius)
* temp.max (modeled maximum daily temperature; Celsius)

`temp.mod` = `temp.doy` + `temp.anom`.

These three (five) columns are the final output, and can be used to assess actual estimated temperature (`temp.mod`); typical seasonal conditions (`temp.doy`), such as to assess general long-term behavior; and departure from typical seasonal conditions (`temp.anom`), such as to assess heat waves or other extremes.

#### Full Steps (Pre-Trained Model)

1. Install `tidyverse` and `fields`.
2. Import the model code (`source("model.R")`).
3. Prepare prediction data (see below for default setup). If you prepare it yourself, the data frame must have the columns `id`, `lon`, `lat`, `day`, `date`, `elevation`, `water`, `shrubland`, `grassland`, `barren`, `lst`, `humidity` (see above for units and meaning).
4. `load()` the file containing the pre-trained model, e.g., `load("model.rda")`, or otherwise have a model function in memory.
5. Call the model function with the input dataframe as the argument, e.g. `result <- model(data)`. The output dataframe will add the above three or five columns to the input dataframe. Done.

Complete working example, assuming Earth Engine data outputs in `GEEData/` directory (see `demo.Rmd`):

```
source("model.R")
load("model.rda")
data <- map_dfr(
  list.files("GEEData/", full.names = T),
  \(fn) read_csv(fn, col_types="cDddddddddd") %>%
    mutate(day = as.integer(format(date, "%j"))) %>%
    filter(lst > -273)  # filters errors
)
pred <- model(data)
```

### Retrieving Prediction Data

Suitable prediction data can be retrieved in any way, but a script for retrieving the data from Google Earth Engine is included.  In `eeretrieval.py`, the function `getAllTimeseries` will retrieve data for specified points and dates to Google Drive.  Points are specified according to the example format in `datapts.py`: a list of `[[[longitude, latitude], "point id"], ...]`.

To run `eeretrieval.py`:

1. Specify your GEE cloud project name in the line `ee.Initialize(project='YOUR_CLOUD_PROJECT_HERE')`.
2. Provide a list of points in the format above, e.g. `datapts = [[[-87.3989084,34.28537757], '02450250'], [[-87.1702791,33.75760677], '02453500'], [[-86.8911032,33.5239948], '02458502']]` for the first few USGS gages used for development. Site IDs are used only to identify locations and have no internal meaning. If you are interested in multiple points along certain rivers, you may find the script below for generating evenly-spaced points to be useful.
3. Specify your dates of interest in the format `[("year", "time", "start", "end"), ...]`. The script contains an example for automatically doing so, which you can likely just modify to specify your years of interest.
4. Run `getAllTimeseries(points, times, basename, folder, prt=True, wait=...)`.  Basename is the base file name and folder is the Google Drive output directory. `wait` specifies pause in seconds between jobs to avoid exceeding the maximum of 3,000 simultaneous jobs; each daily job usually takes about 1 minute per 1,000 sites all told, so for a 20-year run (~7,000 jobs), a wait time of about 20 seconds per 1,000 sites will work. If you are running fewer than 3,000 jobs (about 8 years), you do not need a wait time.
5. After the jobs run, download the output data and use for prediction as above.

Complete working example (see `eeretrieval.py`):

```
import ee
from eeretrieval import getAllTimeseries

ee.Initialize(project='YOUR_CLOUD_PROJECT_HERE')

datapts = [[[-87.3989084,34.28537757], '02450250'], [[-87.1702791,33.75760677], '02453500'], [[-86.8911032,33.5239948], '02458502']]


times = [(str(x), str(y+1),
  str(dt.date(x, 1, 1) + dt.timedelta(y)),
  str(dt.date(x, 1, 1) + dt.timedelta(y+1)))
  for x in range(2001, 2024)
  for y in range(365)]
getAllTimeseries(datapts, times, "AllData", "AFolder", prt=True, wait=20)
```

### Generating Evenly-Spaced Points

You'll need a list of points of interest (USGS gage IDs, COMIDs, or coordinates - all in the same format), your desired maximum number of points per river, and a corresponding list of base IDs (and can optionally specify the resolution; the default is 1,000).

ID format:

- USGS gage: `"USGS-<ID>"` with `site_type="usgs"`
- COMID: just the COMID string with `site_type="comid"`
- Coordinates: (lon, lat) in decimal degrees with `site_type="coordinates"`

Then, just run `points_above.points_above_all([<site ids>], <site type>, <max points per river>, <resolution>, [<base ids>])` to get a list in the suitable format.

Example for Clear Creek and Sagehen Creek:

`points_above_all(["USGS-10343500", "USGS-06719505"], "usgs", 50, 1000, ["Sagehen", "ClearCreek"])`


### Training a Model

Using all default options, `full.schema()` returns a function for training a model.  This function's argument is a data frame similar to above, but with one additional column: `temperature`, the observed temperature (daily mean, Celsius) for training (also `temperature.max` if you wish to train a model for daily maximum temperature as well).  `full.schema()(training data)` then returns a function like `model` above.  A typical approach to model training is for `id` to be set to training ga(u)ge IDs, such as the USGS stream temperature gage network, so that spatial observations can readily be linked to observed temperatures.  Prediction data and observations are retrieved independently of the training observations, but an example for USGS gages is included in `demo.Rmd`.

`full.schema` provides some additional options.

First, by default the model is trained with the two kriging models included, but the `full.schema` function itself is a generic framework for linking a seasonality (`sche`) and anomaly (`ma`) component, and these can be specified separately, e.g. to test a different approach without changing any of the infrastructure.

The third argument to `full.schema` is `rtn.model`, which is `FALSE` by default.  If it is set to `TRUE`, then, instead of returning a prediction function, `full.schema(rtn.model=TRUE)(training data)` will return a list of model components.  This is useful to examine model behaviors directly, such as extracting coefficients from the kriging models.  This is shown in `demo.Rmd`.

If the fourth argument to `full.schema`, `use.max`, is `TRUE` (default `FALSE`), then the model will be trained to predict both mean and maximum temperature, requiring a `temperature.max` column in the training data.

## Citation and More Information

If you use this model in your research, please cite Philippus, Corona, Schneider, Rust, and Hogue, (2025), "Satellite-Based Spatial-Statistical Modeling of Daily Stream Water Temperatures at the CONUS Scale", Journal of Hydrology, https://doi.org/10.1016/j.jhydrol.2025.133321.  This paper also contains a detailed description of model design and performance characteristics.

The seasonality component is based on the "three-sine" stream annual temperature cycle function described in: Philippus, Corona, and Hogue, (2024), "Improved annual temperature cycle function for stream seasonal thermal regimes", JAWRA, https://doi.org/10.1111/1752-1688.13228.

Full training datasets, pre-trained models, and a knitted model validation notebook PDF are available in CUAHSI HydroShare: https://doi.org/10.4211/hs.a8b243957f7946e388d10ab206990675.

Source code is on GitHub: https://github.com/mines-ciroh/tempest2
