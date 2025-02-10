# TempEst2
TempEst 2/SCHEMA is "TEMPerature ESTimation, version 2, using Seasonal Conditions Historical Estimation with Modeled daily Anomaly". It estimates stream water temperature at a point, for streams of any size, using a data-driven model based on satellite remote sensing data.

This repository contains the model source code (`model.R`), Google Earth Engine data retrieval script (`eeretrieval.py`) and example data collection points (`datapts.py`), and R Notebooks demonstrating model usage (`demo.Rmd`) and a large set of model performance tests (`validation.Rmd`), used to make the performance claims reproducible.  In addition, the HydroShare data (see below) includes knitted output of both Notebooks (`demo.pdf`, `validation.pdf`) and the full training/testing dataset (`AllData.csv`, `ExtData.csv` including maximum temperatures).

## Quick Guide

Dependencies in R: `tidyverse`, `fields`. Additional dependencies for Google Earth Engine data retrieval: Python with the GEE API configured.

### Using a Trained Model

A pre-trained model is available in Releases as `model.rda`.  `load("model.rda")` will load a function called `model`.  The argument to `model` is a data frame with columns:

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

### Retrieving Prediction Data

Suitable prediction data can be retrieved in any way, but a script for retrieving the data from Google Earth Engine is included.  In `eeretrieval.py`, the function `getAllTimeseries` will retrieve data for specified points and dates to Google Drive.  Points are specified according to the example format in `datapts.py`: a list of `[[[longitude, latitude], "point id"], ...]`.

### Training a Model

Using all default options, `full.schema()` returns a function for training a model.  This function's argument is a data frame similar to above, but with one additional column: `temperature`, the observed temperature (daily mean, Celsius) for training.  `full.schema()(training data)` then returns a function like `model` above.  A typical approach to model training is for `id` to be set to training ga(u)ge IDs, such as the USGS stream temperature gage network, so that spatial observations can readily be linked to observed temperatures.  Prediction data and observations are retrieved independently of the training observations, but an example for USGS gages is included in `demo.Rmd`.

`full.schema` provides some additional options.

First, by default the model is trained with the two kriging models included, but the `full.schema` function itself is a generic framework for linking a seasonality (`sche`) and anomaly (`ma`) component, and these can be specified separately, e.g. to test a different approach without changing any of the infrastructure.

The third argument to `full.schema` is `rtn.model`, which is `FALSE` by default.  If it is set to `TRUE`, then, instead of returning a prediction function, `full.schema(rtn.model=TRUE)(training data)` will return a list of model components.  This is useful to examine model behaviors directly, such as extracting coefficients from the kriging models.  This is shown in `demo.Rmd`.

If the fourth argument to `full.schema`, `use.max`, is `TRUE` (default `FALSE`), then the model will be trained to predict both mean and maximum temperature, requiring a `temperature.max` column in the training data.

## Citation and More Information

If you use this model in your research, please cite Philippus, Corona, Schneider, Rust, and Hogue, "Satellite-Based Spatial-Statistical Modeling of Daily Stream Water Temperatures at the CONUS Scale" (in review).  This paper also contains a detailed description of model design and performance characteristics.

The seasonality component is based on the "three-sine" stream annual temperature cycle function described in: Philippus, Corona, and Hogue, (2024), "Improved annual temperature cycle function for stream seasonal thermal regimes", JAWRA, https://doi.org/10.1111/1752-1688.13228.

Full training datasets, pre-trained models, and a knitted model validation notebook PDF are available in CUAHSI HydroShare: https://doi.org/10.4211/hs.a8b243957f7946e388d10ab206990675.
