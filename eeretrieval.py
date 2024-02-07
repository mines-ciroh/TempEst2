"""
This script automatically runs Google Earth Engine jobs to retrieve specified
MODIS-based data.

Example usage is commented out at the end of the script and an example data 
points list is included (datapts.py).  Most likely, the relevant function to use
will be `getAllTimeseries`, which retrieves full data for a list of points and
dates, saving output CSVs to Google Drive.  After joining the CSVs together,
these can be used as direct inputs to a trained model in order to generate
predictions.
"""

import ee
from time import sleep
import datetime as dt

# Import a data points list with format:
# [[[lon, lat], gage ID]], e.g.
# datapts =  [[[-74.9658333,39.41916667], "01412000"]]
from datapts import datapts

# Trigger the authentication flow.
# ee.Authenticate()

# Initialize the library.
ee.Initialize()

"""# Temperature Retrieval Setup

Retrieve TempEst inputs, adjusted to use MODIS LST rather than Ermida et al. model.

Required inputs: specific humidity (use NLDAS), trees, builtup, water, near-point LST, far away LST, latitude, longitude, and elevation.  Outputs also include ID, year, timestep, start, and end.  Note that if bulk running jobs is
not a logistical barrier, I can run one job per timestep.
"""

def getBuffered(dataset, band, pt, buffer, start=None, end=None):
  area = ee.Geometry.Point(pt).buffer(buffer)
  data = dataset.filter(ee.Filter.date(start, end)) if start is not None else dataset
  result = data\
          .select(band)\
          .reduce(ee.Reducer.mean())\
          .reduceRegion(ee.Reducer.mean(), area, 10)\
          .getNumber(band + "_mean" if start is not None else "mean")
  return ee.Algorithms.If(
      result,
      result,
      ee.Number(-999)
  )


DEM = ee.Image("MERIT/DEM/v1_0_3")
humidity = ee.ImageCollection('NASA/NLDAS/FORA0125_H002')
landcover = ee.ImageCollection("ESA/WorldCover/v100").first()
LST = ee.ImageCollection("MODIS/061/MOD11A1")
precip = ee.ImageCollection('ECMWF/ERA5/DAILY').select('total_precipitation')
Ecoregions = ee.FeatureCollection('EPA/Ecoregions/2013/L3')


def getEcoregion(pt):
    # pt --> ee.Geometry.Point
    # returns ecoregion string
    filtered = Ecoregions.filterBounds(pt)
    return ee.Algorithms.If(
        ee.Number(filtered.size()).gt(0),
        filtered.first().getString("na_l1name"),
        ee.String("NA"))


def getEcoregions(points):
    # points --> list of [[lon, lat], id]
    # returns feature collection
    return ee.FeatureCollection([
        ee.Feature(None,
                   {
                       "id": pt[1],
                       "lat": pt[0][1],
                       "lon": pt[0][0],
                       "ecoregion": getEcoregion(ee.Geometry.Point(pt[0]))
                       }
                   )
        for pt in points
        ])

# Land cover
lctypes = {
  "10": "trees",
  "20": "shrubland",
  "30": "grassland",
  "40": "cropland",
  "50": "builtup",
  "60": "barren",
  "70": "snow",
  "80": "water",
  "90": "wetland",
  "95": "mangrove",
  "100": "moss"
}
lcn = ee.Dictionary(lctypes)
defaultlc = ee.Dictionary({x: 0 for x in lctypes.values()})

def getAbundances(point, datacol=landcover, buffer=1000):
  geom = ee.Geometry.Point(point).buffer(buffer)
  data = ee.Dictionary(datacol.reduceRegion(
    ee.Reducer.frequencyHistogram(),
    geom
    ).get("Map"))
  total = data.values().reduce(ee.Reducer.sum())
  return defaultlc.combine(data.rename(
      data.keys(),
      data.keys().map(lambda st: lcn.getString(st))
  )).map(lambda _, v: ee.Number(v).divide(total))

# Full data retrieval
# specific humidity (use NLDAS), trees, builtup, water, near-point LST, far away
# LST, latitude, longitude, and elevation
def getPtData(pt, id, year, time, start, end, inner=500, outer=1500):
  return ee.Dictionary({
      "id": id,
      "lat": pt[1],
      "lon": pt[0],
      "date": start,
      "elevation": ee.Number(getBuffered(DEM, "dem", pt, inner)),
      "lst": ee.Number(getBuffered(LST, "LST_Day_1km", pt, inner, start, end)).multiply(0.02).subtract(273.15),
      "humidity": ee.Number(getBuffered(humidity, "specific_humidity", pt, inner, start, end))
      # "precip": ee.Number(getBuffered(precip, "total_precipitation",
      #                                   pt, inner, start, end))
  }).combine(getAbundances(pt, buffer=inner))


def getPrecip(pt, id, year, time, start, end, inner=500, outer=1500):
    return ee.Dictionary({
        "id": id,
        "date": start,
        "precip": ee.Number(getBuffered(precip, "total_precipitation",
                                        pt, inner, start, end))
        })

def mkFeature(dict):
  return ee.Feature(None, dict)

def mkFc(output):
  return ee.FeatureCollection(output.map(mkFeature))

cols = ["id", "date", "lat", "lon", "elevation", "lst", "humidity", "shrubland",
        "grassland",
        "barren", "water"]

def getTimeseries(pt, id, times, basename, folder):
  # times: [("year", "time", "start", "end")]
  for ts in times:
    result = mkFc(ee.List([getPtData(pt, id, ts[0], ts[1], ts[2], ts[3])]))
    ee.batch.Export.table.toDrive(
        collection = result,
        description = basename + ts[0] + "_" + ts[1],
        folder = folder,
        fileFormat = "CSV",
        selectors = cols
    ).start()


def runEcoregions(points, folder, filename):
    result = getEcoregions(points)
    ee.batch.Export.table.toDrive(
        collection = result,
        description = filename,
        folder = folder,
        fileFormat = "CSV",
        selectors = ["id", "lat", "lon", "ecoregion"]
        ).start()

def getAllData(pts, year, time, start, end, inner=500, outer=1500,
               retrieve=getPtData):
  # pts -> [(lon, lat), "id"]
  return ee.List([
      retrieve(pt[0], pt[1], year, time, start, end, inner, outer)
      for pt in pts
  ])


def getAllTimeseries(pts, times, basename, folder, prt=False, wait=None,
                     retrieve=getPtData, cols=cols):
    # pts: [[[lon, lat], "point id"]]
    # times: [("year", "time", "start", "end")]
    # basename: Baseline file name that will be incremented for CSVs
    # folder: Google Drive folder to put results; I recommend making this folder
    #   manually first, because otherwise it will make several redundant folders.
    # prt: print progress
    # wait: how long to wait in between running jobs.  The goal of this is to
    #   avoid exceeding 3000 jobs, which is the maximum allowed by GEE.
    #   For 1300 points, a wait of 60 seconds works well.
    # retrieve: Retrieval function for a point.  See default.
    # cols: Which columns to store in the CSV.  See default.
  for ts in times:
    if prt:
        print("Running %s %s" % (ts[0], ts[1]))
    try:
        result = mkFc(getAllData(pts, ts[0], ts[1], ts[2], ts[3], retrieve=retrieve))
        ee.batch.Export.table.toDrive(
            collection = result,
            description = basename + ts[0] + "_" + ts[1],
            folder = folder,
            fileFormat = "CSV",
            selectors = cols
        ).start()
        if wait is not None:
            sleep(wait)
    except KeyboardInterrupt:
        print("Manually interrupted")
        raise KeyboardInterrupt
        break
    except Exception as err:
        print(err)
        print("Failed %s %s; retrying after delay" % (ts[0], ts[1]))
        sleep(wait * 10)
        try:
            result = mkFc(getAllData(pts, ts[0], ts[1], ts[2], ts[3], retrieve=retrieve))
            ee.batch.Export.table.toDrive(
                collection = result,
                description = basename + ts[0] + "_" + ts[1],
                folder = folder,
                fileFormat = "CSV",
                selectors = cols
            ).start()
            if wait is not None:
                sleep(wait)
        except Exception as err:
            print(err)
            print("Failed %s %s; not retrying" % (ts[0], ts[1]))

# Note that max jobs are 3000.

times = [(str(x), str(y+1),
          str(dt.date(x, 1, 1) + dt.timedelta(y)),
          str(dt.date(x, 1, 1) + dt.timedelta(y+1)))
           for x in range(2001, 2024)
           for y in range(365)]

# 20 seconds works for 1300 points.  Scale accordingly.  The goal of the wait
# is to avoid exceeding 3000 jobs.
if __name__ == "__main__":
    # getAllTimeseries(datapts,
    #     times, "AllData", "AFolder", prt=True, wait=10)
    runEcoregions(datapts, "", "Ecoregions")
