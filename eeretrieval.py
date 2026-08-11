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
# from datapts import datapts
from points_above import points_above_all, full_network

# Trigger the authentication flow.
# ee.Authenticate()

# Initialize the library.  It is now necessary to specify a Google Cloud project (this is still free for academic use).
ee.Initialize(project='earthenginestuff')

"""# Temperature Retrieval Setup

Retrieve TempEst inputs, adjusted to use MODIS LST rather than Ermida et al. model.

Required inputs: specific humidity (use NLDAS), trees, builtup, water, near-point LST, far away LST, latitude, longitude, and elevation.  Outputs also include ID, year, timestep, start, and end.  Note that if bulk running jobs is
not a logistical barrier, I can run one job per timestep.
"""

def getBuffered(dataset, band, pt, buffer, start=None, end=None):
  area = ee.Geometry.Point(pt).buffer(buffer)
  data = dataset.filter(ee.Filter.date(start, end)) if start is not None else dataset
  key = band + "_mean" if start is not None else "mean"
  result = data\
          .select(band)\
          .reduce(ee.Reducer.mean())\
          .reduceRegion(ee.Reducer.mean(), area, 10)
  result = ee.Algorithms.If(
      result.contains(key),
      result.getNumber(key),
      -999
      )
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
Drought = ee.ImageCollection("GRIDMET/DROUGHT");


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
      "humidity": ee.Number(getBuffered(humidity, "specific_humidity", pt, inner, start, end)),
      "drought": ee.Number(getBuffered(Drought, "spei1y", pt, inner, start, end))
      # "precip": ee.Number(getBuffered(precip, "total_precipitation",
      #                                   pt, inner, start, end))
  }).combine(getAbundances(pt, buffer=inner))

# droughtBands = ['pdsi', 'z', 'eddi30d', 'eddi90d', 'eddi180d', 'eddi1y',
#                 'eddi2y', 'eddi5y', 'spi30d', 'spi90d', 'spi180d', 'spi1y',
#                 'spi2y', 'spi5y', 'spei30d', 'spei90d', 'spei180d', 'spei1y',
#                 'spei2y', 'spei5y']
droughtBands = ['pdsi', 'z', 'eddi30d', 'eddi90d', 'eddi1y',
                'spei30d', 'spei90d', 'spei180d', 'spei1y']
droughtCols = ["id", "date"] + droughtBands

def getPtDrought(pt, id, year, time, start, end, inner=500, outer=1500):
  return ee.Dictionary({
      "id": id,
      "date": start
      } | {
          bandName: ee.Number(getBuffered(Drought, bandName, pt, inner, start, end))
          for bandName in droughtBands
          })


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

cols = ["id", "date", "lat", "lon", "elevation", "lst", "humidity", "drought", "shrubland",
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

def fullTimeseries(pts, times, basename, folder, retrieve=getPtData, cols=cols):
    # Retrieve full timeseries in a single pass per point. Note that a ~20-year
    # period of record may exceed the allowable batch size, so try subdividing
    # the times. 5-10 years seems to work.
    for pt in pts:
        desc = pt[1] + basename
        print(f"Running {desc}")
        result = mkFc(ee.List([
            retrieve(pt[0], pt[1], ts[0], ts[1], ts[2], ts[3])
            for ts in times
            ]))
        ee.batch.Export.table.toDrive(
            collection = result,
            description = desc,
            folder = folder,
            fileFormat = "CSV",
            selectors = cols
            ).start()



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
        break
    except Exception as err:
        print(err)
        print("Failed %s %s; retrying after delay" % (ts[0], ts[1]))
        sleep(wait * 10 if wait is not None else 180)
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

times_drt = [(str(x), str(y+1),
          str(dt.date(x, 1, 1) + dt.timedelta(y)),
          str(dt.date(x, 1, 1) + dt.timedelta(y+7)))
           for x in range(2001, 2027)
           for y in range(0, 365, 7)]
times = [(str(x), str(y+1),
          str(dt.date(x, 1, 1) + dt.timedelta(y)),
          str(dt.date(x, 1, 1) + dt.timedelta(y+7)))
           for x in range(2001, 2027)
           for y in range(0, 366)]

# 20 seconds works for 1300 points.  Scale accordingly.  The goal of the wait
# is to avoid exceeding 3000 jobs.
if __name__ == "__main__":
    maxN = 10
    # for i in range(maxN):
    # sites = ["USGS-09180000", "USGS-09315000", "USGS-10129900", "USGS-08279500", "USGS-09504950",
    #          "USGS-10261500", "USGS-14103000", "USGS-06040050", "USGS-06285100"]
    # id_bases = ["DoloresRiver", "GreenRiver", "SilverCreek", "RioGrande",
    #             "VerdeRiver", "MojaveRiver", "DeschutesRiver", "MadisonRiver",
    #             "ShoshoneRiver"]
    # datapts = points_above_all(sites, "usgs", 100, 1000, id_bases)
    # The UCRB has a total of 2,673,605 sites that come up on NLDI. We can handle about
    # 10,000 of those, or 1/250. This can be 1/10 points in 1/25 reaches, though in practice
    # 1/20 reaches suffices.
    co_pour = "USGS-09380000"
    datapts = full_network(co_pour, "usgs", "UCRB", 1000, reach_fraction = 1/25,
                           site_fraction = 1/10, False)
    for i in range(maxN):
        getAllTimeseries(datapts[i::maxN],
            times, f"TE{i}_", "DroughtAll", prt=True, wait=30)
    # getAllTimeseries(datapts,
    #     times_drt, "Drought", "DroughtRivers", prt=True, wait=30, retrieve=getPtDrought,
    #     cols=droughtCols)
    # runEcoregions(datapts, "", "Ecoregions")
    # Point-focused version...
    # maxN = 3
    # step = len(times) // maxN
    # for i in range(maxN):
    #     ts = times[(step*i):(step*(i+1))] if i < (maxN - 1) else times[(step*i):]
    #     fullTimeseries(datapts, ts, f"_Part{i}", "Montana")
