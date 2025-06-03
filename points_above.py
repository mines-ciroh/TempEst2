# Written by Daniel Philippus beginning June 2, 2025.
# This is a quick utility script for identifying 1-km-spaced points on the mainstem of
# a given stream starting from a given location, then printing in TE2 format.

# Input: USGS gage IDs, COMIDs, or coordinate pairs
# Output: [[[lon, lat], "id"]]

import numpy as np
import shapely as shp
import geopandas as gpd
from pynhd import NLDI
nldi = NLDI()

# This is the HRRR projection for convenience.
projstr = "+proj=lcc +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"

# Sagehen Cr.
test_id = "USGS-10343500"
test_co = (-120.2379547, 39.43154577)


def points_above_site(site, site_type, id_base, N, resolution=1000):
    """
    Retrieve *approximately* evenly-spaced points upstream of the given site.

    Parameters
    ----------
    site : string or tuple of floats
        Site location. Either a USGS/COMID string or tuple of
        (longitude, latitude) identifying the downstream point. COMID
        strings are just the ID, not any prefix. USGS should be USGS-<id>.
    site_type : str
        The type of site identifier. "usgs", "comid" or "coordinates".
    N : int
        The number of points to select.
    resolution : int, optional
        Resolution in meters. The default is 1000.
    id_base : str
        Base name for point IDs.

    Returns
    -------
    list of [[longitude, latitude], point_id].

    """
    dist = N * resolution  # convert to km
    if site_type == "coordinates":
        nav = nldi.navigate_byloc(site, "upstreamMain", "flowlines", distance=dist)
    else:
        source_typ = "nwissite" if site_type == "usgs" else "comid"
        nav = nldi.navigate_byid(source_typ, site, "upstreamMain", "flowlines",
                                 distance=dist)
    # Nav is in order from bottom to top. However, each LineString within
    # nav goes from top to bottom.
    # Ideally, we'd like to use GeoPandas' interpolate, but that doesn't work
    # in long geometries - it just returns the endpoints. Instead, we have to
    # decide where to put things.
    nav = nav.to_crs(projstr)  # need to be in meters to be accurate
    lengths = nav.length.cumsum()  # cumulative lengths from the bottom
    steps = np.arange(0, sum(nav.length), resolution)
    result = []
    for step in steps:
        len_row = lengths[lengths > step].head(1)
        # How far it actually gets
        row_max = len_row.iloc[0]
        excess = row_max - step
        row_index = len_row.index
        geom = nav.loc[row_index].geometry.iloc[0]
        if excess > resolution:
            # Negative indices go from the end. Solves the direction problem.
            increments = -np.arange(1, excess, resolution)
            multi = True
        else:
            # Just go "the excess" in from the start (top)
            increments = excess
            multi = False
        points = shp.line_interpolate_point(geom, increments)
        if multi:
            result += points
        else:
            result.append(points)
    # Now we should have a list of points in order from bottom to top.
    res = gpd.GeoDataFrame({
        "site_id": [id_base + "_" + str(i) for i in range(len(result))]
        },
        geometry=gpd.GeoSeries(result)).to_crs(4326) # return to degrees
    return res


def points_above_all(sites, site_type, N, resolution=0.01, id_bases=None):
    """
    Retrieve *approximately* evenly-spaced points upstream of the given sites.

    Parameters
    ----------
    sites : list of strings or tuples
        List of sites. These are either USGS or COMID strings or tuples of
        (longitude, latitude) identifying the downstream point. COMID
        strings are just the ID, not any prefix. USGS should be USGS-<id>.
    site_type : str
        The type of site identifiers. "usgs", "comid" or "coordinates".
    N : int
        The number of points to select.
    resolution : float
        Resolution in decimal degrees. 0.01 is approximately 1 km.
    id_base : list of str, optional
        Base name for point IDs, same length as sites. The default is None.

    Returns
    -------
    List of [[longitude, latitude], point_id].
    
    Notes
    -----
    point_id will be id_base + "_" + index from the bottom. If id_base is
    unspecified, it is derived from the site_ids.

    """
    if id_bases is None:
        if site_type == "coordinates":
            id_bases = [str(lon) + "_" + str(lat)
                        for (lon, lat) in sites]
        else:
            id_bases = [site_type + "_" + site for site in sites]
    all_sites = [points_above_site(site,
                                   site_type,
                                   id_bases[i],
                                   N,
                                   resolution)
                 for (i, site) in enumerate(sites)]
    return [site_row for
            site_set in all_sites for
            site_row in site_set]
    
