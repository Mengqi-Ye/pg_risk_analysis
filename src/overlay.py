import geopandas as gpd
import pandas as pd
from osgeo import ogr,gdal
import os
import xarray as xr
import rasterio
import numpy as np
import pyproj
from pygeos import from_wkb,from_wkt
import pygeos
from tqdm import tqdm
from shapely.wkb import loads
from pathlib import Path
import glob
from shapely.geometry import mapping
pd.options.mode.chained_assignment = None
from rasterio.mask import mask


def overlay_hazard_assets(df_ds,assets):
    """[summary]

    Args:
        df_ds ([type]): [description]
        assets ([type]): [description]

    Returns:
        [type]: [description]
    """
    #overlay 
    hazard_tree = pygeos.STRtree(df_ds.geometry.values)
    if (pygeos.get_type_id(assets.iloc[0].geometry) == 3) | (pygeos.get_type_id(assets.iloc[0].geometry) == 6):
        return  hazard_tree.query_bulk(assets.geometry,predicate='intersects')    
    else:
        return  hazard_tree.query_bulk(assets.buffered,predicate='intersects')
    