import geopandas as gpd
import pandas as pd
from osgeo import ogr,gdal
import os,sys
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

# load from other py files within pg_risk_analysis
from extract import open_pg_data
from damage import assess_damage_pg

gdal.SetConfigOption("OSM_CONFIG_FILE", os.path.join('..',"osmconf.ini"))

# change paths to make it work on your own machine
data_path = os.path.join('C:\\','data','pg_risk_analysis')
tc_path = os.path.join(data_path,'tc_netcdf')
fl_path = os.path.join(data_path,'GLOFRIS')
osm_data_path = os.path.join('C:\\','data','country_osm')
pg_data_path = os.path.join(data_path,'pg_data')
vul_curve_path = os.path.join(data_path,'vulnerability_curves')


def country_analysis_pg(country_code,hazard_type): #
    """_summary_

    Args:
        country_code (_type_): _description_
        hazard_type (str, optional): _description_. Defaults to 'OSM'.

    Returns:
        _type_: _description_
    """
    # extract infrastructure data from OSM
    pg_infra = open_pg_data(country_code)

    # assess damage to wind storms
    pg_damage_infra = assess_damage_pg(country_code,pg_infra,hazard_type)

    return pg_damage_infra


if __name__ == "__main__":
    
    pg_damage_infra = country_analysis_pg(sys.argv[1],sys.argv[2]) #country_code, hazard_type