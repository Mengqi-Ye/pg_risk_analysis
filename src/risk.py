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
from extract import extract_osm_infrastructure
from damage import assess_damage_osm

gdal.SetConfigOption("OSM_CONFIG_FILE", os.path.join('..',"osmconf.ini"))

# change paths to make it work on your own machine
data_path = os.path.join('C:\\','data','pg_risk_analysis')
tc_path = os.path.join(data_path,'tc_netcdf')
fl_path = os.path.join(data_path,'GLOFRIS')
osm_data_path = os.path.join('C:\\','data','country_osm')
pg_data_path = os.path.join(data_path,'pg_data')

def country_analysis_osm(country_code,exposure_data='OSM'): #
    """_summary_

    Args:
        country_code (_type_): _description_
        exposure_data (str, optional): _description_. Defaults to 'OSM'.

    Returns:
        _type_: _description_
    """    
    if exposure_data == 'OSM':
        # extract infrastructure data from OSM
        ctry_power_infra = extract_osm_infrastructure(country_code,osm_data_path)
    
        # assess damage to wind storms
        #climate_models = ['_CMCC-CM2-VHR4','_CNRM-CM6-1-HR','_EC-Earth3P-HR','_HadGEM3-GC31-HM']
        ctry_damage_infra = assess_damage_osm(country_code,ctry_power_infra)
    
        return ctry_damage_infra
    """
    elif exposure_data == 'PG':
        # extract power grid data
        ctry_power_infra = extract_pg_data(country_code,pg_type)

        # extract wind data
        df_ds = extract_wind_data()
    
        # assess damage to wind storms
        climate_models = ['_CMCC-CM2-VHR4','_CNRM-CM6-1-HR','_EC-Earth3P-HR','_HadGEM3-GC31-HM']
        ctry_damage_infra = assess_damage_infrastructure(country_code,ctry_power_infra,climate_models)
    
        return ctry_damage_infra
        """

if __name__ == "__main__":

    ctry_damage_infra = country_analysis_osm(sys.argv[1]) 