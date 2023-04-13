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
from utils import reproject,set_paths
from extract import extract_osm_infrastructure,open_pg_data
from hazard import clip_flood_data
from damage import assess_damage_osm,assess_damage_pg

gdal.SetConfigOption("OSM_CONFIG_FILE", os.path.join('..',"osmconf.ini"))

# change paths to make it work on your own machine

def country_analysis_osm(country_code,hazard_type):
    """_summary_

    Args:
        country_code (_type_): _description_
        hazard_type (str, optional): _description_. Defaults to 'OSM'.

    Returns:
        _type_: _description_
    """    
    # set paths
    data_path,tc_path,fl_path,osm_data_path,pg_data_path,vul_curve_path,output_path = set_paths()


    # extract infrastructure data from OSM
    osm_power_infra = extract_osm_infrastructure(country_code,osm_data_path)

    if hazard_type == 'fl':
        # extract flood hazard data
        clip_flood_data(country_code)

    # assess damage
    osm_damage_infra = assess_damage_osm(country_code,osm_power_infra,hazard_type)
    
    if hazard_type=='tc':
        climate_models = ['','_CMCC-CM2-VHR4','_CNRM-CM6-1-HR','_EC-Earth3P-HR','_HadGEM3-GC31-HM']
    elif hazard_type=='fl':
        climate_models = ['historical','rcp8p5']
        
    for i in range(len(osm_damage_infra)):
        for climate_model in climate_models:
            with pd.ExcelWriter(os.path.join(output_path,'{}_{}_damage_{}'.format(country_code,climate_model,i)+'.xlsx')) as writer:
                osm_damage_infra[i][climate_model].drop(['asset'], axis=1).to_excel(writer)
    
    return osm_damage_infra


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