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
#import rioxarray

# load from other py files within pg_risk_analysis
from utils import reproject


gdal.SetConfigOption("OSM_CONFIG_FILE", os.path.join('..',"osmconf.ini"))

# change paths to make it work on your own machine
data_path = os.path.join('C:\\','data','pg_risk_analysis')
tc_path = os.path.join(data_path,'tc_netcdf')
fl_path = os.path.join(data_path,'GLOFRIS')
osm_data_path = os.path.join('C:\\','data','country_osm')
pg_data_path = os.path.join(data_path,'pg_data')


##### ##### ##### ##### ##### ##### ##### #####  
##### ##### ##### STORM DATA  ##### ##### ##### 
##### ##### ##### ##### ##### ##### ##### #####  

def open_storm_data(climate_model):
    """_summary_

    Args:
        climate_model (_type_): _description_

    Returns:
        _type_: _description_
    """    
    
    with xr.open_dataset(os.path.join(tc_path,'STORM_FIXED_RETURN_PERIODS{}_WP.nc'.format(climate_model))) as ds:
        """
        TC climate model:
            CMCC-CM2-VHR4
            CNRM-CM6-1-HR
            EC-Earth3P-HR
            HadGEM3-GC31-HM
        """
        
        # get the mean values
        df_ds = ds['mean'].to_dataframe().unstack(level=2).reset_index()

        # create geometry values and drop lat lon columns
        df_ds['geometry'] = [pygeos.points(x) for x in list(zip(df_ds['lon'],df_ds['lat']))]
        df_ds = df_ds.drop(['lat','lon'],axis=1,level=0)
        #print(df_ds)
        
        #rename columns to return periods
        return_periods = ['1_{}{}'.format(int(x),climate_model) for x in ds['rp']]
        df_ds.columns = ['1_{}{}'.format(int(x),climate_model) for x in list(df_ds.columns.get_level_values(1))[:-1]]+['geometry']     
        df_ds['geometry'] = pygeos.buffer(df_ds.geometry,radius=0.1/2,cap_style='square').values
        df_ds['geometry'] = reproject(df_ds)
            
        # drop all non values to reduce size
        #if climate_model == '':
        #    df_ds = df_ds.loc[~df_ds['1_10000'].isna()].reset_index(drop=True)
        
        df_ds = df_ds.loc[~df_ds['1_10000{}'.format(climate_model)].isna()].reset_index(drop=True)
        df_ds = df_ds.fillna(0)
        df_ds = df_ds[['1_{}{}'.format(int(x),climate_model) for x in [10,50,100,500,1000]]+['geometry']]
        #print(df_ds)

    return df_ds

# load hazard data 
def extract_wind_data():
    """_summary_

    Returns:
        _type_: _description_
    """   

    climate_models = ['_CMCC-CM2-VHR4','_CNRM-CM6-1-HR','_EC-Earth3P-HR','_HadGEM3-GC31-HM']
    df_ds = {}
    for climate_model in climate_models:
        df_ds_cl = open_storm_data(climate_model)
        df_ds[climate_model] = df_ds_cl
    
    return df_ds


##### ##### ##### ##### ##### ##### ##### #####  
##### ##### ##### FLOOD DATA  ##### ##### ##### 
##### ##### ##### ##### ##### ##### ##### #####     

def clip_flood_data(country_code,time_period='HIST'):
    """_summary_

    Args:
        country_code (_type_): _description_
        time_period (str, optional): _description_. Defaults to 'HIST'.
    """    
    # load country geometry file and create geometry to clip
    ne_countries = gpd.read_file('C:\\Data\\natural_earth\\ne_10m_admin_0_countries.shp') 
    geometry = ne_countries.loc[ne_countries['ISO_A3']==country_code].geometry.values[0]
    geoms = [mapping(geometry)]
    
    rps = ['0010','0050','0100','0500','1000']
    for rp in rps: 
        input_file = os.path.join(fl_path,time_period,'global',
                                  'inuncoast_historical_nosub_hist_rp{}_0.tif'.format(rp))

        # load raster file and save clipped version
        with rasterio.open(input_file) as src:
            out_image, out_transform = mask(src, geoms, crop=True) #rasterio.mask.mask(src, geoms, crop=True)
            out_meta = src.meta

            out_meta.update({"driver": "GTiff",
                     "height": out_image.shape[1],
                     "width": out_image.shape[2],
                     "transform": out_transform})

            file_path = os.path.join(fl_path,time_period,'country',
                                            '_'.join([country_code]+input_file.split('_')[3:]))
            
            with rasterio.open(file_path, "w", **out_meta) as dest:
                dest.write(out_image) 

def open_flood_data(country_code,time_period='HIST'):
    """_summary_

    Args:
        country_code (_type_): _description_
        time_period (str, optional): _description_. Defaults to 'HIST'.

    Returns:
        _type_: _description_
    """    
    # THIS STILL NEEDS TO BE TESTED WITH GLOFRIS DATA
    with xr.open_dataset(os.path.join(fl_path,'HIST/inuncoast_historical_nosub_hist_rp0500_0.nc')) as ds: #, engine="rasterio"
        df_ds = ds.to_dataframe().reset_index()
        df_ds['geometry'] = pygeos.points(df_ds.x,y=df_ds.y)
        df_ds = df_ds.rename(columns={'band_data': 'hazard_intensity'})
        df_ds = df_ds.drop(['band','x', 'y','spatial_ref'], axis=1)
        df_ds = df_ds.dropna()
        df_ds = df_ds.reset_index(drop=True)
        df_ds.geometry= pygeos.buffer(df_ds.geometry,radius=20/2,cap_style='square').values
        df_ds['geometry'] = reproject(df_ds)

    return df_ds