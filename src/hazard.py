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


# load from other py files within pg_risk_analysis
from utils import reproject,set_paths

gdal.SetConfigOption("OSM_CONFIG_FILE", os.path.join('..',"osmconf.ini"))

##### ##### ##### ##### ##### ##### ##### #####  
##### ##### ##### STORM DATA  ##### ##### ##### 
##### ##### ##### ##### ##### ##### ##### #####  

def open_storm_data(country_code):
    """
    This function loads STORM data for a given country code, clips it based on the country geometry,
    and combines data from different basins and climate models.

    Args:
    - country_code (str): a 3-letter ISO code of the country of interest

    Returns:
    - df_ds (dict): a dictionary containing STORM data for different climate models, organized by basin
    """
    # set paths
    data_path,tc_path,fl_path,osm_data_path,pg_data_path,vul_curve_path,output_path,ne_path = set_paths()

    # list of available climate models
    climate_models = ['','_CMCC-CM2-VHR4','_CNRM-CM6-1-HR','_EC-Earth3P-HR','_HadGEM3-GC31-HM']

    # dictionary of basins for each country
    country_basin = {
        "BRN": ["WP"],
        "KHM": ["WP"],
        "CHN": ["WP", "NI"],
        "IDN": ["SI", "SP", "NI", "WP"],
        "JPN": ["WP"],
        "LAO": ["WP"],
        "MYS": ["WP", "NI"],
        "MNG": ["WP", "NI"],
        "MMR": ["NI", "WP"],
        "PRK": ["WP"],
        "PHL": ["WP"],
        "SGP": ["WP"],
        "KOR": ["WP"],
        "TWN": ["WP"],
        "THA": ["WP", "NI"],
        "VNM": ["WP"]
    }

    # load country geometry file and create geometry to clip
    ne_countries = gpd.read_file(ne_path)
    bbox = ne_countries.loc[ne_countries['ISO_A3']==country_code].geometry.buffer(1).values[0].bounds

    df_ds = {}
    for climate_model in climate_models:
        concat_prep = []

        #combine STORM data from different basins
        if "WP" in country_basin[country_code]:
            WP = load_storm_data(climate_model,'WP',bbox)
            concat_prep.append(WP)
        if "SP" in country_basin[country_code]:
            SP = load_storm_data(climate_model,'SP',bbox)
            concat_prep.append(SP)
        if "NI" in country_basin[country_code]:            
            NI = load_storm_data(climate_model,'NI',bbox)
            concat_prep.append(NI)            
        if "SI" in country_basin[country_code]:       
            SI = load_storm_data(climate_model,'SI',bbox)
            concat_prep.append(SI)            
                   
        df_ds_cl = pd.concat(concat_prep, keys=country_basin[country_code])
        df_ds_cl = df_ds_cl.reset_index(drop=True)
        df_ds[climate_model] = df_ds_cl

    return df_ds

def load_storm_data(climate_model,basin,bbox):
    """
    Load storm data from a NetCDF file and process it to return a pandas DataFrame.

    Parameters:
    - climate_model (str): name of the climate model
    - basin (str): name of the basin
    - bbox (tuple): bounding box coordinates in the format (minx, miny, maxx, maxy)
    - ne_crs (str): CRS string of the North-East projection

    Returns:
    - df_ds (pd.DataFrame): pandas DataFrame with interpolated wind speeds for different return periods and geometry column
    """
    # set paths
    data_path,tc_path,fl_path,osm_data_path,pg_data_path,vul_curve_path,output_path,ne_path = set_paths()

    filename = os.path.join(tc_path, f'STORM_FIXED_RETURN_PERIODS{climate_model}_{basin}.nc')
    
    # load data from NetCDF file
    with xr.open_dataset(filename) as ds:
        
        # convert data to WGS84 CRS
        ds.rio.write_crs(4326, inplace=True)
        ds = ds.rio.clip_box(minx=bbox[0], miny=bbox[1], maxx=bbox[2], maxy=bbox[3])
        
        # get the mean values
        df_ds = ds['mean'].to_dataframe().unstack(level=2).reset_index()

        # create geometry values and drop lat lon columns
        df_ds['geometry'] = [pygeos.points(x) for x in list(zip(df_ds['lon'], df_ds['lat']))]
        df_ds = df_ds.drop(['lat', 'lon'], axis=1, level=0)
        
        # interpolate wind speeds of 1, 2, and 5-yr return period
        ## rename columns to return periods (must be integer for interpolating)
        df_ds_geometry = pd.DataFrame()
        df_ds_geometry['geometry'] = df_ds['geometry']
        df_ds = df_ds.drop(['geometry'], axis=1, level=0)
        df_ds = df_ds['mean']
        df_ds.columns = [int(x) for x in ds['mean']['rp']]
        df_ds[1] = np.nan
        df_ds[2] = np.nan
        df_ds[5] = np.nan
        df_ds[25] = np.nan
        df_ds[250] = np.nan
        df_ds = df_ds.reindex(sorted(df_ds.columns), axis=1)
        df_ds = df_ds.interpolate(method='linear', axis=1, limit_direction='both')
        df_ds['geometry'] = df_ds_geometry['geometry']
        df_ds = df_ds[[1, 2, 5, 10, 25, 50, 100, 250, 500, 1000, 'geometry']]
        
        # rename columns to return periods
        df_ds.columns = ['1_{}{}'.format(int(x), climate_model) for x in [1, 2, 5, 10, 25, 50, 100, 250, 500, 1000]] +['geometry']     
        df_ds['geometry'] = pygeos.buffer(df_ds.geometry, radius=0.1/2, cap_style='square').values
        
        # reproject the geometry column to the specified CRS
        df_ds['geometry'] = reproject(df_ds)
            
        # drop all non values to reduce size
        #df_ds = df_ds.loc[~df_ds['1_10000{}'.format(climate_model)].isna()].reset_index(drop=True)
        df_ds = df_ds.fillna(0)

    return df_ds

##### ##### ##### ##### ##### ##### ##### #####  
##### ##### ##### FLOOD DATA  ##### ##### ##### 
##### ##### ##### ##### ##### ##### ##### ##### 

def load_flood_data(country_code,climate_model):

    # set paths
    data_path,tc_path,fl_path,osm_data_path,pg_data_path,vul_curve_path,output_path,ne_path = set_paths()
<<<<<<< HEAD

    # load country geometry file and create geometry to clip
    ne_countries = gpd.read_file(ne_path)
    geometry = ne_countries.loc[ne_countries['ISO_A3']==country_code].geometry.values[0]
    geoms = [mapping(geometry)]
    
    #climate_model: historical, rcp4p5, rcp8p5; time_period: hist, 2030, 2050, 2080
    rps = ['0001','0002','0005','0010','0025','0050','0100','0250','0500','1000']
    climate_models = ['historical','rcp8p5']
    
    #"/scistor/ivm/data_catalogue/open_street_map/pg_risk_analysis/GLOFRIS/global/inuncoast_historical_nosub_hist_rp0001_0.tif"

    for rp in rps:
        #global input_file
        for climate_model in climate_models:
            if climate_model=='historical':
                input_file = os.path.join(fl_path,'global',
                                          'inuncoast_{}_nosub_hist_rp{}_0.tif'.format(climate_model,rp)) 
 
            elif climate_model=='rcp8p5':
                input_file = os.path.join(fl_path,'global',
                                          'inuncoast_{}_nosub_2030_rp{}_0.tif'.format(climate_model,rp))
            
            # load raster file and save clipped version
            with rasterio.open(input_file) as src:
                out_image, out_transform = mask(src, geoms, crop=True)
                out_meta = src.meta

                out_meta.update({"driver": "GTiff",
                         "height": out_image.shape[1],
                         "width": out_image.shape[2],
                         "transform": out_transform})

                if 'scistor' in fl_path:
                    file_path = os.path.join(fl_path,'country','_'.join([country_code]+input_file.split('_')[6:]))
                else:
                    file_path = os.path.join(fl_path,'country','_'.join([country_code]+input_file.split('_')[3:]))

                with rasterio.open(file_path, "w", **out_meta) as dest:
                    dest.write(out_image)

def load_flood_data(country_code,climate_model):

    # set paths
    data_path,tc_path,fl_path,osm_data_path,pg_data_path,vul_curve_path,output_path,ne_path = set_paths()
=======
>>>>>>> e8e251aee62cccc762e09b76df9b70d4b89adb3d
     
    rps = ['0001','0002','0005','0010','0025','0050','0100','0250','0500','1000']
    collect_df_ds = []
    
    if climate_model=='historical':
        print('Loading historical coastal flood data ...')
        for rp in rps:
            #for file in files:
            file_path = os.path.join(fl_path,'country','{}_{}_nosub_hist_rp{}_0.tif'.format(country_code,climate_model,rp))
            with xr.open_dataset(file_path) as ds: #, engine="rasterio"
                df_ds = ds.to_dataframe().reset_index()
                df_ds['geometry'] = pygeos.points(df_ds.x,y=df_ds.y)
                df_ds = df_ds.rename(columns={'band_data': 'rp'+rp}) #rename to return period
                
                # move from meters to centimeters
                df_ds['rp'+rp] = (df_ds['rp'+rp]*100)         
                df_ds = df_ds.drop(['band','x', 'y','spatial_ref'], axis=1)
                df_ds = df_ds.dropna()
                df_ds = df_ds.reset_index(drop=True)
                df_ds.geometry= pygeos.buffer(df_ds.geometry,radius=0.00833/2,cap_style='square').values  #?????????????????????????
                df_ds['geometry'] = reproject(df_ds)
                collect_df_ds.append(df_ds)

        df_all = collect_df_ds[0].merge(collect_df_ds[1]).merge(collect_df_ds[2]).merge(collect_df_ds[3]).merge(collect_df_ds[4])\
                 .merge(collect_df_ds[5]).merge(collect_df_ds[6]).merge(collect_df_ds[7]).merge(collect_df_ds[8]).merge(collect_df_ds[9])
        df_all = df_all.loc[df_all['rp1000']>0].reset_index(drop=True)

    elif climate_model=='rcp8p5':
        print('Loading future coastal flood data ...')
        for rp in rps:
            #for file in files:
            file_path = os.path.join(fl_path,'country','{}_{}_nosub_2030_rp{}_0.tif'.format(country_code,climate_model,rp))
            with xr.open_dataset(file_path) as ds: #, engine="rasterio"
                df_ds = ds.to_dataframe().reset_index()
                df_ds['geometry'] = pygeos.points(df_ds.x,y=df_ds.y)
                df_ds = df_ds.rename(columns={'band_data': 'rp'+rp}) #rename to return period
                df_ds['rp'+rp] = (df_ds['rp'+rp]*100)
                df_ds = df_ds.drop(['band','x', 'y','spatial_ref'], axis=1)
                df_ds = df_ds.dropna()
                df_ds = df_ds.reset_index(drop=True)
                df_ds.geometry= pygeos.buffer(df_ds.geometry,radius=0.00833/2,cap_style='square').values
                df_ds['geometry'] = reproject(df_ds)
                collect_df_ds.append(df_ds)

        df_all = collect_df_ds[0].merge(collect_df_ds[1]).merge(collect_df_ds[2]).merge(collect_df_ds[3]).merge(collect_df_ds[4])\
                 .merge(collect_df_ds[5]).merge(collect_df_ds[6]).merge(collect_df_ds[7]).merge(collect_df_ds[8]).merge(collect_df_ds[9])

        df_all = df_all.loc[df_all['rp1000']>0].reset_index(drop=True)
    return df_all

def open_flood_data(country_code):
    
    climate_models = ['historical','rcp8p5']
    df_ds = {}
    for climate_model in climate_models:
        df_ds_sc = load_flood_data(country_code,climate_model)
<<<<<<< HEAD
=======

>>>>>>> e8e251aee62cccc762e09b76df9b70d4b89adb3d
        df_ds[climate_model] = df_ds_sc
    
    return df_ds