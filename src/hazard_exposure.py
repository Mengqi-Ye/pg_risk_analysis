import geopandas as gpd
import pandas as pd
from osgeo import ogr,gdal
import os
import xarray
import numpy as np
import pyproj
from pygeos import from_wkb,from_wkt
import pygeos
from tqdm import tqdm
from pathlib import Path

import extract_osm_data

from multiprocessing import Pool,cpu_count

def define_paths():
    data_path = os.path.join('C:\\','data','pg_risk_analysis')
    netcdf_path = os.path.join(data_path,'tc_netcdf')
    osm_data_path = os.path.join('C:\\','Data','country_osm')

    return data_path,netcdf_path,osm_data_path

def reproject(df_ds,current_crs="epsg:4326",approximate_crs = "epsg:3857"):
    """[summary]

    Args:
        df_ds ([type]): [description]
        current_crs (str, optional): [description]. Defaults to "epsg:3857".
        approximate_crs (str, optional): [description]. Defaults to "epsg:4326".

    Returns:
        [type]: [description]
    """    

    geometries = df_ds['geometry']
    coords = pygeos.get_coordinates(geometries)
    transformer=pyproj.Transformer.from_crs(current_crs, approximate_crs,always_xy=True)
    new_coords = transformer.transform(coords[:, 0], coords[:, 1])
    
    return pygeos.set_coordinates(geometries.copy(), np.array(new_coords).T) 

def open_hazard_data(climate_model,hazard_type='storm'):
    
    data_path,netcdf_path,osm_data_path = define_paths()

    if hazard_type == 'storm':
        with xarray.open_dataset(os.path.join(netcdf_path,'STORM_FIXED_RETURN_PERIODS_WP{}.nc'.format(climate_model))) as ds:

            # get the mean values
            df_ds = ds['mean'].to_dataframe().unstack(level=2).reset_index()

            # create geometry values and drop lat lon columns
            df_ds['geometry'] = [pygeos.points(x) for x in list(zip(df_ds['lon'],df_ds['lat']))]
            df_ds = df_ds.drop(['lat','lon'],axis=1,level=0)

            #rename columns to return periods
            return_periods = ['1_{}{}'.format(int(x),climate_model) for x in ds['rp']]
            df_ds.columns = ['1_{}{}'.format(int(x),climate_model) for x in list(df_ds.columns.get_level_values(1))[:-1]]+['geometry']     
            df_ds['geometry'] = pygeos.buffer(df_ds.geometry,radius=0.1/2,cap_style='square').values
            df_ds['geometry'] = reproject(df_ds)
            # drop all nan values to reduce size
            df_ds = df_ds.loc[~df_ds['1_10000'].isna()].reset_index(drop=True)
            df_ds = df_ds.fillna(0)
                        
    elif hazard_type == 'flood':
        
        # THIS STILL NEEDS TO BE TESTED WITH GLOFRIS DATA
        with xr.open_dataset(os.path.join(flood_path,'SFINCS',hazard_file), engine="rasterio") as ds:
            df_ds = ds.to_dataframe().reset_index()
            df_ds['geometry'] = pygeos.points(df_ds.x,y=df_ds.y)
            df_ds = df_ds.rename(columns={'band_data': 'hazard_intensity'})
            df_ds = df_ds.drop(['band','x', 'y','spatial_ref'], axis=1)
            df_ds = df_ds.dropna()
            df_ds = df_ds.reset_index(drop=True)
            df_ds.geometry= pygeos.buffer(df_ds.geometry,radius=20/2,cap_style='square').values
            df_ds['geometry'] = reproject(df_ds)

    return df_ds

def load_curves_maxdam(data_path): 
    """[summary]

    Args:
        data_path ([type]): [description]

    Returns:
        [type]: [description]
    """
    # load curves and maximum damages as separate inputs
    curves = pd.read_excel(data_path,sheet_name='flooding_curves',skiprows=8,index_col=[0])
    maxdam=pd.read_excel(data_path,sheet_name='flooding_curves',index_col=[0]).iloc[:5]
    
    curves.columns = maxdam.columns

    #transpose maxdam so its easier work with the dataframe
    maxdam = maxdam.T

    #interpolate the curves to fill missing values
    curves = curves.interpolate()
   
    return curves,maxdam

def buffer_assets(assets,buffer_size=100):
    """[summary]

    Args:
        assets ([type]): [description]
        buffer_size (int, optional): [description]. Defaults to 100.

    Returns:
        [type]: [description]
    """    
    assets['buffered'] = pygeos.buffer(assets.geometry.values,buffer_size)# assets.geometry.progress_apply(lambda x: pygeos.buffer(x,buffer_size)) 

    return assets

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
    
def get_damage_per_asset_per_rp(asset,df_ds,assets,curves,maxdam,return_period,country):
    """[summary]

    Args:
        asset ([type]): [description]
        df_ds ([type]): [description]
        assets ([type]): [description]
        grid_size (int, optional): [description]. Defaults to 90.

    Returns:
        [type]: [description]
    """    


    get_hazard_points = df_ds.iloc[asset[1]['hazard_point'].values].reset_index()
    get_hazard_points = get_hazard_points.loc[pygeos.intersects(get_hazard_points.geometry.values,assets.iloc[asset[0]].geometry)]

    asset_type = assets.iloc[asset[0]].asset
    asset_geom = assets.iloc[asset[0]].geometry

    if asset_type in ['plant','substation']:
        maxdam_asset = maxdam.loc[asset_type].MaxDam/pygeos.area(asset_geom)
    else:
        maxdam_asset = maxdam.loc[asset_type].MaxDam


    hazard_intensity = curves[asset_type].index.values
    fragility_values = curves[asset_type].values
    
    if len(get_hazard_points) == 0:
        return asset[0],0
    else:
        
        if pygeos.get_type_id(asset_geom) == 1:
            get_hazard_points['overlay_meters'] = pygeos.length(pygeos.intersection(get_hazard_points.geometry.values,asset_geom))
            return asset[0],np.sum((np.interp(get_hazard_points[return_period].values,hazard_intensity,fragility_values))*get_hazard_points.overlay_meters*maxdam_asset)
        
        elif  pygeos.get_type_id(asset_geom) == 3:
            get_hazard_points['overlay_m2'] = pygeos.area(pygeos.intersection(get_hazard_points.geometry.values,asset_geom))
            return asset[0],get_hazard_points.apply(lambda x: np.interp(x[return_period], hazard_intensity, fragility_values)*maxdam_asset*x.overlay_m2,axis=1).sum()     
        
        else:
            return asset[0],np.sum((np.interp(get_hazard_points[return_period].values,hazard_intensity,fragility_values))*maxdam_asset)     


if __name__ == '__main__':       
    overlay_tc('LAO')
