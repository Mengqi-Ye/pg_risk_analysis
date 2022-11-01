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
from hazard import extract_wind_data
from overlay import overlay_hazard_assets

gdal.SetConfigOption("OSM_CONFIG_FILE", os.path.join('..',"osmconf.ini"))

# change paths to make it work on your own machine
data_path = os.path.join('C:\\','data','pg_risk_analysis')
tc_path = os.path.join(data_path,'tc_netcdf')
fl_path = os.path.join(data_path,'GLOFRIS')
osm_data_path = os.path.join('C:\\','data','country_osm')
pg_data_path = os.path.join(data_path,'pg_data')

def load_curves_maxdam(data_path,hazard='wind'): 
    """[summary]

    Args:
        data_path ([type]): [description]

    Returns:
        [type]: [description]
    """
    
    if hazard == 'wind':
        sheet_name = 'flooding_curves'
    elif hazard == 'flood':
        sheet_name = 'flooding_curves'
    
    # load curves and maximum damages as separate inputs
    curves = pd.read_excel(data_path,sheet_name=sheet_name,skiprows=8,index_col=[0])
    maxdam=pd.read_excel(data_path,sheet_name=sheet_name,index_col=[0]).iloc[:5]
    
    curves.columns = maxdam.columns

    #transpose maxdam so its easier work with the dataframe
    maxdam = maxdam.T

    #interpolate the curves to fill missing values
    curves = curves.interpolate()
   
    return curves,maxdam

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

    # find the exact hazard overlays:
    get_hazard_points = df_ds.iloc[asset[1]['hazard_point'].values].reset_index()
    get_hazard_points = get_hazard_points.loc[pygeos.intersects(get_hazard_points.geometry.values,assets.iloc[asset[0]].geometry)]

    asset_type = assets.iloc[asset[0]].asset
    asset_geom = assets.iloc[asset[0]].geometry

    if asset_type in ['plant','substation','generator']:
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
            return asset[0],np.sum((np.interp(get_hazard_points[return_period].values,hazard_intensity,
                                fragility_values))*get_hazard_points.overlay_meters*maxdam_asset)
        
        elif  pygeos.get_type_id(asset_geom) == 3:
            get_hazard_points['overlay_m2'] = pygeos.area(pygeos.intersection(get_hazard_points.geometry.values,asset_geom))
            return asset[0],get_hazard_points.apply(lambda x: np.interp(x[return_period], hazard_intensity, 
                            fragility_values)*maxdam_asset*x.overlay_m2,axis=1).sum()     
        
        else:
            return asset[0],np.sum((np.interp(get_hazard_points[return_period].values,hazard_intensity,fragility_values))*maxdam_asset)

##### ##### ##### ##### ##### ##### ##### #####  
##### ##### ##### OSM DAMAGE  ##### ##### ##### 
##### ##### ##### ##### ##### ##### ##### #####  

def assess_damage_osm(country_code,ctry_power_infra):
    """_summary_

    Args:
        country_code (_type_): _description_
        ctry_power_infra (_type_): _description_

    Returns:
        _type_: _description_
    """    
    # load curves and maxdam
    curves,maxdam = load_curves_maxdam(data_path=os.path.join('..','data','infra_vulnerability_data.xlsx'))
    curves['line'] = 1 # remove this when things work!
    
    # read infrastructure data:
    power_lines,power_poly,power_points = ctry_power_infra
    
    # read wind data
    climate_models = ['_CMCC-CM2-VHR4','_CNRM-CM6-1-HR','_EC-Earth3P-HR','_HadGEM3-GC31-HM']
    df_ds = extract_wind_data()
    
    # calculate damaged lines in loop by country_code and climate_model
    damaged_lines = {}
    for climate_model in climate_models:
        return_periods = ['1_10{}'.format(climate_model),'1_50{}'.format(climate_model),
                          '1_100{}'.format(climate_model),'1_500{}'.format(climate_model),'1_1000{}'.format(climate_model)]

        overlay_lines = pd.DataFrame(overlay_hazard_assets(df_ds[climate_model],
                                                           power_lines).T,columns=['asset','hazard_point'])
        collect_line_damages = []
        for asset in tqdm(overlay_lines.groupby('asset'),total=len(overlay_lines.asset.unique()),
                          desc='polyline damage calculation for {} {}'.format(country_code,climate_model)):
            for return_period in return_periods:
                collect_line_damages.append([return_period,get_damage_per_asset_per_rp(asset,df_ds[climate_model],
                                                                                       power_lines,
                                                                                       curves,
                                                                                       maxdam,
                                                                                       return_period,
                                                                                       country_code)])

        collect_line_damages = [(line[0],line[1][0],line[1][1]) for line in collect_line_damages]
        damaged_lines_country = power_lines.merge(pd.DataFrame(collect_line_damages,
                                                                                 columns=['return_period','index','damage']),
                                                                    left_index=True,right_on='index')
        damaged_lines_country = damaged_lines_country.drop(['buffered'],axis=1)
        damaged_lines[climate_model] = damaged_lines_country

    # calculate damaged polygons in loop by country_code and climate_model
    damaged_poly = {}
    for climate_model in climate_models:
        return_periods = ['1_10{}'.format(climate_model),'1_50{}'.format(climate_model),
                          '1_100{}'.format(climate_model),'1_500{}'.format(climate_model),'1_1000{}'.format(climate_model)]

        overlay_poly = pd.DataFrame(overlay_hazard_assets(df_ds[climate_model],power_poly).T,
                                                          columns=['asset','hazard_point'])
        collect_poly_damages = []
        for asset in tqdm(overlay_poly.groupby('asset'),
                                total=len(overlay_poly.asset.unique()),desc='polygon damage calculation for {} {}'.format(country_code,climate_model)):
            for return_period in return_periods:
                collect_poly_damages.append([return_period,get_damage_per_asset_per_rp(asset,df_ds[climate_model],
                                                                                       power_poly,
                                                                                       curves,maxdam,
                                                                                       return_period,
                                                                                       country_code)])

        collect_poly_damages = [(line[0],line[1][0],line[1][1]) for line in collect_poly_damages]
        damaged_poly_country = power_poly.merge(pd.DataFrame(collect_poly_damages,
                                                columns=['return_period','index','damage']),left_index=True,right_on='index')

        damaged_poly[climate_model] = damaged_poly_country
            
    # calculate damaged points in loop by country_code and climate_model
    damaged_points = {}
    for climate_model in climate_models:
        return_periods = ['1_10{}'.format(climate_model),'1_50{}'.format(climate_model),
                          '1_100{}'.format(climate_model),'1_500{}'.format(climate_model),'1_1000{}'.format(climate_model)]

        overlay_points = pd.DataFrame(overlay_hazard_assets(df_ds[climate_model],power_points).T,
                                      columns=['asset','hazard_point'])
        collect_point_damages = []
        for asset in tqdm(overlay_points.groupby('asset'),
                                    total=len(overlay_points.asset.unique()),desc='point damage calculation for {} {}'.format(country_code,climate_model)):
            for return_period in return_periods:
                collect_point_damages.append([return_period,get_damage_per_asset_per_rp(asset,
                                                                                        df_ds[climate_model],
                                                                                        power_points,
                                                                                        curves,
                                                                                        maxdam,
                                                                                        return_period,
                                                                                        country_code)])

        collect_point_damages = [(line[0],line[1][0],line[1][1]) for line in collect_point_damages]
        damaged_points_country = power_points.merge(pd.DataFrame(collect_point_damages,
                                                        columns=['return_period','index','damage']),left_index=True,right_on='index')
        damaged_points_country = damaged_points_country.drop(['buffered'],axis=1)
        damaged_points[climate_model] = damaged_points_country
                                      
    return damaged_lines,damaged_poly,damaged_points


##### ##### ##### ##### ##### ##### ##### #####  
##### ##### ##### GOV DAMAGE  ##### ##### ##### 
##### ##### ##### ##### ##### ##### ##### #####  

def assess_damage_pg(country_code,pg_data_country):
    """_summary_

    Args:
        country_code (_type_): _description_
        pg_data_country (_type_): _description_

    Returns:
        _type_: _description_
    """    
    # load curves and maxdam
    curves,maxdam = load_curves_maxdam(data_path=os.path.join('..','data','infra_vulnerability_data.xlsx'))
    curves['line'] = 1 # remove this when things work!
    
    # read infrastructure data:
    #power_lines,power_poly,power_points = ctry_power_infra
    
    # read wind data
    climate_models = ['_CMCC-CM2-VHR4','_CNRM-CM6-1-HR','_EC-Earth3P-HR','_HadGEM3-GC31-HM']
    df_ds = extract_wind_data()
    
    # calculate damaged lines in loop by country_code and climate_model
    damaged_lines = {}
    for climate_model in climate_models:
        return_periods = ['1_10{}'.format(climate_model),'1_50{}'.format(climate_model),
                          '1_100{}'.format(climate_model),'1_500{}'.format(climate_model),'1_1000{}'.format(climate_model)]

        overlay_lines = pd.DataFrame(overlay_hazard_assets(df_ds[climate_model],
                                                           pg_data_country).T,columns=['asset','hazard_point'])
        collect_line_damages = []
        for asset in tqdm(overlay_lines.groupby('asset'),total=len(overlay_lines.asset.unique()),
                          desc='polyline damage calculation for {}'.format(country_code,climate_model)):
            for return_period in return_periods:
                collect_line_damages.append([return_period,get_damage_per_asset_per_rp(asset,df_ds[climate_model],
                                                                                       pg_data_country,
                                                                                       curves,
                                                                                       maxdam,
                                                                                       return_period,
                                                                                       country_code)])

        collect_line_damages = [(line[0],line[1][0],line[1][1]) for line in collect_line_damages]
        damaged_lines_country = pg_data_country.merge(pd.DataFrame(collect_line_damages,
                                                                                 columns=['return_period','index','damage']),
                                                                    left_index=True,right_on='index')
        damaged_lines_country = damaged_lines_country.drop(['buffered'],axis=1)
        damaged_lines[climate_model] = damaged_lines_country

    return damaged_lines