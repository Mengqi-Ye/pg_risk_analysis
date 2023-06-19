# FIX ERROR IN IDN AND MYS
import os,sys
os.environ['USE_PYGEOS'] = '0'
import geopandas as gpd
import pandas as pd
from osgeo import ogr,gdal
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
import rioxarray
import matplotlib.pyplot as plt
from scipy import integrate
from collections.abc import Iterable
import openpyxl
from openpyxl import Workbook

import warnings
warnings.filterwarnings("ignore")

from scipy import integrate
from hazard import open_storm_data, open_flood_data
from utils import overlay_hazard_assets,set_paths,flatten
from extract import extract_osm_infrastructure,extract_pg_infrastructure
from damage import assess_damage_osm,assess_damage_pg,load_curves_maxdam,get_damage_per_asset_per_rp

data_path,tc_path,fl_path,osm_data_path,pg_data_path,vul_curve_path,output_path,ne_path = set_paths()

# load curves and maxdam
curves,maxdam = load_curves_maxdam('MYS',vul_curve_path,'tc')

# read infrastructure data:
osm_lines,osm_poly,osm_points = extract_osm_infrastructure('MYS',osm_data_path)
print(type(osm_points))

#calculate damaged lines/polygons/points in loop by climate_model
damaged_points = {}

# read wind data
climate_models = ['']
df_ds = open_storm_data('MYS')
    
for climate_model in climate_models:
    return_periods = ['1_1{}'.format(climate_model),'1_2{}'.format(climate_model),
                      '1_5{}'.format(climate_model),'1_10{}'.format(climate_model),
                      '1_25{}'.format(climate_model),'1_50{}'.format(climate_model),'1_100{}'.format(climate_model),
                      '1_250{}'.format(climate_model),'1_500{}'.format(climate_model),'1_1000{}'.format(climate_model)]

    #assess damage for points
    overlay_points = pd.DataFrame(overlay_hazard_assets(df_ds[''],osm_points).T,
                                  columns=['asset','hazard_point'])

    if len(overlay_points) == 0:
        damaged_points[climate_model] = pd.DataFrame()

    else:
        collect_point_damages = []
        for asset in tqdm(overlay_points.groupby('asset'),total=len(overlay_points.asset.unique()),
                          desc='point damage calculation for {} {} ({})'.format('MYS','tc','')):
            for return_period in return_periods:
                collect_point_damages.append(get_damage_per_asset_per_rp(asset,
                                                                        df_ds[''],
                                                                        osm_points,
                                                                        curves,
                                                                        maxdam,
                                                                        return_period,
                                                                        'MYS'))

        get_asset_type_point = dict(zip(osm_points.index,osm_points.asset))

        collect_point_damages = [[item for item in sublist if not isinstance(item, int)] for sublist in collect_point_damages]

        results = pd.DataFrame([item for sublist in collect_point_damages
                                for item in sublist],columns=['rp','asset','curve','meandam','lowerdam','upperdam'])
        
        print(get_asset_type_point)
        print(get_asset_type_point.keys())
        
        if '_' in get_asset_type_point:
            print("Key '_' exists")
        else:
            print("Key '_' doesn't exit")
        
        if '_' in results.columns:
            print("'_'' is present as a column in 'results_asset' DataFrame.")
        else:
            print("'_'' is not present as a column in 'results_asset' DataFrame.")



        results['asset_type'] = results.asset.apply(lambda x : get_asset_type_point[x])

        damaged_points[''] = results.groupby(['rp','curve','asset_type']).sum().drop(['asset'], axis=1).reset_index()
