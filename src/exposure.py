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
from scipy import integrate

# load from other py files within pg_risk_analysis
from utils import reproject,set_paths,buffer_assets,overlay_hazard_assets
from extract import extract_osm_infrastructure
from hazard import open_storm_data,open_flood_data

gdal.SetConfigOption("OSM_CONFIG_FILE",os.path.join('/scistor/ivm/mye500/projects/pg_risk_analysis/osmconf.ini'))


def save_exposure(country_code,hazard_type):
    # set paths
    data_path,tc_path,fl_path,osm_data_path,pg_data_path,vul_curve_path,output_path,ne_path = set_paths()
    
    # read infrastructure data:
    osm_lines,osm_poly,osm_points = extract_osm_infrastructure(country_code,osm_data_path)

    if hazard_type=='tc':
        # read wind data
        climate_models = [''] #,'_CMCC-CM2-VHR4','_CNRM-CM6-1-HR','_EC-Earth3P-HR','_HadGEM3-GC31-HM'
        df_ds = open_storm_data(country_code)
        
        # remove assets that will not have any damage
        osm_lines = osm_lines.loc[osm_lines.asset != 'cable'].reset_index(drop=True)
        osm_lines['asset'] = osm_lines['asset'].replace(['minor_line'], 'line')
        osm_poly = osm_poly.loc[osm_poly.asset != 'plant'].reset_index(drop=True)
            
    elif hazard_type=='fl':
        # read flood data
        climate_models = ['historical'] #,'rcp8p5'
        df_ds = open_flood_data(country_code)
                
    for climate_model in climate_models:
        if hazard_type=='tc':
            return_periods = ['1_1{}'.format(climate_model),'1_2{}'.format(climate_model),'1_5{}'.format(climate_model),
                              '1_10{}'.format(climate_model),'1_25{}'.format(climate_model),'1_50{}'.format(climate_model),
                              '1_100{}'.format(climate_model),'1_250{}'.format(climate_model),'1_500{}'.format(climate_model),
                              '1_1000{}'.format(climate_model)]
            
        elif hazard_type == 'fl':
            return_periods = ['rp0001','rp0002','rp0005','rp0010','rp0025','rp0050','rp0100','rp0250','rp0500','rp1000'] 
            
        for return_period in return_periods:
            # assess damage for lines
            #print(df_ds[climate_model][[return_period,'geometry']])

            df_hazard = df_ds[climate_model][[return_period,'geometry']]
            df_hazard = df_hazard[~df_hazard.eq(0).any(axis=1)]
            #print(df_hazard)
            overlay_lines = pd.DataFrame(overlay_hazard_assets(df_hazard,osm_lines).T,
                                         columns=['asset','hazard_point'])

            overlay_lines['geometry'] = None

            for index, row in overlay_lines.iterrows():
                hazard_point = row['hazard_point']
                if hazard_point not in df_hazard.index:
                    overlay_lines = overlay_lines.drop(index)
                else:
                    geometry = df_hazard.loc[hazard_point, 'geometry']
                    overlay_lines.at[index, 'geometry'] = geometry

            #overlay_lines.to_excel(os.path.join(output_path,f'{country_code}_overlay_lines_{climate_model}.xlsx'))

            # assess damage for polygons
            if len(osm_poly) > 0:
                overlay_poly = pd.DataFrame(overlay_hazard_assets(df_hazard,osm_poly).T,
                                        columns=['asset','hazard_point'])
            else:
                overlay_poly = pd.DataFrame()

            overlay_poly['geometry'] = None

            for index, row in overlay_poly.iterrows():
                hazard_point = row['hazard_point']
                if hazard_point not in df_hazard.index:
                    overlay_poly = overlay_poly.drop(index)
                else:
                    geometry = df_hazard.loc[hazard_point, 'geometry']
                    overlay_poly.at[index, 'geometry'] = geometry
                    
            #assess damage for points
            overlay_points = pd.DataFrame(overlay_hazard_assets(df_hazard,osm_points).T,
                                          columns=['asset','hazard_point'])

            overlay_points['geometry'] = None

            for index, row in overlay_points.iterrows():
                hazard_point = row['hazard_point']
                if hazard_point not in df_hazard.index:
                    overlay_points = overlay_points.drop(index)
                else:
                    geometry = df_hazard.loc[hazard_point, 'geometry']
                    overlay_points.at[index, 'geometry'] = geometry
        
            df = pd.concat([overlay_lines,overlay_poly,overlay_points])

            # 根据hazard_point计算每个hazard_point对应的asset个数
            hazard_counts = df.groupby('hazard_point')['asset'].nunique().reset_index()
            hazard_counts.columns = ['hazard_point', 'asset_count']

            # 从原始DataFrame中获取每个hazard_point对应的geometry
            hazard_geometry = df[['hazard_point', 'geometry']].drop_duplicates()

            osm_exposure = pd.merge(hazard_counts, hazard_geometry, on='hazard_point')

            osm_exposure.to_excel(os.path.join(output_path,'exposure',f'{country_code}_osm_exposure_{hazard_type}_{return_period}.xlsx'))

    return osm_exposure


if __name__ == "__main__":
    
    osm_exposure = save_exposure(sys.argv[1],sys.argv[2]) #country_code,hazard_type