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

##### ##### ##### ##### ##### ##### ##### #####  
##### ##### #####   OSM RISK  ##### ##### ##### 
##### ##### ##### ##### ##### ##### ##### #####  
def country_analysis_osm(country_code,hazard_type): 
    # set paths
    data_path,tc_path,fl_path,osm_data_path,pg_data_path,vul_curve_path,output_path = set_paths()
    
    # assess damage to hazard_type
    osm_damage_infra = assess_damage_osm(country_code,hazard_type)
    
    if hazard_type=='tc':
        climate_models = ['','_CMCC-CM2-VHR4','_CNRM-CM6-1-HR','_EC-Earth3P-HR','_HadGEM3-GC31-HM']
    
    elif hazard_type=='fl':
        climate_models = ['historical','rcp8p5']

    line_risk = {}
    plant_risk = {}
    substation_risk = {}
    tower_risk = {}
    pole_risk = {}
    
    for i in range(len(osm_damage_infra)):
        for climate_model in climate_models:
            if climate_model in osm_damage_infra[i]:
                df = osm_damage_infra[i][climate_model]

                with pd.ExcelWriter(os.path.join(output_path,'damage','{}_osm_{}_{}_damage_{}'.format(country_code,hazard_type,climate_model,i)+'.xlsx')) as writer:
                    df.to_excel(writer)

                df['rp'] = df['rp'].replace(['rp0001','rp0002','rp0005','rp0010','rp0025','rp0050','rp0100','rp0250','rp0500','rp1000'],
                                            [1,0.5,0.2,0.1,0.04,0.02,0.01,0.004,0.002,0.001])

                #assess risk for power lines
                if i == 0:
                    loss_list = df.meandam.values.tolist()
                    RPS = df.rp.values.tolist()
                    line_risk[climate_model] = integrate.simps(y=loss_list[::-1], x=RPS[::-1])

                #assess risk for power plants and substations
                elif i == 1:
                    loss_list = df.loc[df['asset_type'] == 'plant']
                    loss_list = loss_list.meandam.values.tolist()
                    RPS = df.loc[df['asset_type'] == 'plant']
                    RPS = RPS.rp.values.tolist()
                    plant_risk[climate_model] = integrate.simps(y=loss_list[::-1], x=RPS[::-1])

                    loss_list = df.loc[df['asset_type'] == 'substation']
                    loss_list = loss_list.meandam.values.tolist()
                    RPS = df.loc[df['asset_type'] == 'substation']
                    RPS = RPS.rp.values.tolist()
                    substation_risk[climate_model] = integrate.simps(y=loss_list[::-1], x=RPS[::-1])

                #assess risk for power towers and power poles
                elif i == 2:
                    loss_list = df.loc[df['asset_type'] == 'power_tower']
                    loss_list = loss_list.meandam.values.tolist()
                    RPS = df.loc[df['asset_type'] == 'power_tower']
                    RPS = RPS.rp.values.tolist()
                    tower_risk[climate_model] = integrate.simps(y=loss_list[::-1], x=RPS[::-1])

                    loss_list = df.loc[df['asset_type'] == 'power_pole']
                    loss_list = loss_list.meandam.values.tolist()
                    RPS = df.loc[df['asset_type'] == 'power_pole']
                    RPS = RPS.rp.values.tolist()
                    pole_risk[climate_model] = integrate.simps(y=loss_list[::-1], x=RPS[::-1])
            
            else:
                print("No risk: {} {}".format(country_code,climate_model))
                
    return line_risk,plant_risk,substation_risk,tower_risk,pole_risk



##### ##### ##### ##### ##### ##### ##### #####  
##### ##### #####   GOV RISK  ##### ##### ##### 
##### ##### ##### ##### ##### ##### ##### #####  
def country_analysis_pg()


##### ##### ##### ##### ##### ##### ##### #####  
##### ##### ##### SAVE RESULT ##### ##### ##### 
##### ##### ##### ##### ##### ##### ##### #####  
def risk_output(country_code,hazard_type,infra_type):
    
    if infra_type == 'osm':
        total_risk = pd.DataFrame(country_analysis_osm(country_code,hazard_type))
        total_risk.to_excel(os.path.join(output_path,'risk','{}_{}_{}_risk'.format(country_code,infra_type,hazard_type)+'.xlsx'))
    
    elif infra_type == 'gov':
        total_risk = pd.DataFrame(country_analysis_pg(country_code,hazard_type))
        total_risk.to_excel(os.path.join(output_path,'risk','{}_{}_{}_risk'.format(country_code,infra_type,hazard_type)+'.xlsx'))
    
    return total_risk


if __name__ == "__main__":
    
    osm_damage_infra = risk_output(sys.argv[1],sys.argv[2],sys.argv[3]) #country_code, hazard_type,infra_type