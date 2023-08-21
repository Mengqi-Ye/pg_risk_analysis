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
from utils import reproject,set_paths
from extract import extract_osm_infrastructure
from damage_exception_without_scale import assess_damage_osm

gdal.SetConfigOption("OSM_CONFIG_FILE", os.path.join('..',"osmconf.ini"))

# change paths to make it work on your own machine
# Used for China and Japan

##### ##### ##### ##### ##### ##### ##### #####  
##### ##### #####   OSM RISK  ##### ##### ##### 
##### ##### ##### ##### ##### ##### ##### #####

def country_analysis_osm(country_code,hazard_type,climate_model,geo_type):
    # set paths
    data_path,tc_path,fl_path,osm_data_path,pg_data_path,vul_curve_path,output_path,ne_path = set_paths()
    output_path = os.path.join(output_path,'output_without_scaling')
    
    if climate_model == 'present':
        climate_model = ''

    else:
        climate_model = climate_model
        
    # extract infrastructure data from OSM
    osm_power_infra = extract_osm_infrastructure(country_code,osm_data_path)
    
    # assess damage to hazard_type
    osm_damage_infra = assess_damage_osm(country_code,osm_power_infra,hazard_type,climate_model,geo_type)
    
    line_risk = {}
    plant_risk = {}
    substation_risk = {}
    tower_risk = {}
    pole_risk = {}

    df = osm_damage_infra[climate_model]

    curve_code_substation = ['W2_1_1','W2_1_2','W2_1_3','W2_1_4','W2_1_5','W2_1_6','W2_2_1','W2_2_2','W2_2_3','W2_2_4','W2_2_5','W2_2_6',
                             'W2_3_1','W2_3_2','W2_3_3','W2_3_4','W2_3_5','W2_3_6','W2_4_1','W2_4_2','W2_4_3','W2_4_4','W2_4_5','W2_4_6',
                             'W2_5_1','W2_5_2','W2_5_3','W2_5_4','W2_5_5','W2_5_6','W2_6_1','W2_6_2','W2_6_3','W2_6_4','W2_6_5','W2_6_6',
                             'W2_7_1','W2_7_2','W2_7_3','W2_7_4','W2_7_5','W2_7_6']

    curve_code_tower = ['W3_1','W3_2','W3_3','W3_4','W3_5','W3_6','W3_7','W3_8','W3_9','W3_10','W3_11','W3_12','W3_13','W3_14','W3_15',
                        'W3_16','W3_17','W3_18','W3_19','W3_20','W3_21','W3_22','W3_23','W3_24','W3_25','W3_26','W3_27','W3_28']

    curve_code_pole = ['W4_1','W4_2','W4_3','W4_4','W4_5','W4_6','W4_7','W4_8','W4_9','W4_10','W4_11','W4_12',
                       'W4_13','W4_14','W4_15','W4_16','W4_17','W4_18','W4_19','W4_20','W4_21','W4_22','W4_23',
                       'W4_24','W4_25','W4_26','W4_27','W4_28','W4_29','W4_30','W4_31','W4_32','W4_33','W4_34',
                       'W4_35','W4_36','W4_37','W4_38','W4_39','W4_40','W4_41','W4_42','W4_43','W4_44','W4_45',
                       'W4_46','W4_47','W4_48','W4_49','W4_50','W4_51','W4_52','W4_53','W4_54','W4_55','W4_56']

    curve_code_line = ['W5_1_1','W5_1_2','W5_1_3','W5_1_4','W5_1_5','W5_1_6','W5_1_7','W5_1_8','W5_1_9','W5_1_10','W5_1_11','W5_1_12',
                       'W5_2_1','W5_2_2','W5_2_3','W5_2_4','W5_2_5','W5_2_6','W5_2_7','W5_2_8','W5_2_9','W5_2_10','W5_2_11','W5_2_12',
                       'W5_3_1','W5_3_2','W5_3_3','W5_3_4','W5_3_5','W5_3_6','W5_3_7','W5_3_8','W5_3_9','W5_3_10','W5_3_11','W5_3_12',
                       'W5_4_1','W5_4_2','W5_4_3','W5_4_4','W5_4_5','W5_4_6','W5_4_7','W5_4_8','W5_4_9','W5_4_10','W5_4_11','W5_4_12',
                       'W5_5_1','W5_5_2','W5_5_3','W5_5_4','W5_5_5','W5_5_6','W5_5_7','W5_5_8','W5_5_9','W5_5_10','W5_5_11','W5_5_12',
                       'W5_6_1','W5_6_2','W5_6_3','W5_6_4','W5_6_5','W5_6_6','W5_6_7','W5_6_8','W5_6_9','W5_6_10','W5_6_11','W5_6_12',
                       'W5_7_1','W5_7_2','W5_7_3','W5_7_4','W5_7_5','W5_7_6','W5_7_7','W5_7_8','W5_7_9','W5_7_10','W5_7_11','W5_7_12']
            
    #assess risk for power lines
    if geo_type == 'lines':
        if len(df) == 0:
            print("No {}_{} risk of infra_type 'line' in {}".format(hazard_type,climate_model,country_code))

        else:
            with pd.ExcelWriter(os.path.join(output_path,'damage','{}_osm_{}{}_damage_{}'.format(country_code,hazard_type,climate_model,geo_type)+'.xlsx')) as writer:
                df.to_excel(writer)

            df['rp'] = df['rp'].replace(['1_1{}'.format(climate_model),'1_2{}'.format(climate_model),'1_5{}'.format(climate_model),
                                        '1_10{}'.format(climate_model),'1_25{}'.format(climate_model),'1_50{}'.format(climate_model),
                                        '1_100{}'.format(climate_model),'1_250{}'.format(climate_model),'1_500{}'.format(climate_model),
                                        '1_1000{}'.format(climate_model)],
                                        [1,0.5,0.2,0.1,0.04,0.02,0.01,0.004,0.002,0.001])

        for curve_code in curve_code_line:
            loss_list = df.loc[df['curve'] == curve_code]
            loss_list = loss_list.sort_values(by='rp',ascending=False)
            if len(loss_list) == 0:
                print("No risk of power lines ...")
            
            else:
                loss_list_mean = loss_list.meandam.values.tolist()
                loss_list_lower = loss_list.lowerdam.values.tolist()
                loss_list_upper = loss_list.upperdam.values.tolist()
                RPS = loss_list.loc[loss_list['curve'] == curve_code]
                RPS = RPS.rp.values.tolist()
                line_risk[climate_model,curve_code] = {
                    'mean_risk': integrate.simps(y=loss_list_mean[::-1], x=RPS[::-1]),
                    'lower_risk': integrate.simps(y=loss_list_lower[::-1], x=RPS[::-1]),
                    'upper_risk': integrate.simps(y=loss_list_upper[::-1], x=RPS[::-1])
                }

        return pd.DataFrame(line_risk)

    #assess risk for power substations
    elif geo_type == 'poly':         
        if len(df) == 0:
            print("No {}_{} risk of infra_type 'poly' in {}".format(hazard_type,climate_model,country_code))

        else:
            with pd.ExcelWriter(os.path.join(output_path,'damage','{}_osm_{}{}_damage_{}'.format(country_code,hazard_type,climate_model,geo_type)+'.xlsx')) as writer:
                df.to_excel(writer)

            df['rp'] = df['rp'].replace(['1_1{}'.format(climate_model),'1_2{}'.format(climate_model),'1_5{}'.format(climate_model),
                                        '1_10{}'.format(climate_model),'1_25{}'.format(climate_model),'1_50{}'.format(climate_model),
                                        '1_100{}'.format(climate_model),'1_250{}'.format(climate_model),'1_500{}'.format(climate_model),
                                        '1_1000{}'.format(climate_model)],
                                        [1,0.5,0.2,0.1,0.04,0.02,0.01,0.004,0.002,0.001])
            
        for curve_code in curve_code_substation:
            loss_list = df.loc[df['curve'] == curve_code]
            loss_list = loss_list.sort_values(by='rp',ascending=False)
            if len(loss_list) == 0:
                print("No risk of substations ...")
            
            else:
                loss_list_mean = loss_list.meandam.values.tolist()
                loss_list_lower = loss_list.lowerdam.values.tolist()
                loss_list_upper = loss_list.upperdam.values.tolist()
                RPS = loss_list.loc[loss_list['curve'] == curve_code]
                RPS = RPS.rp.values.tolist()
                substation_risk[climate_model,curve_code] = {
                    'mean_risk': integrate.simps(y=loss_list_mean[::-1], x=RPS[::-1]),
                    'lower_risk': integrate.simps(y=loss_list_lower[::-1], x=RPS[::-1]),
                    'upper_risk': integrate.simps(y=loss_list_upper[::-1], x=RPS[::-1])
                }

        return pd.DataFrame(substation_risk)

    #assess risk for power towers and power poles
    elif geo_type == 'points':
        if len(df) == 0:
            print("No {}_{} risk of infra_type 'points' in {}".format(hazard_type,climate_model,country_code))

        else:
            with pd.ExcelWriter(os.path.join(output_path,'damage','{}_osm_{}{}_damage_{}'.format(country_code,hazard_type,climate_model,geo_type)+'.xlsx')) as writer:
                df.to_excel(writer)

            df['rp'] = df['rp'].replace(['1_1{}'.format(climate_model),'1_2{}'.format(climate_model),'1_5{}'.format(climate_model),
                                        '1_10{}'.format(climate_model),'1_25{}'.format(climate_model),'1_50{}'.format(climate_model),
                                        '1_100{}'.format(climate_model),'1_250{}'.format(climate_model),'1_500{}'.format(climate_model),
                                        '1_1000{}'.format(climate_model)],
                                        [1,0.5,0.2,0.1,0.04,0.02,0.01,0.004,0.002,0.001])
            
        for curve_code in curve_code_tower:
            loss_list = df.loc[df['curve'] == curve_code]
            loss_list = loss_list.sort_values(by='rp',ascending=False)
            if len(loss_list) == 0:
                print("No risk of power towers ...")

            else:
                loss_list_mean = loss_list.meandam.values.tolist()
                loss_list_lower = loss_list.lowerdam.values.tolist()
                loss_list_upper = loss_list.upperdam.values.tolist()
                RPS = loss_list.loc[loss_list['curve'] == curve_code]
                RPS = RPS.rp.values.tolist()
                tower_risk[climate_model,curve_code] = {
                    'mean_risk': integrate.simps(y=loss_list_mean[::-1], x=RPS[::-1]),
                    'lower_risk': integrate.simps(y=loss_list_lower[::-1], x=RPS[::-1]),
                    'upper_risk': integrate.simps(y=loss_list_upper[::-1], x=RPS[::-1])
                }

        for curve_code in curve_code_pole:
            loss_list = df.loc[df['curve'] == curve_code]
            loss_list = loss_list.sort_values(by='rp',ascending=False)
            if len(loss_list) == 0:
                print("No risk of power poles ...")

            else:                    
                loss_list_mean = loss_list.meandam.values.tolist()
                loss_list_lower = loss_list.lowerdam.values.tolist()
                loss_list_upper = loss_list.upperdam.values.tolist()
                RPS = loss_list.loc[loss_list['curve'] == curve_code]
                RPS = RPS.rp.values.tolist()
                pole_risk[climate_model,curve_code] = {
                    'mean_risk': integrate.simps(y=loss_list_mean[::-1], x=RPS[::-1]),
                    'lower_risk': integrate.simps(y=loss_list_lower[::-1], x=RPS[::-1]),
                    'upper_risk': integrate.simps(y=loss_list_upper[::-1], x=RPS[::-1])
                }

        return pd.DataFrame(tower_risk),pd.DataFrame(pole_risk)


##### ##### ##### ##### ##### ##### ##### #####  
##### ##### ##### SAVE RESULT ##### ##### ##### 
##### ##### ##### ##### ##### ##### ##### #####
def risk_output(country_code,hazard_type,infra_type,climate_model,geo_type):
    # set paths
    data_path,tc_path,fl_path,osm_data_path,pg_data_path,vul_curve_path,output_path,ne_path = set_paths()
    output_path = os.path.join(output_path,'output_without_scaling')

    if climate_model == 'present':
        climate_model = ''

    else:
        climate_model = climate_model
    
    if climate_model == '':
        writer = pd.ExcelWriter(os.path.join(output_path,'risk','{}_{}_{}_{}_{}_risk'.format(country_code,infra_type,hazard_type,'present',geo_type)+'.xlsx'),
                                engine='openpyxl')
    else:
        writer = pd.ExcelWriter(os.path.join(output_path,'risk','{}_{}_{}{}_{}_risk'.format(country_code,infra_type,hazard_type,climate_model,geo_type)+'.xlsx'),
                                engine='openpyxl')
        
    if geo_type == 'lines':
        line_risk = country_analysis_osm(country_code,hazard_type,climate_model,geo_type)
        if len(line_risk) != 0:
            line_risk[climate_model].to_excel(writer, sheet_name='line_risk')

    elif geo_type == 'poly':
        substation_risk = country_analysis_osm(country_code,hazard_type,climate_model,geo_type)
        if len(substation_risk) != 0:
            substation_risk[climate_model].to_excel(writer, sheet_name='substation_risk')

    elif geo_type == 'points':
        tower_risk,pole_risk = country_analysis_osm(country_code,hazard_type,climate_model,geo_type)    
        if len(tower_risk) != 0:
            tower_risk[climate_model].to_excel(writer, sheet_name='tower_risk')
        if len(pole_risk) != 0:
            pole_risk[climate_model].to_excel(writer, sheet_name='pole_risk')

    # save the Excel file
    if writer.sheets:
        writer.save()
            

if __name__ == "__main__":
    
    save_risk = risk_output(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5]) #country_code,hazard_type,infra_type,climate_model,geo_type
