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
from extract_points import extract_osm_infrastructure
from damage_points_without_scale import assess_damage_osm

gdal.SetConfigOption("OSM_CONFIG_FILE", os.path.join('..',"osmconf.ini"))

# change paths to make it work on your own machine

##### ##### ##### ##### ##### ##### ##### #####  
##### ##### #####   OSM RISK  ##### ##### ##### 
##### ##### ##### ##### ##### ##### ##### #####

def country_analysis_osm(country_code,hazard_type,climate_model,geo_type,i):
    # set paths
    data_path,tc_path,fl_path,osm_data_path,pg_data_path,vul_curve_path,output_path,ne_path = set_paths()
    output_path = os.path.join(output_path,'output_without_scaling')

    if climate_model == 'present':
        climate_model = ''

    else:
        climate_model = climate_model
    
    tower_risk = {}
    pole_risk = {}
    
    curve_code_tower = ['W3_1','W3_2','W3_3','W3_4','W3_5','W3_6','W3_7','W3_8','W3_9','W3_10','W3_11','W3_12','W3_13','W3_14','W3_15',
                            'W3_16','W3_17','W3_18','W3_19','W3_20','W3_21','W3_22','W3_23','W3_24','W3_25','W3_26','W3_27','W3_28']

    curve_code_pole = ['W4_1','W4_2','W4_3','W4_4','W4_5','W4_6','W4_7','W4_8','W4_9','W4_10','W4_11','W4_12',
                       'W4_13','W4_14','W4_15','W4_16','W4_17','W4_18','W4_19','W4_20','W4_21','W4_22','W4_23',
                       'W4_24','W4_25','W4_26','W4_27','W4_28','W4_29','W4_30','W4_31','W4_32','W4_33','W4_34',
                       'W4_35','W4_36','W4_37','W4_38','W4_39','W4_40','W4_41','W4_42','W4_43','W4_44','W4_45',
                       'W4_46','W4_47','W4_48','W4_49','W4_50','W4_51','W4_52','W4_53','W4_54','W4_55','W4_56']
    
    # extract infrastructure data from OSM
    #osm_power_infra = extract_osm_infrastructure(country_code,osm_data_path)
    osm_points = extract_osm_infrastructure(country_code,osm_data_path)

    # split osm_points into 10 parts
    num_parts = 10
    split_points = np.array_split(osm_points,num_parts)
    
    tower_risk_df = pd.DataFrame()
    pole_risk_df = pd.DataFrame()
    
    i = int(i)

    # calculate damage for split data
    for points in split_points[i-1:i]:
        osm_damage_infra = assess_damage_osm(country_code,points.reset_index(drop=True),hazard_type,climate_model,geo_type='points')
    
        # assess damage to hazard_type
        #osm_damage_infra = assess_damage_osm(country_code,osm_power_infra,hazard_type,climate_model,geo_type)

        df = osm_damage_infra[climate_model]

        #assess risk for power towers and power poles
        if len(df) == 0:
            print("No {}_{} risk of infra_type 'points' in {}".format(hazard_type,climate_model,country_code))

        else:
            # with pd.ExcelWriter(os.path.join(output_path,'damage','{}_osm_{}{}_damage_{}'.format(country_code,hazard_type,climate_model,geo_type)+'.xlsx')) as writer:
            #     df.to_excel(writer)

            df['rp'] = df['rp'].replace(['1_1{}'.format(climate_model),'1_2{}'.format(climate_model),'1_5{}'.format(climate_model),
                                        '1_10{}'.format(climate_model),'1_25{}'.format(climate_model),'1_50{}'.format(climate_model),
                                        '1_100{}'.format(climate_model),'1_250{}'.format(climate_model),'1_500{}'.format(climate_model),
                                        '1_1000{}'.format(climate_model)],
                                        [1,0.5,0.2,0.1,0.04,0.02,0.01,0.004,0.002,0.001])

        for curve_code in curve_code_tower:
            loss_list = df.loc[df['curve'] == curve_code]
            loss_list = loss_list.sort_values(by='rp',ascending=False)
            if len(loss_list) == 0:
                print("No risk of power towers in {}...".format(i))

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
                print("No risk of power poles in {}...".format(i))

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
                
        tower_risk_df = pd.concat([tower_risk_df, pd.DataFrame(tower_risk)])
        pole_risk_df = pd.concat([pole_risk_df, pd.DataFrame(pole_risk)])

    return tower_risk_df,pole_risk_df


##### ##### ##### ##### ##### ##### ##### #####  
##### ##### ##### SAVE RESULT ##### ##### ##### 
##### ##### ##### ##### ##### ##### ##### #####
def risk_output(country_code,hazard_type,infra_type,climate_model,geo_type,i):
    # set paths
    data_path,tc_path,fl_path,osm_data_path,pg_data_path,vul_curve_path,output_path,ne_path = set_paths()
    output_path = os.path.join(output_path,'output_without_scaling')

    if climate_model == 'present':
        climate_model = ''

    else:
        climate_model = climate_model
    
    if climate_model == '':
        writer = pd.ExcelWriter(os.path.join(output_path,'risk','{}_{}_{}_{}_{}_risk_{}'.format(country_code,infra_type,hazard_type,'present',geo_type,i)+'.xlsx'),
                                engine='openpyxl')
    else:
        writer = pd.ExcelWriter(os.path.join(output_path,'risk','{}_{}_{}{}_{}_risk_{}'.format(country_code,infra_type,hazard_type,climate_model,geo_type,i)+'.xlsx'),
                                engine='openpyxl')
        

    tower_risk,pole_risk = country_analysis_osm(country_code,hazard_type,climate_model,'points',i)    
    if len(tower_risk) != 0:
        tower_risk[climate_model].to_excel(writer, sheet_name='tower_risk')
    if len(pole_risk) != 0:
        pole_risk[climate_model].to_excel(writer, sheet_name='pole_risk')

    # save the Excel file
    if writer.sheets:
        writer.close()
            

if __name__ == "__main__":
    
    save_risk = risk_output(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6]) #country_code,hazard_type,infra_type,climate_model,geo_type,i
