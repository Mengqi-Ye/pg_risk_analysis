import pickle
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
from collections.abc import Iterable
#import rioxarray

# load from other py files within pg_risk_analysis
from hazard import open_storm_data, open_flood_data
from utils import overlay_hazard_assets,set_paths,flatten


def load_curves_maxdam(country_code,vul_curve_path,hazard_type):
    """Load vulnerability curves and maximum damages for a specific country and hazard type.

    Args:
        country_code (str): Country code for the desired country.
        vul_curve_path (str): Path to the input vulnerability curves file.
        hazard_type (str): Type of hazard ('tc' for tropical cyclone, 'fl' for coastal flooding).

    Returns:
        tuple: A tuple containing two pandas DataFrames:
               - curves: Vulnerability curves.
               - maxdam: Maximum damages.
    """

    # dictionary of GDP per capita ratio for each country
    gdp_ratio = {
        "BRN": {"ratio_usa": 0.5201},
        "KHM": {"ratio_usa": 0.0240},
        "CHN": {"ratio_usa": 0.1772},
        "IDN": {"ratio_usa": 0.0647},
        "JPN": {"ratio_usa": 0.5912},
        "LAO": {"ratio_usa": 0.0434},
        "MYS": {"ratio_usa": 0.1775},
        "MNG": {"ratio_usa": 0.0703},
        "MMR": {"ratio_usa": 0.0276},
        "PRK": {"ratio_usa": 0.0106},
        "PHL": {"ratio_usa": 0.0547},
        "SGP": {"ratio_usa": 1.0091},
        "KOR": {"ratio_usa": 0.5367},
        "TWN": {"ratio_usa": 0.4888},
        "THA": {"ratio_usa": 0.1034},
        "VNM": {"ratio_usa": 0.0573},
        "HKG": {"ratio_usa": 0.7091},
        "MAC": {"ratio_usa": 0.5913}}
    
    if hazard_type == 'tc':
        sheet_name = 'wind_curves'
        
        # load curves and maximum damages as separate inputs
        curves = pd.read_excel(vul_curve_path,sheet_name=sheet_name,skiprows=11)
        
        # dictionary of design wind speeds for each country
        design_wind_speed = {
            "BRN": {"dws": 32},
            "KHM": {"dws": 32},
            "CHN": {"dws": 52},
            "IDN": {"dws": 32},
            "JPN": {"dws": 52},
            "LAO": {"dws": 32},
            "MYS": {"dws": 32},
            "MNG": {"dws": 0},
            "MMR": {"dws": 39},
            "PRK": {"dws": 39},
            "PHL": {"dws": 52},
            "SGP": {"dws": 32},
            "KOR": {"dws": 52},
            "TWN": {"dws": 60},
            "THA": {"dws": 39},
            "VNM": {"dws": 44}}
        dws = design_wind_speed.get(country_code, {}).get("dws", None)
        
        # shift design wind speed of all curves to 60 m/s
        scaling_factor = dws / 60

        curves = curves.apply(lambda x: x * scaling_factor if pd.api.types.is_numeric_dtype(x) else x)
        curves = curves.set_index('Wind speed (m/s)')
        
    elif hazard_type == 'fl':
        sheet_name = 'flooding_curves'    
        
        # load curves and maximum damages as separate inputs
        curves = pd.read_excel(vul_curve_path,sheet_name=sheet_name,skiprows=11,index_col=[0])

    maxdam = pd.read_excel(vul_curve_path,sheet_name=sheet_name,index_col=[0],header=[0,1]).iloc[:8]
    curves.columns = maxdam.columns
    
    #interpolate the curves to fill missing values
    curves = curves.interpolate()
    
    #transpose maxdam so its easier work with the dataframe
    maxdam = maxdam.T
    
    ratio_usa = gdp_ratio.get(country_code, {}).get("ratio_usa", None)

    if ratio_usa is not None:
        print(f"The ratio_usa for {country_code} is {ratio_usa}")
    else:
        print(f"No ratio_usa found for {country_code}")
        
    maxdam['MaxDam'] = maxdam['MaxDam'] * ratio_usa
    maxdam['LowerDam'] = maxdam['LowerDam'] * ratio_usa
    maxdam['UpperDam'] = maxdam['UpperDam'] * ratio_usa

    return curves,maxdam
    
    
def get_damage_per_asset_per_rp(asset,df_ds,assets,curves,maxdam,return_period,country):
    """
    Calculates the damage per asset per return period based on asset type, hazard curves and maximum damage

    Args:
        asset (tuple): Tuple with two dictionaries, containing the asset index and the hazard point index of the asset
        df_ds (pandas.DataFrame): A pandas DataFrame containing hazard points with a 'geometry' column
        assets (geopandas.GeoDataFrame): A GeoDataFrame containing asset geometries and asset type information
        curves (dict): A dictionary with the asset types as keys and their corresponding hazard curves as values
        maxdam (pandas.DataFrame): A pandas DataFrame containing the maximum damage for each asset type
        return_period (str): The return period for which the damage should be calculated
        country (str): The country for which the damage should be calculated

    Returns:
        list or tuple: Depending on the input, the function either returns a list of tuples with the asset index, the curve name and the calculated damage, or a tuple with None,
        None, None if no hazard points are found
    """
    
    # find the exact hazard overlays:
    get_hazard_points = df_ds.iloc[asset[1]['hazard_point'].values].reset_index()
    get_hazard_points = get_hazard_points.loc[pygeos.intersects(get_hazard_points.geometry.values,assets.iloc[asset[0]].geometry)]

    asset_type = assets.iloc[asset[0]].asset
    asset_geom = assets.iloc[asset[0]].geometry

    if asset_type in ['plant','substation','generator']:
        # if plant,substation are points, do not calculate the area
        if pygeos.area(asset_geom) == 0:
            maxdam_asset = maxdam.loc[asset_type].MaxDam
            lowerdam_asset = maxdam.loc[asset_type].LowerDam
            upperdam_asset = maxdam.loc[asset_type].UpperDam
        else:
            maxdam_asset = maxdam.loc[asset_type].MaxDam/pygeos.area(asset_geom)
            lowerdam_asset = maxdam.loc[asset_type].LowerDam/pygeos.area(asset_geom)
            upperdam_asset = maxdam.loc[asset_type].UpperDam/pygeos.area(asset_geom)
    else:
        maxdam_asset = maxdam.loc[asset_type].MaxDam
        lowerdam_asset = maxdam.loc[asset_type].LowerDam
        upperdam_asset = maxdam.loc[asset_type].UpperDam

    hazard_intensity = curves[asset_type].index.values
    
    if isinstance(curves[asset_type],pd.core.series.Series):
        fragility_values = curves[asset_type].values.flatten()
        only_one = True
        curve_name = curves[asset_type].name
    elif len(curves[asset_type].columns) == 1:
        fragility_values = curves[asset_type].values.flatten()      
        only_one = True   
        curve_name = curves[asset_type].columns[0]
    else:
        fragility_values = curves[asset_type].values#.T[0]
        maxdam_asset = maxdam_asset.values#[0]
        only_one = False

    if len(get_hazard_points) == 0:
        if only_one:
            return [return_period,asset[0],curve_name,0,0,0]
        else:
            return [return_period,asset[0],curves[asset_type].columns[0],0,0,0]
            
    else:
        if only_one:    
            # run the calculation as normal when the asset just has a single curve
            if pygeos.get_type_id(asset_geom) == 1:            
                get_hazard_points['overlay_meters'] = pygeos.length(pygeos.intersection(get_hazard_points.geometry.values,asset_geom))
                return [return_period,asset[0],curve_name,np.sum((np.interp(get_hazard_points[return_period].values,hazard_intensity,
                                                             fragility_values))*get_hazard_points.overlay_meters*maxdam_asset),
                                                          np.sum((np.interp(get_hazard_points[return_period].values,hazard_intensity,
                                                             fragility_values))*get_hazard_points.overlay_meters*lowerdam_asset),
                                                          np.sum((np.interp(get_hazard_points[return_period].values,hazard_intensity,
                                                             fragility_values))*get_hazard_points.overlay_meters*upperdam_asset)]

            elif (pygeos.get_type_id(asset_geom) == 3) | (pygeos.get_type_id(asset_geom) == 6) :
                get_hazard_points['overlay_m2'] = pygeos.area(pygeos.intersection(get_hazard_points.geometry.values,asset_geom))
                return [return_period,asset[0],curve_name,get_hazard_points.apply(lambda x: np.interp(x[return_period],hazard_intensity, 
                                                                  fragility_values)*maxdam_asset*x.overlay_m2,axis=1).sum(),
                                                          get_hazard_points.apply(lambda x: np.interp(x[return_period],hazard_intensity, 
                                                                  fragility_values)*lowerdam_asset*x.overlay_m2,axis=1).sum(),
                                                          get_hazard_points.apply(lambda x: np.interp(x[return_period],hazard_intensity, 
                                                                  fragility_values)*upperdam_asset*x.overlay_m2,axis=1).sum()]  

            else:
                return [return_period,asset[0],curve_name,np.sum((np.interp(get_hazard_points[return_period].values,
                                                             hazard_intensity,fragility_values))*maxdam_asset),
                                                          np.sum((np.interp(get_hazard_points[return_period].values,
                                                             hazard_intensity,fragility_values))*lowerdam_asset),
                                                          np.sum((np.interp(get_hazard_points[return_period].values,
                                                             hazard_intensity,fragility_values))*upperdam_asset)]
        else:
            # run the calculation when the asset has multiple curves
            if pygeos.get_type_id(asset_geom) == 1:            
                get_hazard_points['overlay_meters'] = pygeos.length(pygeos.intersection(get_hazard_points.geometry.values,asset_geom))
            elif (pygeos.get_type_id(asset_geom) == 3) | (pygeos.get_type_id(asset_geom) == 6) :
                get_hazard_points['overlay_m2'] = pygeos.area(pygeos.intersection(get_hazard_points.geometry.values,asset_geom))
            
            collect_all = []
            for iter_,curve_ids in enumerate(curves[asset_type].columns):
                if pygeos.get_type_id(asset_geom) == 1:
                    collect_all.append([return_period,asset[0],curves[asset_type].columns[iter_],
                                        np.sum((np.interp(get_hazard_points[return_period].values,
                                                          hazard_intensity,fragility_values.T[iter_]))*get_hazard_points.overlay_meters*maxdam_asset[iter_]),
                                        np.sum((np.interp(get_hazard_points[return_period].values,
                                                          hazard_intensity,fragility_values.T[iter_]))*get_hazard_points.overlay_meters*lowerdam_asset[iter_]),
                                        np.sum((np.interp(get_hazard_points[return_period].values,
                                                          hazard_intensity,fragility_values.T[iter_]))*get_hazard_points.overlay_meters*upperdam_asset[iter_])])
                                   
                elif (pygeos.get_type_id(asset_geom) == 3) | (pygeos.get_type_id(asset_geom) == 6) :
                    collect_all.append([return_period,asset[0],curves[asset_type].columns[iter_],
                                        get_hazard_points.apply(lambda x: np.interp(x[return_period], hazard_intensity,
                                                                                    fragility_values.T[iter_])*maxdam_asset[iter_]*x.overlay_m2,axis=1).sum(),
                                        get_hazard_points.apply(lambda x: np.interp(x[return_period], hazard_intensity,
                                                                                    fragility_values.T[iter_])*lowerdam_asset[iter_]*x.overlay_m2,axis=1).sum(),
                                        get_hazard_points.apply(lambda x: np.interp(x[return_period], hazard_intensity,
                                                                                    fragility_values.T[iter_])*upperdam_asset[iter_]*x.overlay_m2,axis=1).sum()])

                else:
                    collect_all.append([return_period,asset[0],curves[asset_type].columns[iter_],
                                        np.sum((np.interp(get_hazard_points[return_period].values,
                                                          hazard_intensity,fragility_values.T[iter_]))*maxdam_asset[iter_]),
                                        np.sum((np.interp(get_hazard_points[return_period].values,
                                                          hazard_intensity,fragility_values.T[iter_]))*lowerdam_asset[iter_]),
                                        np.sum((np.interp(get_hazard_points[return_period].values,
                                                          hazard_intensity,fragility_values.T[iter_]))*upperdam_asset[iter_])])
            return collect_all



##### ##### ##### ##### ##### ##### ##### #####  
##### ##### ##### OSM DAMAGE  ##### ##### ##### 
##### ##### ##### ##### ##### ##### ##### #####
def assess_damage_osm(country_code,osm_power_infra,hazard_type): #NEW VERSION
    # set paths
    data_path,tc_path,fl_path,osm_data_path,pg_data_path,vul_curve_path,output_path,ne_path = set_paths()

    # load curves and maxdam
    curves,maxdam = load_curves_maxdam(country_code,vul_curve_path,hazard_type)
    
    # read infrastructure data:
    osm_lines,osm_poly,osm_points = osm_power_infra
    
    #calculate damaged lines/polygons/points in loop by climate_model
    damaged_lines = {}
    damaged_poly = {}
    damaged_points = {}

    if hazard_type=='tc':
        # read wind data
        climate_models = ['','_CMCC-CM2-VHR4','_CNRM-CM6-1-HR','_EC-Earth3P-HR','_HadGEM3-GC31-HM']
        df_ds = open_storm_data(country_code)

        # remove assets that will not have any damage
        osm_lines = osm_lines.loc[osm_lines.asset != 'cable'].reset_index(drop=True)
        osm_lines['asset'] = osm_lines['asset'].replace(['minor_line'], 'line')
        osm_poly = osm_poly.loc[osm_poly.asset != 'plant'].reset_index(drop=True)            
    
    elif hazard_type=='fl':
        # read flood data
        climate_models = ['historical','rcp8p5']
        df_ds = open_flood_data(country_code)
        
    for climate_model in climate_models:
        if hazard_type=='tc':
            return_periods = ['1_1{}'.format(climate_model),'1_2{}'.format(climate_model),'1_5{}'.format(climate_model),'1_10{}'.format(climate_model),
                              '1_25{}'.format(climate_model),'1_50{}'.format(climate_model),'1_100{}'.format(climate_model),
                              '1_250{}'.format(climate_model),'1_500{}'.format(climate_model),'1_1000{}'.format(climate_model)]
        elif hazard_type == 'fl':
            return_periods = ['rp0001','rp0002','rp0005','rp0010','rp0025','rp0050','rp0100','rp0250','rp0500','rp1000']     
    
        # assess damage for lines
        overlay_lines = pd.DataFrame(overlay_hazard_assets(df_ds[climate_model],osm_lines).T,
                                     columns=['asset','hazard_point'])

        if len(overlay_lines) == 0:
            damaged_lines[climate_model] = pd.DataFrame()

        else:
            collect_line_damages = []
            for asset in tqdm(overlay_lines.groupby('asset'),total=len(overlay_lines.asset.unique()),
                              desc='polyline damage calculation for {} {} ({})'.format(country_code,hazard_type,climate_model)):
                for return_period in return_periods:
                    collect_line_damages.append(get_damage_per_asset_per_rp(asset,
                                                                            df_ds[climate_model],
                                                                            osm_lines,
                                                                            curves,
                                                                            maxdam,
                                                                            return_period,
                                                                            country_code))

            get_asset_type_line = dict(zip(osm_lines.index,osm_lines.asset))
            
            if hazard_type == 'tc':
                results = pd.DataFrame([item for sublist in collect_line_damages
                                        for item in sublist],columns=['rp','asset','curve','meandam','lowerdam','upperdam'])

                results['asset_type'] = results.asset.apply(lambda x : get_asset_type_line[x])

                damaged_lines[climate_model] = results.groupby(['rp','curve','asset_type']).sum().drop(['asset'], axis=1).reset_index()

            elif hazard_type == 'fl':
                #results = pd.DataFrame(collect_line_damages,columns=['rp','asset','curve','meandam','lowerdam','upperdam'])
                results = pd.DataFrame(np.array(list(flatten(collect_line_damages))).reshape(
                    int(len(list(flatten(collect_line_damages)))/6), 6),
                                       columns=['rp','asset','curve','meandam','lowerdam','upperdam'])
                results['asset'] = results['asset'].astype(int)
                results[['meandam','lowerdam','upperdam']] = results[['meandam','lowerdam','upperdam']].astype(float)
                
                results['asset_type'] = results.asset.apply(lambda x : get_asset_type_line[x])

                #sum damage of line, cable, and minor_line
                results['curve'] = results['curve'].replace(['cable', 'minor_line'], 'line')
                results['asset_type'] = results['asset_type'].replace(['cable', 'minor_line'], 'line')

                damaged_lines[climate_model] = results.groupby(['rp','curve','asset_type']).sum().drop(['asset'], axis=1).reset_index()
                
        # assess damage for polygons
        if len(osm_poly) > 0:
            overlay_poly = pd.DataFrame(overlay_hazard_assets(df_ds[climate_model],osm_poly).T,
                                    columns=['asset','hazard_point'])
        else:
            overlay_poly = pd.DataFrame()

        if len(overlay_poly) == 0:
            damaged_poly[climate_model] = pd.DataFrame()

        else:
            collect_poly_damages = []
            for asset in tqdm(overlay_poly.groupby('asset'),total=len(overlay_poly.asset.unique()),
                              desc='polygon damage calculation for {} {} ({})'.format(country_code,hazard_type,climate_model)):
                for return_period in return_periods:
                    collect_poly_damages.append(get_damage_per_asset_per_rp(asset,
                                                                            df_ds[climate_model],
                                                                            osm_poly,
                                                                            curves,
                                                                            maxdam,
                                                                            return_period,
                                                                            country_code))

            get_asset_type_poly = dict(zip(osm_poly.index,osm_poly.asset))
            
            results = pd.DataFrame([item for sublist in collect_poly_damages 
                                    for item in sublist],columns=['rp','asset','curve','meandam','lowerdam','upperdam'])

            results['asset_type'] = results.asset.apply(lambda x : get_asset_type_poly[x])    

            damaged_poly[climate_model] = results.groupby(['rp','curve','asset_type']).sum().drop(['asset'], axis=1).reset_index()
                

        #assess damage for points
        overlay_points = pd.DataFrame(overlay_hazard_assets(df_ds[climate_model],osm_points).T,
                                      columns=['asset','hazard_point'])

        if len(overlay_points) == 0:
            damaged_points[climate_model] = pd.DataFrame()

        else:
            collect_point_damages = []
            for asset in tqdm(overlay_points.groupby('asset'),total=len(overlay_points.asset.unique()),
                              desc='point damage calculation for {} {} ({})'.format(country_code,hazard_type,climate_model)):
                for return_period in return_periods:
                    
                    # check_error = get_damage_per_asset_per_rp(asset,df_ds[climate_model],osm_points,curves,maxdam,return_period,country_code)
                    # with open(os.path.join(output_path,'get_damage_per_asset_per_rp_{}_{}.pkl'.format(country_code,climate_model)), 'wb') as f:
                    #     pickle.dump(check_error,f)

                    collect_point_damages.append(get_damage_per_asset_per_rp(asset,
                                                                            df_ds[climate_model],
                                                                            osm_points,
                                                                            curves,
                                                                            maxdam,
                                                                            return_period,
                                                                            country_code))

            get_asset_type_point = dict(zip(osm_points.index,osm_points.asset))
            
            if hazard_type == 'tc':
                #catch and remove integers in the sublists of collect_point_damages
                collect_point_damages = [[item for item in sublist if not isinstance(item, int)] for sublist in collect_point_damages]
                collect_point_damages = [[item for item in sublist if len(item) == 6] for sublist in collect_point_damages]

                # with open(os.path.join(output_path,'collect_point_damages_{}_{}.pkl'.format(country_code,climate_model)), 'wb') as f:
                #     pickle.dump(collect_point_damages,f)

                results = pd.DataFrame([item for sublist in collect_point_damages
                                        for item in sublist],columns=['rp','asset','curve','meandam','lowerdam','upperdam'])

                results.drop(results.index[(results.iloc[:, 1] == '_')], inplace=True)
                results.drop(results.index[(results.iloc[:, 1] == '3')], inplace=True)

                results['asset_type'] = results.asset.apply(lambda x : get_asset_type_point[x])
                
                damaged_points[climate_model] = results.groupby(['rp','curve','asset_type']).sum().drop(['asset'], axis=1).reset_index()
            
            elif hazard_type == 'fl':
                results = pd.DataFrame(np.array(list(flatten(collect_point_damages))).reshape(
                    int(len(list(flatten(collect_point_damages)))/6), 6),
                                       columns=['rp','asset','curve','meandam','lowerdam','upperdam'])
                results['asset'] = results['asset'].astype(int)
                results[['meandam','lowerdam','upperdam']] = results[['meandam','lowerdam','upperdam']].astype(float)
                
                #return collect_point_damages,get_asset_type_point
                results['asset_type'] = results.asset.apply(lambda x : get_asset_type_point[x])    

                damaged_points[climate_model] = results.groupby(['rp','curve','asset_type']).sum().drop(['asset'], axis=1).reset_index()

    return damaged_lines,damaged_poly,damaged_points



##### ##### ##### ##### ##### ##### ##### #####  
##### ##### ##### GOV DAMAGE  ##### ##### ##### 
##### ##### ##### ##### ##### ##### ##### #####
def assess_damage_pg(country_code,pg_infra,hazard_type):    
    # set paths
    data_path,tc_path,fl_path,osm_data_path,pg_data_path,vul_curve_path,output_path,ne_path = set_paths()

    # load curves and maxdam
    curves,maxdam = load_curves_maxdam(country_code,vul_curve_path,hazard_type)
    
    # read infrastructure data:
    pg_lines,pg_points = pg_infra
    
    #calculate damaged lines/polygons/points in loop by climate_model
    damaged_lines = {}
    damaged_points = {}

    if hazard_type=='tc':
        # read wind data
        climate_models = ['','_CMCC-CM2-VHR4','_CNRM-CM6-1-HR','_EC-Earth3P-HR','_HadGEM3-GC31-HM']
        df_ds = open_storm_data(country_code)
        
        # remove assets that will not have any damage
        pg_points = pg_points.loc[pg_points.asset != 'plant'].reset_index(drop=True)
    
    elif hazard_type == 'fl':
        # read flood data
        climate_models = ['historical','rcp8p5']
        df_ds = open_flood_data(country_code)
        
    for climate_model in climate_models:
        if hazard_type=='tc':
            return_periods = ['1_1{}'.format(climate_model),'1_2{}'.format(climate_model),'1_5{}'.format(climate_model),'1_10{}'.format(climate_model),
                              '1_25{}'.format(climate_model),'1_50{}'.format(climate_model),'1_100{}'.format(climate_model),
                              '1_250{}'.format(climate_model),'1_500{}'.format(climate_model),'1_1000{}'.format(climate_model)]
        elif hazard_type == 'fl':
            return_periods = ['rp0001','rp0002','rp0005','rp0010','rp0025','rp0050','rp0100','rp0250','rp0500','rp1000'] 
            
        # assess damage for lines
        overlay_lines = pd.DataFrame(overlay_hazard_assets(df_ds[climate_model],pg_lines).T,
                                     columns=['asset','hazard_point'])

        if len(overlay_lines) == 0:
            damaged_lines[climate_model] = pd.DataFrame()

        else:
            collect_line_damages = []
            for asset in tqdm(overlay_lines.groupby('asset'),total=len(overlay_lines.asset.unique()),
                              desc='polyline damage calculation for {} {} ({})'.format(country_code,hazard_type,climate_model)):
                for return_period in return_periods:
                    collect_line_damages.append(get_damage_per_asset_per_rp(asset,
                                                                            df_ds[climate_model],
                                                                            pg_lines,
                                                                            curves,
                                                                            maxdam,
                                                                            return_period,
                                                                            country_code))

            get_asset_type_line = dict(zip(pg_lines.index,pg_lines.asset))
            
            if hazard_type=='tc':
                results = pd.DataFrame([item for sublist in collect_line_damages
                                        for item in sublist],columns=['rp','asset','curve','meandam','lowerdam','upperdam'])
                
                results['asset_type'] = results.asset.apply(lambda x : get_asset_type_line[x])

                damaged_lines[climate_model] = results.groupby(['rp','curve','asset_type']).sum().drop(['asset'], axis=1).reset_index()
                
            elif hazard_type == 'fl':
                results = pd.DataFrame(np.array(list(flatten(collect_line_damages))).reshape(
                    int(len(list(flatten(collect_line_damages)))/6), 6),
                                       columns=['rp','asset','curve','meandam','lowerdam','upperdam'])
                results['asset'] = results['asset'].astype(int)
                results[['meandam','lowerdam','upperdam']] = results[['meandam','lowerdam','upperdam']].astype(float)
                
                results['asset_type'] = results.asset.apply(lambda x : get_asset_type_line[x])

        # assess damage for points
        overlay_points = pd.DataFrame(overlay_hazard_assets(df_ds[climate_model],pg_points).T,
                                      columns=['asset','hazard_point'])

        if len(overlay_points) == 0:
            damaged_points[climate_model] = pd.DataFrame()

        else:
            collect_point_damages = []
            for asset in tqdm(overlay_points.groupby('asset'),total=len(overlay_points.asset.unique()),
                              desc='point damage calculation for {} {} ({})'.format(country_code,hazard_type,climate_model)):
                for return_period in return_periods:
                    collect_point_damages.append(get_damage_per_asset_per_rp(asset,
                                                                            df_ds[climate_model],
                                                                            pg_points,
                                                                            curves,
                                                                            maxdam,
                                                                            return_period,
                                                                            country_code))

            get_asset_type_point = dict(zip(pg_points.index,pg_points.asset))
            
            if hazard_type == 'tc':
                results = pd.DataFrame([item for sublist in collect_point_damages
                                        for item in sublist],columns=['rp','asset','curve','meandam','lowerdam','upperdam'])
            elif hazard_type == 'fl':
                results = pd.DataFrame(np.array(list(flatten(collect_point_damages))).reshape(
                    int(len(list(flatten(collect_point_damages)))/6), 6),
                                       columns=['rp','asset','curve','meandam','lowerdam','upperdam'])
                results['asset'] = results['asset'].astype(int)
                results[['meandam','lowerdam','upperdam']] = results[['meandam','lowerdam','upperdam']].astype(float)

            results['asset_type'] = results.asset.apply(lambda x : get_asset_type_point[x])    

            damaged_points[climate_model] = results.groupby(['rp','curve','asset_type']).sum().drop(['asset'], axis=1).reset_index()

                
    return damaged_lines,damaged_points

