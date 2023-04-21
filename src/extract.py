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
import rioxarray
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings("ignore")

# load from other py files within pg_risk_analysis
from utils import reproject,buffer_assets,set_paths

gdal.SetConfigOption("OSM_CONFIG_FILE", os.path.join('..',"osmconf.ini"))

##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### ##### ##### Extract OSM data  ##### ##### ##### 
##### ##### ##### ##### ##### ##### ##### ##### ##### 

def query_b(geoType,keyCol,**valConstraint):
    """
    This function builds an SQL query from the values passed to the retrieve() function.
    Arguments:
         *geoType* : Type of geometry (osm layer) to search for.
         *keyCol* : A list of keys/columns that should be selected from the layer.
         ***valConstraint* : A dictionary of constraints for the values. e.g. WHERE 'value'>20 or 'value'='constraint'
    Returns:
        *string: : a SQL query string.
    """
    query = "SELECT " + "osm_id"
    for a in keyCol: query+= ","+ a  
    query += " FROM " + geoType + " WHERE "
    # If there are values in the dictionary, add constraint clauses
    if valConstraint: 
        for a in [*valConstraint]:
            # For each value of the key, add the constraint
            for b in valConstraint[a]: query += a + b
        query+= " AND "
    # Always ensures the first key/col provided is not Null.
    query+= ""+str(keyCol[0]) +" IS NOT NULL" 
    return query 


def retrieve(osm_path,geoType,keyCol,**valConstraint):
    """
    Function to extract specified geometry and keys/values from OpenStreetMap
    Arguments:
        *osm_path* : file path to the .osm.pbf file of the region 
        for which we want to do the analysis.     
        *geoType* : Type of Geometry to retrieve. e.g. lines, multipolygons, etc.
        *keyCol* : These keys will be returned as columns in the dataframe.
        ***valConstraint: A dictionary specifiying the value constraints.  
        A key can have multiple values (as a list) for more than one constraint for key/value.  
    Returns:
        *GeoDataFrame* : a geopandas GeoDataFrame with all columns, geometries, and constraints specified.    
    """
    driver=ogr.GetDriverByName('OSM')
    data = driver.Open(osm_path)
    query = query_b(geoType,keyCol,**valConstraint)
    sql_lyr = data.ExecuteSQL(query)
    features =[]
    # cl = columns 
    cl = ['osm_id'] 
    for a in keyCol: cl.append(a)
    if data is not None:
        print('query is finished, lets start the loop')
        for feature in tqdm(sql_lyr,desc='extract'):
            #try:
            if feature.GetField(keyCol[0]) is not None:
                geom1 = (feature.geometry().ExportToWkt())
                #print(geom1)
                geom = from_wkt(feature.geometry().ExportToWkt()) 
                if geom is None:
                    continue
                # field will become a row in the dataframe.
                field = []
                for i in cl: field.append(feature.GetField(i))
                field.append(geom)   
                features.append(field)
            #except:
            #    print("WARNING: skipped OSM feature")   
    else:
        print("ERROR: Nonetype error when requesting SQL. Check required.")    
    cl.append('geometry')                   
    if len(features) > 0:
        return pd.DataFrame(features,columns=cl)
    else:
        print("WARNING: No features or No Memory. returning empty GeoDataFrame") 
        return pd.DataFrame(columns=['osm_id','geometry'])

def power_polyline(osm_path):
    """
    Function to extract all energy linestrings from OpenStreetMap  
    Arguments:
        *osm_path* : file path to the .osm.pbf file of the region 
        for which we want to do the analysis.        
    Returns:
        *GeoDataFrame* : a geopandas GeoDataFrame with specified unique energy linestrings.
    """
    df = retrieve(osm_path,'lines',['power','voltage'])
    
    df = df.reset_index(drop=True).rename(columns={'power': 'asset'})
    
    #print(df) #check infra keys
    
    return df.reset_index(drop=True)

def power_polygon(osm_path): # check with joel, something was wrong here with extracting substations
    """
    Function to extract energy polygons from OpenStreetMap  
    Arguments:
        *osm_path* : file path to the .osm.pbf file of the region 
        for which we want to do the analysis.        
    Returns:
        *GeoDataFrame* : a geopandas GeoDataFrame with specified unique energy linestrings.
    """
    df = retrieve(osm_path,'multipolygons',['other_tags']) 
    
    df = df.loc[(df.other_tags.str.contains('power'))]   #keep rows containing power data         
    df = df.reset_index(drop=True).rename(columns={'other_tags': 'asset'})     
    
    df['asset'].loc[df['asset'].str.contains('"power"=>"substation"', case=False)]  = 'substation' #specify row
    df['asset'].loc[df['asset'].str.contains('"power"=>"plant"', case=False)] = 'plant' #specify row
    
    df = df.loc[(df.asset == 'substation') | (df.asset == 'plant')]
            
    return df.reset_index(drop=True) 

def electricity(osm_path):
    """
    Function to extract building polygons from OpenStreetMap    
    Arguments:
        *osm_path* : file path to the .osm.pbf file of the region 
        for which we want to do the analysis.        
    Returns:
        *GeoDataFrame* : a geopandas GeoDataFrame with all unique building polygons.    
    """
    df = retrieve(osm_path,'multipolygons',['power'])
    
    df = df.reset_index(drop=True).rename(columns={'power': 'asset'})
    
    #df = df[df.asset!='generator']
    df['asset'].loc[df['asset'].str.contains('"power"=>"substation"', case=False)]  = 'substation' #specify row
    df['asset'].loc[df['asset'].str.contains('"power"=>"plant"', case=False)] = 'plant' #specify row
    
    #print(df)  #check infra keys
    
    df = df.loc[(df.asset == 'substation') | (df.asset == 'plant')]
    
    return df.reset_index(drop=True)

def retrieve_poly_subs(osm_path, w_list, b_list):
    """
    Function to extract electricity substation polygons from OpenStreetMap
    Arguments:
        *osm_path* : file path to the .osm.pbf file of the region
        for which we want to do the analysis.
        *w_list* :  white list of keywords to search in the other_tags columns
        *b_list* :  black list of keywords of rows that should not be selected
    Returns:
        *GeoDataFrame* : a geopandas GeoDataFrame with specified unique substation.
    """
    df = retrieve(osm_path,'multipolygons',['other_tags'])
    df = df[df.other_tags.str.contains('substation', case=False, na=False)]
    #df = df.loc[(df.other_tags.str.contains('substation'))]
    df = df[~df.other_tags.str.contains('|'.join(b_list))]
    #df = df.reset_index(drop=True).rename(columns={'other_tags': 'asset'})
    df['asset']  = 'substation' #specify row
    #df = df.loc[(df.asset == 'substation')] #specify row
    return df.reset_index(drop=True)

def power_point(osm_path):
    """
    Function to extract energy points from OpenStreetMap  
    Arguments:
        *osm_path* : file path to the .osm.pbf file of the region 
        for which we want to do the analysis.        
    Returns:
        *GeoDataFrame* : a geopandas GeoDataFrame with specified unique energy linestrings.
    """   
    df = retrieve(osm_path,'points',['other_tags']) 
    df = df.loc[(df.other_tags.str.contains('power'))]  #keep rows containing power data       
    df = df.reset_index(drop=True).rename(columns={'other_tags': 'asset'})     
        
    df['asset'].loc[df['asset'].str.contains('"power"=>"tower"', case=False)]  = 'power_tower' #specify row
    df['asset'].loc[df['asset'].str.contains('"power"=>"pole"', case=False)] = 'power_pole' #specify row
    #df['asset'].loc[df['asset'].str.contains('"utility"=>"power"', case=False)] = 'power_tower' #specify row
    
    df = df.loc[(df.asset == 'power_tower') | (df.asset == 'power_pole')]
            
    return df.reset_index(drop=True)

def extract_osm_infrastructure(country_code,osm_data_path):
    """_summary_

    Args:
        country_code (_type_): _description_
        osm_data_path (_type_): _description_

    Returns:
        _type_: _description_
    """

    # lines
    osm_path = os.path.join(osm_data_path,'{}.osm.pbf'.format(country_code))
    power_lines_country = power_polyline(osm_path)
    power_lines_country['geometry'] = reproject(power_lines_country)
    power_lines_country = buffer_assets(power_lines_country.loc[power_lines_country.asset.isin(
        ['cable','minor_cable','line','minor_line'])],buffer_size=100).reset_index(drop=True)
    
    # polygons
    osm_path = os.path.join(osm_data_path,'{}.osm.pbf'.format(country_code))
    power_poly_country = electricity(osm_path)
    power_poly_country['geometry'] = reproject(power_poly_country)
    
    # points
    osm_path = os.path.join(osm_data_path,'{}.osm.pbf'.format(country_code))
    power_points_country = power_point(osm_path)
    power_points_country['geometry'] = reproject(power_points_country)
    power_points_country = buffer_assets(power_points_country.loc[power_points_country.asset.isin(
        ['power_tower','power_pole'])],buffer_size=100).reset_index(drop=True)


    return power_lines_country,power_poly_country,power_points_country


##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### ##### ##### Extract Gov data  ##### ##### ##### 
##### ##### ##### ##### ##### ##### ##### ##### ##### 

def extract_pg_data(country_code,pg_type):
    """_summary_

    Args:
        country_code (_type_): _description_
        pg_type (_type_): _description_

    Returns:
        _type_: _description_
    """

    # set paths
    data_path,tc_path,fl_path,osm_data_path,pg_data_path,vul_curve_path,output_path,ne_path = set_paths()
    
    files = [x for x in os.listdir(pg_data_path)  if country_code in x ]
    
    if pg_type=='line':
        for file in files: 
            file_path = os.path.join(pg_data_path,'{}_{}.gpkg'.format(country_code,pg_type))

            pg_data_country = gpd.read_file(file_path)
            pg_data_country = pd.DataFrame(pg_data_country.copy())
            #print(pg_data_country.head())
            pg_data_country.geometry = pygeos.from_shapely(pg_data_country.geometry)
            pg_data_country['geometry'] = reproject(pg_data_country)

        pg_data_country = buffer_assets(pg_data_country.loc[pg_data_country.asset.isin(['line'])],
                                        buffer_size=100).reset_index(drop=True)

    elif pg_type=='point':
        for file in files:
            file_path = os.path.join(pg_data_path,'{}_{}.gpkg'.format(country_code,pg_type))
                
            pg_data_country = gpd.read_file(file_path)
            pg_data_country = pd.DataFrame(pg_data_country.copy())
            pg_data_country.geometry = pygeos.from_shapely(pg_data_country.geometry)
            pg_data_country['geometry'] = reproject(pg_data_country)
            #print(pg_data_country)

        pg_data_country = buffer_assets(pg_data_country.loc[pg_data_country.asset.isin(['plant_point','substation_point','power_tower','power_pole'])],
                                        buffer_size=100).reset_index(drop=True)

    return pg_data_country

def open_pg_data(country_code):
    pg_lines = extract_pg_data(country_code,'line')
    pg_points = extract_pg_data(country_code,'point')

    return pg_lines,pg_points