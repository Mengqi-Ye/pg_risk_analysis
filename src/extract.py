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
from utils import reproject,buffer_assets

gdal.SetConfigOption("OSM_CONFIG_FILE", os.path.join('..',"osmconf.ini"))

# change paths to make it work on your own machine
data_path = os.path.join('C:\\','data','pg_risk_analysis')
tc_path = os.path.join(data_path,'tc_netcdf')
fl_path = os.path.join(data_path,'GLOFRIS')
osm_data_path = os.path.join('C:\\','data','country_osm')
pg_data_path = os.path.join(data_path,'pg_data')

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
        print('OSM query is finished: create Dataframe with assets:')
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

def power_polygon(osm_path):
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
    
    #print(df)
    
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
    power_poly_country = power_polygon(osm_path)
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
    pg_path = os.path.join(pg_data_path,'{}_{}.gpkg'.format(country_code,pg_type)) #e.g.,LAO_line
    pg_data_country = gpd.read_file(os.path.join(pg_path))
    
    pg_data_country = pd.DataFrame(pg_data_country.copy())
    #print(pg_data_country.head())
    pg_data_country.geometry = pygeos.from_shapely(pg_data_country.geometry)
    pg_data_country['geometry'] = reproject(pg_data_country)
    
    if pg_type == 'line':
        pg_data_country = buffer_assets(pg_data_country.loc[pg_data_country.asset.isin(['line'])],buffer_size=100).reset_index(drop=True)
        return pg_data_country
    
    elif pg_type == 'point':
        pg_data_country = buffer_assets(pg_data_country.loc[pg_data_country.asset.isin(['point'])],buffer_size=100).reset_index(drop=True)
        return pg_data_country

    else:
        return pg_data_country