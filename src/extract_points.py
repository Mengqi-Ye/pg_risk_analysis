import geopandas as gpd
import pandas as pd
from osgeo import ogr,gdal
import os
import xarray as xr
import numpy as np
from pygeos import from_wkb,from_wkt
import pygeos
from tqdm import tqdm
from shapely.wkb import loads
from pathlib import Path
from shapely.geometry import mapping
pd.options.mode.chained_assignment = None
from rasterio.mask import mask
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

def power_point(osm_path):
    """
    Function to extract energy points from OpenStreetMap  
    Arguments:
        *osm_path* : file path to the .osm.pbf file of the region 
        for which we want to do the analysis.        
    Returns:
        *GeoDataFrame* : a geopandas GeoDataFrame with specified unique energy points.
    """   
    df = retrieve(osm_path,'points',['other_tags'])
    df = df.loc[(df.other_tags.str.contains('power'))]  #keep rows containing power data       
    df = df.reset_index(drop=True).rename(columns={'other_tags': 'asset'})
        
    df['asset'].loc[df['asset'].str.contains('"power"=>"tower"', case=False)]  = 'power_tower' #specify row
    df['asset'].loc[df['asset'].str.contains('"power"=>"pole"', case=False)] = 'power_pole' #specify row
    
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
    # set paths
    data_path,tc_path,fl_path,osm_data_path,pg_data_path,vul_curve_path,output_path,ne_path = set_paths()
    
    # points
    osm_path = os.path.join(osm_data_path,'{}.osm.pbf'.format(country_code))
    osm_points = power_point(osm_path)
    osm_points['geometry'] = reproject(osm_points)
    osm_points = buffer_assets(osm_points.loc[osm_points.asset.isin(
        ['power_tower','power_pole'])],buffer_size=100).reset_index(drop=True)
    
    return osm_points