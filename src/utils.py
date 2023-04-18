import os
import pygeos
import pyproj
import numpy as np

def set_paths(local=False): # set local=False if run on cluster

    if local == True:
        data_path = os.path.join('C:\\','Data','pg_risk_analysis')
        tc_path = os.path.join(data_path,'tc_netcdf')
        fl_path = os.path.join(data_path,'GLOFRIS')
        osm_data_path = os.path.join('C:\\','Data','country_osm')
        pg_data_path = os.path.join(data_path,'pg_data')
        vul_curve_path = os.path.join(data_path,'vulnerability_curves','input_vulnerability_data.xlsx')
        output_path = os.path.join('C:\\','projects','pg_risk_analysis','output')
    else:
        base_path = "/scistor/ivm/data_catalogue/open_street_map"
        data_path = os.path.join(base_path,'pg_risk_analysis')
        tc_path = os.path.join(data_path,'tc_netcdf')
        fl_path = os.path.join(data_path,'GLOFRIS')
        osm_data_path = os.path.join(base_path,'country_osm')
        pg_data_path = os.path.join(data_path,'pg_data')
        vul_curve_path = os.path.join(data_path,'vulnerability_curves','input_vulnerability_data.xlsx')
        #output_path = os.path.join('C:\\','projects','pg_risk_analysis','output')
        output_path = os.path.join("/scistor/ivm/mye500/projects/output")

    return data_path,tc_path,fl_path,osm_data_path,pg_data_path,vul_curve_path,output_path


def reproject(df_ds, current_crs="epsg:4326", approximate_crs="epsg:3857"):
    """
    Reproject a GeoPandas DataFrame from one CRS to another.

    Parameters
    ----------
    df_ds : geopandas.GeoDataFrame
        The input GeoPandas DataFrame to reproject.
    current_crs : str, optional
        The current CRS of the input geometry column, by default "epsg:4326".
    approximate_crs : str, optional
        The target CRS to reproject to, by default "epsg:3857".

    Returns
    -------
    geopandas.GeoSeries
        The reprojected geometry column as a GeoPandas GeoSeries.
    """

    # Extract the input geometries as a numpy array of coordinates
    geometries = df_ds['geometry']
    coords = pygeos.get_coordinates(geometries)

    # Transform the coordinates using pyproj
    transformer = pyproj.Transformer.from_crs(current_crs, approximate_crs, always_xy=True)
    new_coords = transformer.transform(coords[:, 0], coords[:, 1])

    # Create a new GeoSeries with the reprojected coordinates
    return pygeos.set_coordinates(geometries.copy(), np.array(new_coords).T)

def buffer_assets(assets, buffer_size=100):
    """
    Create a buffer of a specified size around the geometries in a GeoDataFrame.
    
    Args:
        assets (GeoDataFrame): A GeoDataFrame containing geometries to be buffered.
        buffer_size (int, optional): The distance in the units of the GeoDataFrame's CRS to buffer the geometries.
            Defaults to 100.
    
    Returns:
        GeoDataFrame: A new GeoDataFrame with an additional column named 'buffered' containing the buffered
            geometries.
    """
    # Create a buffer of the specified size around the geometries
    assets['buffered'] = pygeos.buffer(assets.geometry.values, buffer_size)
    
    return assets


def overlay_hazard_assets(df_ds, assets):
    """
    Overlay a set of assets with a hazard dataset and return the subset of assets that intersect with
    one or more hazard polygons or lines.
    
    Args:
        df_ds (GeoDataFrame): A GeoDataFrame containing the hazard dataset.
        assets (GeoDataFrame): A GeoDataFrame containing the assets to be overlaid with the hazard dataset.
    
    Returns:
        ndarray: A numpy array of integers representing the indices of the hazard geometries that intersect with
            the assets. If the assets have a 'buffered' column, the buffered geometries are used for the overlay.
    """
    hazard_tree = pygeos.STRtree(df_ds.geometry.values)
    if (pygeos.get_type_id(assets.iloc[0].geometry) == 3) | (pygeos.get_type_id(assets.iloc[0].geometry) == 6):
        return  hazard_tree.query_bulk(assets.geometry,predicate='intersects')    
    else:
        return  hazard_tree.query_bulk(assets.buffered,predicate='intersects')