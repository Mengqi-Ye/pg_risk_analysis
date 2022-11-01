import pygeos
import pyproj
import numpy as np

def reproject(df_ds,current_crs="epsg:4326",approximate_crs = "epsg:3857"):
    """[summary]

    Args:
        df_ds ([type]): [description]
        current_crs (str, optional): [description]. Defaults to "epsg:3857".
        approximate_crs (str, optional): [description]. Defaults to "epsg:4326".

    Returns:
        [type]: [description]
    """    

    geometries = df_ds['geometry']
    coords = pygeos.get_coordinates(geometries)
    transformer=pyproj.Transformer.from_crs(current_crs, approximate_crs,always_xy=True)
    new_coords = transformer.transform(coords[:, 0], coords[:, 1])
    
    return pygeos.set_coordinates(geometries.copy(), np.array(new_coords).T) 

def buffer_assets(assets,buffer_size=100):
    """[summary]

    Args:
        assets ([type]): [description]
        buffer_size (int, optional): [description]. Defaults to 100.

    Returns:
        [type]: [description]
    """    
    assets['buffered'] = pygeos.buffer(assets.geometry.values,buffer_size)
    return assets