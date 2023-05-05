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

gdal.SetConfigOption("OSM_CONFIG_FILE", os.path.join('..',"osmconf.ini"))

def clip_flood_data(country_code):
    # set paths
    data_path,tc_path,fl_path,osm_data_path,pg_data_path,vul_curve_path,output_path,ne_path = set_paths()

    # load country geometry file and create geometry to clip
    ne_countries = gpd.read_file(ne_path)
    geometry = ne_countries.loc[ne_countries['ISO_A3']==country_code].geometry.values[0]
    geoms = [mapping(geometry)]
    
    #climate_model: historical, rcp4p5, rcp8p5; time_period: hist, 2030, 2050, 2080
    rps = ['0001','0002','0005','0010','0025','0050','0100','0250','0500','1000']
    climate_models = ['historical','rcp8p5']
    
    for rp in rps:
        #global input_file
        for climate_model in climate_models:
            if climate_model=='historical':
                input_file = os.path.join(fl_path,'global',
                                          'inuncoast_{}_nosub_hist_rp{}_0.tif'.format(climate_model,rp)) 
 
            elif climate_model=='rcp8p5':
                input_file = os.path.join(fl_path,'global',
                                          'inuncoast_{}_nosub_2030_rp{}_0.tif'.format(climate_model,rp))

            # load raster file and save clipped version
            with rasterio.open(input_file) as src:
                out_image, out_transform = mask(src, geoms, crop=True)
                out_meta = src.meta

                out_meta.update({"driver": "GTiff",
                         "height": out_image.shape[1],
                         "width": out_image.shape[2],
                         "transform": out_transform})

                if 'scistor' in fl_path:
                    file_path = os.path.join(fl_path,'country','_'.join([country_code]+input_file.split('_')[6:]))
                else:
                    file_path = os.path.join(fl_path,'country','_'.join([country_code]+input_file.split('_')[3:]))

                with rasterio.open(file_path, "w", **out_meta) as dest:
                    dest.write(out_image)
