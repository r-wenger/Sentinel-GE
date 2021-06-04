import osgeo.gdal as gdal
import osgeo.ogr as ogr
import osgeo.osr as osr
from osgeo.gdalconst import *

# Numpy Import
import numpy as np
import numpy.ma as ma

# Others imports
import random
import time
import os, shutil, sys
import zipfile
import datetime
import fileinput

def rasterize_ground_ref(path_raster_reference, path_shapefile_to_raster, attribute_shapefile, output_filename):
    """
    Rasterizing ground truth.
    @params:
        path_raster_reference   	- Required  : raster or satellite image as reference
        path_shapefile_to_raster   	- Required  : shapefile to rasterize
        attribute_shapefile   		- Required  : attribute to rasterize (name of shapefile column)
        output_filename   			- Required  : output file rasterized
    """
    dataset_raster = gdal.Open(path_raster_reference, GA_ReadOnly)
    geo_transform = dataset_raster.GetGeoTransform()

    dataset_shapefile = ogr.Open(path_shapefile_to_raster)
    layer_shapefile = dataset_shapefile.GetLayer()

    nodata_value = 255
    # source_layer = data.GetLayer()

    x_min = geo_transform[0]
    y_max = geo_transform[3]
    x_max = x_min + geo_transform[1] * dataset_raster.RasterXSize
    y_min = y_max + geo_transform[5] * dataset_raster.RasterYSize
    x_res = dataset_raster.RasterXSize
    y_res = dataset_raster.RasterYSize
    pixel_width = geo_transform[1]

    target_ds = gdal.GetDriverByName('GTiff').Create(output_filename, x_res, y_res, 1, gdal.GDT_UInt16)
    target_ds.SetGeoTransform(geo_transform)
    target_ds.SetProjection(dataset_raster.GetProjection())
    band = target_ds.GetRasterBand(1)
    band.SetNoDataValue(nodata_value)
    band.FlushCache()
    band.Fill(nodata_value)

    gdal.RasterizeLayer(target_ds, [1], layer_shapefile, options=["ATTRIBUTE=" + str(attribute_shapefile)])

    target_ds = None


def get_sentinel_tiles(path_shapefile_sentinel, field):
    """
    Returning an array containing all tilenames in the shapefile tiling.
    @params:
        path_shapefile_sentinel   	- Required  : path to the shapefile containing tiling
        field                      	- Required  : column name in the shapefile to extract tile number
    """
    driver = ogr.GetDriverByName("ESRI Shapefile")
    data_source = driver.Open(path_shapefile_sentinel, 0)
    layer = data_source.GetLayer()

    res = []

    for feature in layer:
        _val = feature.GetField(field)
        res.append(_val)

    return res


def download_sentinel(tile, start_date, end_date):
    """
    Download every Sentinel-2 L2A images between two dates.
    @params:
        tile   	                    - Required  : tilename
        start_date                  - Required  : start date
        end_date                    - Required  : end date
    """
    command = "python ./sentinel_download/theia_download.py -t T" + str(tile) + " -c SENTINEL2 -a ./sentinel_download/config_theia.cfg + -d " + str(start_date) + " -f " + str(end_date) + " -m 10"
    #os.system(command)


def unzip(path_to_zip_file, directory_to_extract_to):
    """
    Unzip a file to a precise directory and remove it after process.
    @params:
        path_to_zip_file   	        - Required  : Path to zip file
        directory_to_extract_to     - Required  : Path to the directory where zip file should be extracted
    """
    if not os.path.exists(directory_to_extract_to):
        os.mkdir(directory_to_extract_to)

    command_unzip = 'unzip -o ' + str(path_to_zip_file) + ' -d ' + str(directory_to_extract_to)
    os.system(command_unzip)
    command_del = 'rm ' + str(path_to_zip_file)
    os.system(command_del)

    print(datetime.datetime.now().isoformat() + 'Successfully extracted ' + str(path_to_zip_file))


def get_date_from_image_path(name):
    """
    Return date from Sentinel-2 folder name.
    @params:
        path   	                    - Required  : Folder name containing every image bands.
    """
    date = name.split('_')[1].split('-')[0]

    return date


def calculate_patches_tile(tile, directory_ground_reference, directory_sat_images, directory_dataset):
    return 0


def main():
    cwd = os.getcwd()
    path_shapefile_sentinel = './inputs/S2_Tiling_GE_test.shp'
    field = 'Name'
    directory_images = os.path.join(cwd, 'S2_img')
    directory_ground_reference_raster = os.path.join(cwd, 'ground_reference_raster')
    directory_dataset = os.path.join(cwd, 'dataset')
    path_shp_ground_reference = './inputs/merge_OCSGE.shp'
    attribute_shp_ground_reference = 'cod_n4'

    #Read all S2-Tiles and put them in a, array
    tiles_list = get_sentinel_tiles(path_shapefile_sentinel, field)

    #Download Sentinel for each tile
    start_date = '2018-07-01'
    end_date = '2018-07-30'

    print(datetime.datetime.now().isoformat() + ' Downloading S2 images from ' + str(
        start_date) + ' to ' + str(end_date) + ' for each S2 tiles.')
    for tile in tiles_list:
        download_sentinel(tile, start_date, end_date)
    print(datetime.datetime.now().isoformat() + ' Successfully downloaded S2 images from ' + str(
        start_date) + ' to ' + str(end_date) + '.')

    print(datetime.datetime.now().isoformat() + ' Unzipping S2 images downloaded.')
    for file in os.listdir(cwd):
        file_path = os.path.join(cwd, file)
        if zipfile.is_zipfile(file_path):
            unzip(file_path, directory_images)
    print(datetime.datetime.now().isoformat() + ' Successfully unzipped S2 images.')

    print(datetime.datetime.now().isoformat() + ' Rasterizing ground reference for each S2 tiles.')
    if not os.path.exists(directory_ground_reference_raster):
        os.mkdir(directory_ground_reference_raster)

    for tile in tiles_list:
        for img in os.listdir(directory_images):
            if tile in img:
                path_band_ref = os.path.join(directory_images, img, img + '_SRE_B3.tif')
                output_filename = os.path.join(cwd, directory_ground_reference_raster, tile + '_ground_reference.tif')
                rasterize_ground_ref(path_band_ref, path_shp_ground_reference, attribute_shp_ground_reference,
                                     output_filename)
    print(datetime.datetime.now().isoformat() + ' Successfully rasterized ground reference for each S2 tiles.')

    print(datetime.datetime.now().isoformat() + ' Downloading and calibrating Sentinel-1 data.')

    print(datetime.datetime.now().isoformat() + ' Successfully downloaded and calibrated Sentinel-1 data.')

    print(datetime.datetime.now().isoformat() + ' Calculating patches for each tile.')
    if not os.path.exists(directory_dataset):
        os.mkdir(directory_dataset)

    for tile in tiles_list:
        calculate_patches_tile(tile, directory_ground_reference_raster, directory_images, directory_dataset)
    print(datetime.datetime.now().isoformat() + ' Successfully calculated patches for each tile.')

    print(datetime.datetime.now().isoformat() + ' Checking overlapped patches.')

    print(datetime.datetime.now().isoformat() + ' Successfully checked and deleted overlapped patches')


if __name__ == '__main__':
    main()
