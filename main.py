import osgeo.gdal as gdal
import osgeo.ogr as ogr
from osgeo.gdalconst import *

# Numpy Import
import numpy as np

# Others imports
import os
import zipfile
from datetime import datetime, timedelta
import json
import tifffile as tiff


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


def download_sentinel2(tile, start_date, end_date):
    """
    Download every Sentinel-2 L2A images between two dates.
    @params:
        tile   	                    - Required  : tilename
        start_date                  - Required  : start date
        end_date                    - Required  : end date
    """
    command = "python ./sentinel_download/theia_download.py -t T" + str(tile) + " -c SENTINEL2 -a ./sentinel_download/config_theia.cfg + -d " + str(start_date) + " -f " + str(end_date) + " -m 10"
    os.system(command)


def download_sentinel1(tiles_list, directory_images, step=8):
    """
    Download every Sentinel-2 L2A images between two dates.
    @params:
        tiles_list   	        - Required  : list of Sentinel-2 tiles
        directory_images        - Required  : directory containing Sentinel-2 images
        step                    - Optional  : day interval between Sentinel-2 image to download Sentinel-1 images
    """
    s1_output = './S1_img'
    for tile in tiles_list:
        s1_tile = tile

        for img in os.listdir(directory_images):
            if s1_tile in img:
                _date = get_date_from_image_path(img)
                datetime_object = datetime.strptime(_date, '%Y%m%d')
                s1_first_date = (datetime_object - timedelta(days=int(step))).strftime('%Y-%m-%d')
                s1_last_date = (datetime_object + timedelta(days=int(step))).strftime('%Y-%m-%d')

                replacements = {'{output_path}': s1_output, '{first_date}': s1_first_date, '{last_date}': s1_last_date,
                                '{tiles}': s1_tile}

                with open('./sentinel_1_tiling/S1Processor.cfg') as infile, open(
                        './sentinel_1_tiling/S1Processor_tmp.cfg', 'w') as outfile:
                    for line in infile:
                        for src, target in replacements.iteritems():
                            line = line.replace(src, target)
                        outfile.write(line)

                command = 'python ./sentinel_1_tiling/S1Processor.py ./sentinel_1_tiling/S1Processor_tmp.cfg'
                os.system(command)
                os.system('rm ./sentinel_1_tiling/S1Processor_tmp.cfg')


def check_satellite_data(tiles, directory_s2, directory_s1, step=8, filtered=True):
    """
    Construct JSON file to create database of S1 and S2 images.
    @params:
        tiles   	        - Required  : list of Sentinel-2 tiles
        directory_s2        - Required  : directory containing Sentinel-2 images
        directory_s1        - Required  : directory containing Sentinel-1 images pretreated from s1tiling
        step                - Optional  : day interval between Sentinel-2 image to download Sentinel-1 images
        filtered            - Optional  : using or not filtered images
    """
    data = {}

    for tile in tiles:
        data[tile] = []
        for img_s2 in os.listdir(directory_s2):
            if tile in img_s2:
                date = get_date_from_image_path(img_s2)
                s1_path_vv = get_s1_from_s2_tile(tile, date, directory_s1, 'vv', step=step, filtered=filtered)
                s1_path_vh = get_s1_from_s2_tile(tile, date, directory_s1, 'vh', step=step, filtered=filtered)
                data[tile].append({
                    'dateS2': date,
                    's2_path': os.path.join(directory_s2, img_s2),
                    's1_path_vv': s1_path_vv,
                    's1_path_vh': s1_path_vh
                })

    if os.path.exists('./sentinel_data.json'):
        os.remove('./sentinel_data.json')

    with open('./sentinel_data.json', 'w') as file:
        json.dump(data, file)

    return data


def get_s1_from_s2_tile(tile, date_s2, directory_s1, pol, step=8, filtered=True):
    """
    Extract filename of the Sentinel-1 image
    @params:
        tile   	            - Required  : current tile
        date_s2             - Required  : date of current Sentinel-2
        directory_s1        - Required  : directory containing Sentinel-1 images pretreated from s1tiling
        step                - Optional  : day interval between Sentinel-2 image to download Sentinel-1 images
        filtered            - Optional  : using or not filtered images
    """
    path_folder = os.path.join(directory_s1, str(tile))
    if filtered:
        path_folder = os.path.join(directory_s1, str(tile), 'filtered')

    s1_files = []
    datetime_s2 = datetime.strptime(date_s2, '%Y%m%d')
    s1_first_date = datetime_s2 - timedelta(days=int(step))
    s1_last_date = datetime_s2 + timedelta(days=int(step))

    for file in os.listdir(path_folder):
        if (tile in file) and (pol in file):
            datetime_current_s1 = datetime.strptime(get_date_from_image_path(file, s1=True), '%Y%m%d')
            if s1_first_date <= datetime_current_s1 <= s1_last_date:
                s1_files.append(file)

    _size = None
    res_img = None
    for file in s1_files:
        _img = tiff.imread(os.path.join(path_folder, file))
        if _size is None:
            _size = np.count_nonzero(_img == 0)
            res_img = file
        else:
            if np.count_nonzero(_img == 0) < _size:
                res_img = file

    return os.path.join(path_folder, res_img)


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

    print(datetime.now().isoformat() + 'Successfully extracted ' + str(path_to_zip_file))


def get_date_from_image_path(name, s1=False):
    """
    Return date from Sentinel-2 folder name or Sentinel-1 pretreated name.
    @params:
        path   	                    - Required  : Folder name containing every image bands.
        s1   	                    - Optional  : True if name is Sentinel-1 filename
    """
    if s1:
        date = name.split('_')[5][0:8]
    else:
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
    start_date = '2020-06-15'
    end_date = '2020-06-30'

    print(datetime.now().isoformat() + ' Downloading S2 images from ' + str(
        start_date) + ' to ' + str(end_date) + ' for each S2 tiles.')
    for tile in tiles_list:
        print('')
        #download_sentinel2(tile, start_date, end_date)
    print(datetime.now().isoformat() + ' Successfully downloaded S2 images from ' + str(
        start_date) + ' to ' + str(end_date) + '.')

    print(datetime.now().isoformat() + ' Unzipping S2 images downloaded.')
    for file in os.listdir(cwd):
        file_path = os.path.join(cwd, file)
        if zipfile.is_zipfile(file_path):
            unzip(file_path, directory_images)
    print(datetime.now().isoformat() + ' Successfully unzipped S2 images.')

    print(datetime.now().isoformat() + ' Rasterizing ground reference for each S2 tiles.')
    if not os.path.exists(directory_ground_reference_raster):
        os.mkdir(directory_ground_reference_raster)

    for tile in tiles_list:
        for img in os.listdir(directory_images):
            if tile in img:
                path_band_ref = os.path.join(directory_images, img, img + '_SRE_B3.tif')
                output_filename = os.path.join(cwd, directory_ground_reference_raster, tile + '_ground_reference.tif')
                #rasterize_ground_ref(path_band_ref, path_shp_ground_reference, attribute_shp_ground_reference, output_filename)
    print(datetime.now().isoformat() + ' Successfully rasterized ground reference for each S2 tiles.')

    print(datetime.now().isoformat() + ' Downloading and calibrating Sentinel-1 data.')
    download_sentinel1(tiles_list, directory_images, step=8)
    print(datetime.now().isoformat() + ' Successfully downloaded and calibrated Sentinel-1 data.')

    print(datetime.now().isoformat() + ' Checking Sentinel-1 data.')
    s1_output = './S1_img'
    json = check_satellite_data(tiles_list, directory_images, s1_output)
    print(datetime.now().isoformat() + ' Successfully checked Sentinel-1 data.')

    print(datetime.now().isoformat() + ' Calculating patches for each tile.')
    if not os.path.exists(directory_dataset):
        os.mkdir(directory_dataset)

    for tile in tiles_list:
        calculate_patches_tile(tile, directory_ground_reference_raster, directory_images, directory_dataset)
    print(datetime.now().isoformat() + ' Successfully calculated patches for each tile.')

    print(datetime.now().isoformat() + ' Checking overlapped patches.')

    print(datetime.now().isoformat() + ' Successfully checked and deleted overlapped patches')


if __name__ == '__main__':
    main()
