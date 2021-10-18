import osgeo.gdal as gdal
import osgeo.ogr as ogr
import osgeo.osr as osr
from osgeo.gdalconst import *

# Numpy Import
import numpy as np
import numpy.ma as ma
from skimage import morphology

# Others imports
import os
import zipfile
from datetime import datetime, timedelta
from scipy.ndimage import zoom
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


def change_projection(path_raster_reference, path_shapefile_to_raster, attribute_shapefile_n5, attribute_shapefile_n4, output_filename,
                      attribute_values_n5=[], attribute_values_n4=[]):
    """
    Changing shapefile projection to match with raster.
    @params:
        path_raster_reference   	- Required  : raster or satellite image as reference
        path_shapefile_to_raster   	- Required  : shapefile to rasterize
        attribute_shapefile   		- Required  : attribute to rasterize (name of shapefile column)
        output_filename   			- Required  : output file rasterized
        attribute_values   			- Required  : list of attributes to rasterize
    """
    if (len(attribute_values_n5) == 0):
        raise Exception('Error : length of attributes array is 0')

    tif = gdal.Open(path_raster_reference)
    driver = ogr.GetDriverByName("ESRI Shapefile")

    dataSource = driver.Open(path_shapefile_to_raster, 0)

    layer = dataSource.GetLayer()

    sourceprj = layer.GetSpatialRef()
    targetprj = osr.SpatialReference(wkt=tif.GetProjection())
    transform = osr.CoordinateTransformation(sourceprj, targetprj)

    to_fill = ogr.GetDriverByName("ESRI Shapefile")
    ds = to_fill.CreateDataSource(output_filename)
    outlayer = ds.CreateLayer('', targetprj, ogr.wkbPolygon)
    outlayer.CreateField(ogr.FieldDefn('ID', ogr.OFTInteger))
    outlayer.CreateField(ogr.FieldDefn(attribute_shapefile_n5, ogr.OFTInteger))

    i = 0

    for feature in layer:
        _val = None

        if (feature.GetField(attribute_shapefile_n5) is not None) and (feature.GetField(attribute_shapefile_n5) != 0) :
            if int(feature.GetField(attribute_shapefile_n5)) in attribute_values_n5:
                _val = feature.GetField(attribute_shapefile_n5)
            else:
                _val = 0
            transformed = feature.GetGeometryRef()
            transformed.Transform(transform)

            geom = ogr.CreateGeometryFromWkb(transformed.ExportToWkb())
            defn = outlayer.GetLayerDefn()
            feat = ogr.Feature(defn)
            feat.SetField('ID', i)
            feat.SetField(attribute_shapefile_n5, _val)
            feat.SetGeometry(geom)
            outlayer.CreateFeature(feat)
            i += 1
            feat = None
        else:
            if int(feature.GetField(attribute_shapefile_n4)) in attribute_values_n4:
                _val = feature.GetField(attribute_shapefile_n4)
            else:
                _val = 0
            transformed = feature.GetGeometryRef()
            transformed.Transform(transform)

            geom = ogr.CreateGeometryFromWkb(transformed.ExportToWkb())
            defn = outlayer.GetLayerDefn()
            feat = ogr.Feature(defn)
            feat.SetField('ID', i)
            feat.SetField(attribute_shapefile_n5, _val)
            feat.SetGeometry(geom)
            outlayer.CreateFeature(feat)
            i += 1
            feat = None

    ds = None


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
    command = "python ./sentinel_download/theia_download.py -t T" + str(
        tile) + " -c SENTINEL2 -a ./sentinel_download/config_theia.cfg + -d " + str(start_date) + " -f " + str(
        end_date) + " -m 10"
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


def merge_classes(input_raster, dic_classes, output_filename):
    """
    Merging some classes togerther to get a new typology
    @params:
        input_raster   	- Required  : raster to merge
        dic_classes   	- Required  : dictionnary containing original classes as keys, and new typologies as value
        output_filename - Required  : output path to create the new raster
    """
    ds_raster = gdal.Open(input_raster, GA_ReadOnly)
    array_raster = ds_raster.GetRasterBand(1).ReadAsArray()
    tab_return = None

    for key in dic_classes:
        if tab_return is None:
            tab_return = np.where(array_raster == key, dic_classes[key], 0)
        else:
            tab_return += np.where(array_raster == key, dic_classes[key], 0)

    if -99 in tab_return:
        tab_return = ma.masked_array(tab_return, tab_return == -99)
        for shift in (-3, 3):
            for axis in (0, 1):
                a_shifted = np.roll(tab_return, shift=shift, axis=axis)
                idx = ~a_shifted.mask * tab_return.mask
                tab_return[idx] = a_shifted[idx]

    se_disk = morphology.disk(3)
    se_rectangle = morphology.rectangle(3, 3)
    tab_return = morphology.opening(tab_return, se_disk)
    tab_return = morphology.opening(tab_return, se_rectangle)

    driver = gdal.GetDriverByName("GTiff")
    dst_ds = driver.Create(output_filename, ds_raster.RasterXSize, ds_raster.RasterYSize, 1, gdal.GDT_Byte)
    dst_ds.SetProjection(ds_raster.GetProjection())
    geotransform = ds_raster.GetGeoTransform()
    dst_ds.SetGeoTransform(geotransform)
    dst_ds.GetRasterBand(1).WriteArray(tab_return)
    dst_ds = None


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


def get_coordinates(dataset, x_index_topleft, y_index_topleft):
    """
    Get upper left coordinates using matrix cell number
    @params:
        dataset   		  - Required  : dataset raster
        x_index_topleft   - Required  : x cell index
        y_index_topleft   - Required  : y cell index
    """
    if dataset is not None:
        (upper_left_x, x_size, x_rotation, upper_left_y, y_rotation, y_size) = dataset.GetGeoTransform()

        x_coords_topleft = x_index_topleft * x_size + upper_left_x
        y_coords_topleft = y_index_topleft * y_size + upper_left_y

        return (x_coords_topleft, y_coords_topleft, x_size, y_size)
    else:
        raise Exception('Error : dataset is None')


def resize_bands(dataser_raster):
    gt = dataser_raster.GetGeoTransform()
    res_x = int(gt[1])
    res_y = int(abs(gt[-1]))

    _res_band = None
    if res_x != 10 or res_y != 10:
        ratio = float(res_x / 10)
        to_resample = dataser_raster.GetRasterBand(1).ReadAsArray()
        _res_band = np.array(zoom(to_resample, ratio, order=0))
    else:
        _res_band = np.array(dataser_raster.GetRasterBand(1).ReadAsArray())

    return _res_band


def get_s1_dates(path):
    res = []

    assert len(os.listdir(path)) > 2

    for img in os.listdir(path):
        _c = img.split('_')[5][0:8]
        if _c not in res:
            res.append(_c)

    return res


def create_folders_dataset(path_directory_dataset):
    if os.path.exists(os.path.join(path_directory_dataset, 's1')):
        os.remove(os.path.join(path_directory_dataset, 's1'))
        os.mkdir(os.path.join(path_directory_dataset, 's1'))
    else:
        os.mkdir(os.path.join(path_directory_dataset, 's1'))

    if os.path.exists(os.path.join(path_directory_dataset, 's2')):
        os.remove(os.path.join(path_directory_dataset, 's2'))
        os.mkdir(os.path.join(path_directory_dataset, 's2'))
    else:
        os.mkdir(os.path.join(path_directory_dataset, 's2'))

    if os.path.exists(os.path.join(path_directory_dataset, 'ground_reference')):
        os.remove(os.path.join(path_directory_dataset, 'ground_reference'))
        os.mkdir(os.path.join(path_directory_dataset, 'ground_reference'))
    else:
        os.mkdir(os.path.join(path_directory_dataset, 'ground_reference'))

    if os.path.exists(os.path.join(path_directory_dataset, 'labels')):
        os.remove(os.path.join(path_directory_dataset, 'labels'))
        os.mkdir(os.path.join(path_directory_dataset, 'labels'))
    else:
        os.mkdir(os.path.join(path_directory_dataset, 'labels'))


# Nomenclature des patchs : tile_sensor_date_topX_topY_
def calculate_patches_tile(tile, directory_ground_reference, directory_sat_images_s2, directory_sat_images_s1,
                           directory_dataset, patch_size, method=1):
    bands_S2 = ['B2.tif', 'B3.tif', 'B4.tif', 'B5.tif', 'B6.tif', 'B7.tif', 'B8.tif', 'B8A.tif', 'B11.tif', 'B12.tif']
    polar_S1 = ['vv', 'vh']
    path_gr = None

    for gr in os.listdir(directory_ground_reference):
        if tile in gr:
            path_gr = os.path.join(directory_ground_reference, gr)

    assert path_gr is not None

    ds_gr_raster = gdal.Open(path_gr)
    gr_raster = ds_gr_raster.GetRasterBand(1).ReadAsArray()

    if method == 1:
        step = patch_size
    else:
        step = 1

    for y in range(0, gr_raster.shape[0], step):
        for x in range(0, gr_raster.shape[1], step):
            if method == 2:
                step = 1

            gr_patch = gr_raster[y:y + patch_size, x:x + patch_size]

            if 0 not in gr_patch:
                if method == 2:
                    step = patch_size

                x_coords_topleft, y_coords_topleft, x_size, y_size = get_coordinates(ds_gr_raster,
                                                                                     x, y)
                new_geotrans = (x_coords_topleft, x_size, 0, y_coords_topleft, 0, y_size)

                driver = gdal.GetDriverByName("GTiff")

                for img_dir in os.listdir(directory_sat_images_s2):
                    if tile in img_dir:
                        date_s2 = get_date_from_image_path(img_dir)
                        output_patch_name_s2 = os.path.join(directory_dataset, 's2', str(tile) + '_' + str(date_s2) + '_S2_' +
                                                            str(x) + '_' + str(y) + '.tif')
                        output_patch_name_gr = os.path.join(directory_dataset, 'ground_reference', str(tile) + '_GR_' +
                                                            str(x) + '_' + str(y) + '.tif')

                        dst_ds_img_s2 = driver.Create(output_patch_name_s2, patch_size, patch_size, len(bands_S2),
                                                      gdal.GDT_Float32, options=["COMPRESS=LZW"])
                        dst_ds_img_s2.SetProjection(ds_gr_raster.GetProjection())
                        dst_ds_img_s2.SetGeoTransform(new_geotrans)

                        if not os.path.exists(output_patch_name_gr):
                            dst_ds_mask = driver.Create(output_patch_name_gr, patch_size, patch_size, 1, gdal.GDT_UInt16,
                                                    options=["COMPRESS=LZW"])
                            dst_ds_mask.SetProjection(ds_gr_raster.GetProjection())
                            dst_ds_mask.SetGeoTransform(new_geotrans)
                            dst_ds_mask.GetRasterBand(1).WriteArray(gr_patch)

                        count_s2_bands = 1
                        for band_name in bands_S2:
                            for band in os.listdir(os.path.join(directory_sat_images_s2, img_dir)):
                                if ('SRE' in band) and (band_name in band):
                                    ds_band = gdal.Open(os.path.join(directory_sat_images_s2, img_dir, band))
                                    patch_band = resize_bands(ds_band)[y:y + patch_size, x:x + patch_size]
                                    dst_ds_img_s2.GetRasterBand(count_s2_bands).WriteArray(patch_band)
                                    count_s2_bands += 1

                for folder_tile in os.listdir(directory_sat_images_s1):
                    if str(folder_tile) == str(tile):
                        for year in os.listdir(os.path.join(directory_sat_images_s1, folder_tile)):
                            for month in os.listdir(os.path.join(directory_sat_images_s1, folder_tile, year)):
                                complete_dates_s1 = get_s1_dates(os.path.join(directory_sat_images_s1, folder_tile,
                                                                              year, month))

                                for date in complete_dates_s1:
                                    count_s1_polar = 1
                                    date_s1 = str(year) + str(month) + str(date[-2:])
                                    output_patch_name_s1 = os.path.join(directory_dataset, 's1',
                                                                        str(tile) + '_' + str(date_s1) + '_S1_' +
                                                                        str(x) + '_' + str(y) + '.tif')
                                    dst_ds_img_s1 = driver.Create(output_patch_name_s1, patch_size, patch_size,
                                                                  len(polar_S1),
                                                                  gdal.GDT_Float32, options=["COMPRESS=LZW"])
                                    dst_ds_img_s1.SetProjection(ds_gr_raster.GetProjection())
                                    dst_ds_img_s1.SetGeoTransform(new_geotrans)

                                    for polar in polar_S1:
                                        for img in os.listdir(
                                                os.path.join(directory_sat_images_s1, folder_tile, year, month)):
                                            if (polar in img) and (date in img):
                                                path_radar = os.path.join(directory_sat_images_s1, folder_tile, year,
                                                                          month, img)
                                                ds_band = gdal.Open(path_radar)
                                                patch_band = resize_bands(ds_band)[y:y + patch_size, x:x + patch_size]
                                                dst_ds_img_s1.GetRasterBand(count_s1_polar).WriteArray(patch_band)
                                                count_s1_polar += 1


def get_geometry(raster):
    """
    Calculate geometry of a raster
    @params:
        raster   		  - Required  : dataset raster
    """
    transform = raster.GetGeoTransform()
    pixelWidth = transform[1]
    pixelHeight = transform[5]
    cols = raster.RasterXSize
    rows = raster.RasterYSize

    xLeft = transform[0]
    yTop = transform[3]
    xRight = xLeft + cols * pixelWidth
    yBottom = yTop - rows * pixelHeight

    ring = ogr.Geometry(ogr.wkbLinearRing)
    ring.AddPoint(xLeft, yTop)
    ring.AddPoint(xLeft, yBottom)
    ring.AddPoint(xRight, yTop)
    ring.AddPoint(xRight, yBottom)
    ring.AddPoint(xLeft, yTop)
    rasterGeometry = ogr.Geometry(ogr.wkbPolygon)
    rasterGeometry.AddGeometry(ring)

    return rasterGeometry


def extract_dic(tab_classes):
    """
    Transform an array containing classes to a dictionnary
    @params:
        tab_classes   		  - Required  : array containing classes definition
    """
    dic_res = {}

    for e in tab_classes:
        _tmp = e.split(':')
        dic_res[int(_tmp[0])] = int(_tmp[1])

    return dic_res


def computing_roads(output_folder, shapefile_roads, roads_pix_value, output_merging_tif, output_name_with_roads):
    output_filename_reproj = os.path.join(output_folder, 'roads_reproj.shp')
    change_projection(output_merging_tif, shapefile_roads, 'VAL', None, output_filename_reproj, [99], None)

    output_filename_rasterize = os.path.join(output_folder, 'roads_reproj.tif')
    rasterize_ground_ref(output_merging_tif, output_filename_reproj, 'VAL', output_filename_rasterize)

    dataset_roads = gdal.Open(output_filename_rasterize, GA_ReadOnly)
    dataset_gt = gdal.Open(output_merging_tif, GA_ReadOnly)

    band_roads = dataset_roads.GetRasterBand(1).ReadAsArray()
    band_gt = dataset_gt.GetRasterBand(1).ReadAsArray()

    band_roads[band_roads == 255] = 0
    _tmp = band_gt + band_roads

    _tmp[_tmp >= 90] = roads_pix_value
    _tmp[_tmp == 25] = 5

    shape = _tmp.shape

    output = output_name_with_roads

    geo_transform = dataset_gt.GetGeoTransform()
    driver = gdal.GetDriverByName("GTiff")
    dst_ds = driver.Create(output, shape[1], shape[0], 1, gdal.GDT_UInt16)
    dst_ds.SetProjection(dataset_gt.GetProjection())
    dst_ds.SetGeoTransform(geo_transform)
    dst_ds.GetRasterBand(1).WriteArray(_tmp)

    return output


def check_overlap(path_raster_1, path_raster_2):
    """
    Checking overlap of two rasters
    @params:
        path_raster_1   		  - Required  : path to the first raster
        path_raster_2   		  - Required  : path to the second raster
    """
    raster_1 = gdal.Open(path_raster_1)
    raster_2 = ogr.Open(path_raster_2)

    raster_1_geometry = get_geometry(raster_1)
    raster_2_geometry = get_geometry(raster_2)

    return raster_1_geometry.Intersect(raster_2_geometry)


def merge_two_dicts(x, y):
    """
    Merging two dictionnaries
    @params:
        x   		  - Required  : first dictionnary
        z   		  - Required  : second dictionnary
    """
    z = x.copy()
    z.update(y)
    return z


def extract_info_patches(patch_name, is_gr=False):
    """
    Extract informations from patch name.
    @params:
        patch_name   	- Required  : filename of the patch
        is_gr           - Required  : boolean to tell if the filename is ground reference or not
    """
    split = patch_name.split('_')

    if is_gr:
        tile = split[0]
        tile_x = split[2]
        tile_y = split[3]
        output = (tile, tile_x, tile_y)
    else:
        tile = split[0]
        date_img = split[1]
        sat = split[2]
        tile_x = split[3]
        tile_y = split[4]
        output = (tile, date_img, sat, tile_x, tile_y)

    return output


def get_corresponding_sat(path_folder, tile_x, tile_y):
    """
    Extract every file in the folder corresponding to the coordinates and concatenate them in one single string.
    @params:
        path_folder   	- Required  : filename of the patch
        tile_x          - Required  : boolean to tell if the filename is ground reference or not
        tile_y          - Required  : boolean to tell if the filename is ground reference or not
    """
    output = None

    for file in os.listdir(path_folder):
        string = tile_x + '_' + tile_y
        if string in file:
            if output is None:
                output = file
            else:
                output += ';' + file

    return output


def get_projection_patch(path_patch):
    """
    Get the projection (WKT) of a patch.
    @params:
        path_patch   	- Required  : path to the patch
    """
    dataset = gdal.Open(path_patch, GA_ReadOnly)
    proj = dataset.GetProjection()

    return proj


def list_to_string(list, sep=';'):
    """
    Convert a list to string using a separator between each element.
    @params:
        list   	- Required  : list to convert
        sep     - Required  : separator to put between each element
    """
    res = None

    for e in list:
        if res is None:
            res = str(e)
        else:
            res += sep + str(e)

    return res


def extract_classes_gr(path_gr_patch):
    """
    Extract classes in ground reference patch.
    @params:
        path_gr_patch   	- Required  : filename of the patch
    """
    dataset = gdal.Open(path_gr_patch, GA_ReadOnly)
    array = dataset.GetRasterBand(1).ReadAsArray()

    unique, counts = np.unique(array, return_counts=True)

    return dict(zip(unique, counts))


def create_json_labels(path_dataset):
    """
    Construct JSON file to generate labels for scene classification.
    @params:
        path_ground_reference   	- Required  : path to ground reference folder
        path_output_json            - Required  : path to the json directory
    """
    for file in os.listdir(os.path.join(path_dataset, 'ground_reference')):
        data = {}
        patch_infos = extract_info_patches(os.path.basename(file), is_gr=True)
        data['corresponding_s2'] = get_corresponding_sat(os.path.join(path_dataset, 's2'), patch_infos[1],
                                                         patch_infos[2])
        data['corresponding_s1'] = get_corresponding_sat(os.path.join(path_dataset, 's1'), patch_infos[1],
                                                         patch_infos[2])
        data['projection'] = get_projection_patch(os.path.join(path_dataset, 'ground_reference', file))
        data['labels'] = list_to_string(extract_classes_gr(os.path.join(path_dataset, 'ground_reference', file)).keys())

        json_name = patch_infos[0] + '_' + patch_infos[1] + '_' + patch_infos[2] + '.json'

        if os.path.exists(os.path.join(path_dataset, 'labels', json_name)):
            os.remove(os.path.join(path_dataset, 'labels', json_name))

        with open(os.path.join(path_dataset, 'labels', json_name), 'w') as file:
            json.dump(data, file)


def calculate_stats_dataset(path_dataset):

    return 0


def main():
    cwd = os.getcwd()
    # TODO: Regroupement des classes
    urban_regroup_classes = '11111:1 11112:1 11113:1 11121:1 11122:1 11123:1 11211:2 11212:2 11213:2 11221:2 11222:2 11223:2 11231:2 11232:2 11233:2 11241:2 11242:2 11243:2 11301:2 11302:2 11303:2 11401:3 11402:3 11403:3 12111:3 12112:3 12113:3 12121:3 12122:3 12123:3 12131:3 12132:3 12133:4 13141:3 13142:3 13123:4 12151:3 12152:3 12153:4 12201:3 12202:3 12203:4 13111:3 13112:3 13113:3 13121:3 13122:3 13123:3 13131:3 13132:3 13133:3 13141:3 13142:3 13143:4 13201:1 13202:3 13203:3 13301:2 13302:2 13303:2 13401:3 13402:3 13403:3 14111:3 14112:3 14113:5 14201:3 14202:3 14203:4 14301:3 14302:3 14303:3 15101:3 15102:3 15103:4 16101:3 16102:3 16103:3 17101:1 17102:1 17103:1 12141:4 12142:4 12143:3 14131:5 14132:5 14133:5 14121:-99 14122:-99 14123:-99'
    natural_regroup_classes = '2110:6 2120:6 2210:7 2221:8 2222:8 2223:8 2310:9 2320:10 3110:11 3120:11 3130:11 3140:11 3150:11 3210:12 3220:12 3230:12 3310:13 3320:13 3340:13 4110:14 4120:14 5110:15 5120:15 5130:15'
    path_shapefile_sentinel = './inputs/S2_Tiling_GE_test.shp'
    field = 'Name'
    directory_images_s2 = os.path.join(cwd, 'S2_img')
    directory_images_s1 = os.path.join(cwd, 'S1_img')
    shapefile_roads = './inputs/routes_finales_GE.shp'
    patch_size = 256
    directory_ground_reference_raster = os.path.join(cwd, 'ground_reference_raster')
    directory_tmp = os.path.join(cwd, 'tmp')
    directory_dataset = os.path.join(cwd, 'dataset')
    path_shp_ground_reference = './inputs/merge_grand_est.shp'
    attribute_shp_ground_reference_urban = 'cod_n5'
    attribute_shp_ground_reference_natural = 'cod_n4'

    # Read all S2-Tiles and put them in a, array
    # tiles_list = get_sentinel_tiles(path_shapefile_sentinel, field)
    tiles_list = ['31UGP', '32ULU', '32ULV', '32UMU', '32UMV']
    #tiles_list = ['32ULU']

    # Download Sentinel for each tile
    start_date = '2020-02-01'
    #end_date = '2020-02-10'
    end_date = '2020-09-30'

    ''''print(datetime.now().isoformat() + ' Downloading S2 images from ' + str(
        start_date) + ' to ' + str(end_date) + ' for each S2 tiles.')
    for tile in tiles_list:
        print('')
        download_sentinel2(tile, start_date, end_date)
    print(datetime.now().isoformat() + ' Successfully downloaded S2 images from ' + str(
        start_date) + ' to ' + str(end_date) + '.')

    print(datetime.now().isoformat() + ' Unzipping S2 images downloaded.')
    for file in os.listdir(cwd):
        file_path = os.path.join(cwd, file)
        if zipfile.is_zipfile(file_path):
            unzip(file_path, directory_images_s2)
    print(datetime.now().isoformat() + ' Successfully unzipped S2 images.')

    print(datetime.now().isoformat() + ' Rasterizing ground reference for each S2 tiles.')
    if not os.path.exists(directory_ground_reference_raster):
        os.mkdir(directory_ground_reference_raster)

    if not os.path.exists(directory_tmp):
        os.mkdir(directory_tmp)

    dic_classes_urban = extract_dic(urban_regroup_classes.split(' '))
    dic_classes_natural = extract_dic(natural_regroup_classes.split(' '))

    attributes_merge = []
    for key in dic_classes_urban:
        if dic_classes_urban[key] not in attributes_merge:
            attributes_merge.append(dic_classes_urban[key])

    if -99 in attributes_merge:
        attributes_merge.remove(-99)

    all_classes_urban = []
    for key in dic_classes_urban:
        if key not in all_classes_urban:
            all_classes_urban.append(key)

    all_classes_natural = []
    for key in dic_classes_natural:
        if key not in all_classes_natural:
            all_classes_natural.append(key)

    for tile in tiles_list:
        for img in os.listdir(directory_images_s2):
            if tile in img:
                path_band_ref = os.path.join(directory_images_s2, img, img + '_SRE_B3.tif')
                output_filename_tmp = os.path.join(cwd, directory_tmp,
                                                     tile + '_tmp.tif')
                shapefile_to_raster_tmp = os.path.join(cwd, directory_tmp,
                                                         path_shp_ground_reference.split('/')[2][:-4] + '_tmp.shp')
                output_filename_merge_tmp = os.path.join(cwd, directory_tmp,
                                                     tile + '_merge_tmp.tif')
                output_filename_merge_final = os.path.join(cwd, directory_ground_reference_raster,
                                                     tile + '_ground_reference.tif')

                # Changing ground refenre projection and reclassify
                change_projection(path_band_ref, path_shp_ground_reference, attribute_shp_ground_reference_urban,
                                  attribute_shp_ground_reference_natural, shapefile_to_raster_tmp, all_classes_urban, all_classes_natural)

                # Rasterizing shapefile
                rasterize_ground_ref(path_band_ref, shapefile_to_raster_tmp, attribute_shp_ground_reference_urban,
                                     output_filename_tmp)

                merge_classes(output_filename_tmp, merge_two_dicts(dic_classes_urban, dic_classes_natural), output_filename_merge_tmp)

                #Adding roads on reference data
                pixel_roads_value = 25
                computing_roads(os.path.join(cwd, directory_tmp), shapefile_roads, pixel_roads_value,
                                output_filename_merge_tmp, output_filename_merge_final)

                break

    print(datetime.now().isoformat() + ' Successfully rasterized ground reference for each S2 tiles.')

    print(datetime.now().isoformat() + ' Downloading and calibrating Sentinel-1 data.')
    download_sentinel1(tiles_list, directory_images, step=8)
    print(datetime.now().isoformat() + ' Successfully downloaded and calibrated Sentinel-1 data.')

    print(datetime.now().isoformat() + ' Checking Sentinel-1 data.')
    s1_output = './S1_img'
    json = check_satellite_data(tiles_list, directory_images, s1_output)
    print(datetime.now().isoformat() + ' Successfully checked Sentinel-1 data.')'''

    print(datetime.now().isoformat() + ' Calculating patches for each tile.')
    if not os.path.exists(directory_dataset):
        os.mkdir(directory_dataset)

    #create_folders_dataset(directory_dataset)

    for tile in tiles_list:
        '''calculate_patches_tile(tile, directory_ground_reference_raster, directory_images_s2,
                               directory_images_s1,
                               directory_dataset, patch_size, method=1)'''
        create_json_labels(directory_dataset)
    print(datetime.now().isoformat() + ' Successfully calculated patches for each tile.')

    print(datetime.now().isoformat() + ' Checking overlapped patches.')
    # TODO : checker l'overlap entre les tuiles adjacentes (entre tous les patchs S2 puis supprimer les patchs S1
    #  correspondants)
    print(datetime.now().isoformat() + ' Successfully checked and deleted overlapped patches')


if __name__ == '__main__':
    main()
