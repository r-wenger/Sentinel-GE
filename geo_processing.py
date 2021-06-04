import osgeo.gdal as gdal
import osgeo.ogr as ogr
import osgeo.osr as osr
from osgeo.gdalconst import *

from os.path import join, splitext
import os, shutil, sys

# Numpy Import
import numpy as np
import numpy.ma as ma
from scipy.ndimage import zoom

# Others imports
from skimage import morphology
import random
import time
import datetime


def rasterize_ground_ref(path_raster_reference, path_shapefile_to_raster, attribute_shapefile, output_filename):
    """
    Rasterizing ground truth.
    @params:
        path_raster_reference   	- Required  : raster or satellite image as reference
        path_shapefile_to_raster   	- Required  : shapefile to rasterize
        attribute_shapefile   		- Required  : attribute to rasterize (name of shapefile column)
        output_filename   			- Required  : output file rasterized
    """

    def compute_areas(pathshapefile, attribute_code, attribute_area, output_stats):
        driver = ogr.GetDriverByName("ESRI Shapefile")
        dataSource = driver.Open(pathshapefile, 0)
        layer = dataSource.GetLayer()
        dic = {}

        for feature in layer:
            _local_area = feature.GetField(str(attribute_area))
            _local_code = feature.GetField(str(attribute_code))

            if str(_local_code) in dic:
                dic[str(_local_code)] += float(_local_area)
            else:
                dic[str(_local_code)] = float(_local_area)

            feature = None

        dataSource = None

        with open(output_stats, "w+") as f:
            for key, value in dic.items():
                f.write(str(key) + ',' + str(float(value) / 1000000) + '\n')

    _name_stats_csv = str(splitext(path_shapefile_to_raster)[0] + '_stats_areas.csv')
    compute_areas(path_shapefile_to_raster, attribute_shapefile, 'AREA', _name_stats_csv)

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


def change_projection(path_raster_reference, path_shapefile_to_raster, attribute_shapefile, output_filename,
                      attribute_values=[]):
    """
    Changing shapefile projection to match with raster.
    @params:
        path_raster_reference   	- Required  : raster or satellite image as reference
        path_shapefile_to_raster   	- Required  : shapefile to rasterize
        attribute_shapefile   		- Required  : attribute to rasterize (name of shapefile column)
        output_filename   			- Required  : output file rasterized
        attribute_values   			- Required  : list of attributes to rasterize
    """
    if len(attribute_values) == 0:
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
    outlayer.CreateField(ogr.FieldDefn('AREA', ogr.OFTReal))
    outlayer.CreateField(ogr.FieldDefn(attribute_shapefile, ogr.OFTInteger))

    i = 0

    for feature in layer:
        _val = None

        if feature.GetField(attribute_shapefile) is not None:
            if int(feature.GetField(attribute_shapefile)) in attribute_values:
                _val = feature.GetField(attribute_shapefile)
            else:
                _val = 0
            transformed = feature.GetGeometryRef()
            transformed.Transform(transform)

            geom = ogr.CreateGeometryFromWkb(transformed.ExportToWkb())
            defn = outlayer.GetLayerDefn()
            feat = ogr.Feature(defn)
            feat.SetField('ID', i)
            feat.SetField('AREA', geom.GetArea())
            # TODO : modifications pour suppression emprises scolaires
            if (int(geom.GetArea()) < 300) and (int(str(feature.GetField(attribute_shapefile))[0:-1]) == 1211):
                feat.SetField(attribute_shapefile, '99')
                feat.SetGeometry(geom)
                outlayer.CreateFeature(feat)
            else:
                feat.SetField(attribute_shapefile, _val)
                feat.SetGeometry(geom)
                outlayer.CreateFeature(feat)
            i += 1
            feat = None

    ds = None


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


def clip_raster_wt_shapefile(roi_shape, field, gt_raster, img_raster, output_path):
    def ClipRasterWithPolygon(rasterPath, polyPath, outPath):
        os.system(
            "gdalwarp -q -cutline " + polyPath + " -crop_to_cutline " + " -of GTiff " + rasterPath + " " + outPath)

    def CreateClippingPolygons(inPath, field):
        driverSHP = ogr.GetDriverByName("ESRI Shapefile")
        ds = driverSHP.Open(inPath)

        if ds is None:
            print("Layer not open")
            sys.exit(1)

        lyr = ds.GetLayer()
        spatialRef = lyr.GetSpatialRef()

        i = 0

        for feature in lyr:
            fieldVal = feature.GetField(field)

            path_output = os.path.join(output_path, str(fieldVal) + '_' + str(i) + '.shp')

            if not os.path.exists(path_output):
                outds = driverSHP.CreateDataSource(path_output)
                outlyr = outds.CreateLayer(path_output, srs=spatialRef, geom_type=ogr.wkbPolygon)
                outDfn = outlyr.GetLayerDefn()
                ingeom = feature.GetGeometryRef()
                outFeat = ogr.Feature(outDfn)
                outFeat.SetGeometry(ingeom)
                outlyr.CreateFeature(outFeat)
                i += 1

    CreateClippingPolygons(roi_shape, field)

    for f in os.listdir(output_path):
        if '.shp' in f:
            out_path_gt = os.path.join(output_path, os.path.splitext(f)[0]) + '_gt.tif'
            out_path_img = os.path.join(output_path, os.path.splitext(f)[0]) + '_img.tif'
            ClipRasterWithPolygon(gt_raster, os.path.join(output_path, f), out_path_gt)
            ClipRasterWithPolygon(img_raster, os.path.join(output_path, f), out_path_img)


def computing_roads(output_folder, shapefile_roads, roads_pix_value, output_merging_tif):
    output_filename_reproj = os.path.join(output_folder, 'roads_reproj.shp')
    change_projection(output_merging_tif, shapefile_roads, 'VAL', output_filename_reproj, [99])

    output_filename_rasterize = os.path.join(output_folder, 'roads_reproj.tif')
    rasterize_ground_ref(output_merging_tif, output_filename_reproj, 'VAL', output_filename_rasterize)

    dataset_roads = gdal.Open(output_filename_rasterize, GA_ReadOnly)
    dataset_gt = gdal.Open(output_merging_tif, GA_ReadOnly)

    band_roads = dataset_roads.GetRasterBand(1).ReadAsArray()
    band_gt = dataset_gt.GetRasterBand(1).ReadAsArray()

    band_roads[band_roads == 255] = 0
    _tmp = band_gt + band_roads

    _tmp[_tmp >= 90] = roads_pix_value

    shape = _tmp.shape

    output = os.path.splitext(output_merging_tif)[0] + '_wROADS.tif'

    geo_transform = dataset_gt.GetGeoTransform()
    driver = gdal.GetDriverByName("GTiff")
    dst_ds = driver.Create(output, shape[1], shape[0], 1, gdal.GDT_UInt16)
    dst_ds.SetProjection(dataset_gt.GetProjection())
    dst_ds.SetGeoTransform(geo_transform)
    dst_ds.GetRasterBand(1).WriteArray(_tmp)

    return output


def random_selection_image(path_sat_image, path_ground_reference, output_ref, output_image, learning_method, attributes,
                           sat_lvl, bands_to_compute=['B2', 'B3', 'B4', 'B8'], indexes_to_compute=['NDVI'],
                           convert_to_tif=True):
    """
    Random selection of small images.
    @params:
        path_sat_image   		- Required  : path to the folder containing raster bands
        path_ground_reference   - Required  : path to the raster corresponding to the reference
        output_ref      		- Required  : outpupath to write the reference subset image
        output_image      		- Required  : outpupath to write the satellite subset image containing n bands
        bands_to_compute    	- Required  : bands used to compute. Default : B02, B03, B04, B08
        img_size        		- Required  : size of the subsets. Default 256
        convert_to_tif    		- Required  : if it is needed to convert bands to tif
        thresold_attributes    	- Required  : threshold used to pick random images
    """
    dic_path_tif = {}
    dic_path_all_ands = {}
    path_delete = None

    if sat_lvl == 'SPOT':
        path_delete = os.path.dirname(os.path.abspath(path_sat_image))
    else:
        path_delete = path_sat_image

    delete_files_with_keyword(path_delete, 'index')

    if sat_lvl == 'SPOT':
        dataset_spot = gdal.Open(path_sat_image, GA_ReadOnly)

        folder_spot = os.path.dirname(os.path.abspath(path_sat_image))
        for _band in bands_to_compute:
            band_spot = dataset_spot.GetRasterBand(int(_band)).ReadAsArray()

            shape = band_spot.shape
            outputpathband = os.path.join(folder_spot, 'SPOT_BAND' + str(_band) + '.tif')

            geo_transform = dataset_spot.GetGeoTransform()
            driver = gdal.GetDriverByName("GTiff")
            dst_ds = driver.Create(outputpathband, shape[1], shape[0], 1, gdal.GDT_UInt16)
            dst_ds.SetProjection(dataset_spot.GetProjection())
            dst_ds.SetGeoTransform(geo_transform)
            dst_ds.GetRasterBand(1).WriteArray(band_spot)
            print(_band)
            dic_path_tif[_band] = outputpathband
            dic_path_all_ands[_band] = outputpathband
    else:
        for f in os.listdir(path_sat_image):
            if os.path.isfile(join(path_sat_image, f)):
                _band = str(splitext(f)[0].split('_')[7])
                if (_band in bands_to_compute) and (get_sat_lvl_ref(sat_lvl) in f):
                    _name_local_tif = str(join(path_sat_image, splitext(f)[0] + '.tif'))
                    if convert_to_tif and (not os.path.exists(_name_local_tif)):
                        _command = 'gdal_translate -of GTiff ' + str(join(path_sat_image, f)) + ' ' + _name_local_tif
                    # os.system(_command)
                    dic_path_tif[_band] = _name_local_tif
                    dic_path_all_ands[_band] = _name_local_tif
                else:
                    if (get_sat_lvl_ref(sat_lvl) in f):
                        _tmp_name = str(join(path_sat_image, splitext(f)[0] + '.tif'))
                        dic_path_all_ands[_band] = _tmp_name

    if indexes_to_compute is not None:
        for e in indexes_to_compute:
            if ('NDVI' in e):
                output_name = join(path_sat_image, 'NDVI_index.tif')
                compute_difference(dic_path_all_ands['B8'], dic_path_all_ands['B4'], dic_path_all_ands['B2'],
                                   output_name)
                dic_path_tif['NDVI'] = output_name
            elif ('NDBI' in e):
                output_name = join(path_sat_image, 'NDBI_index.tif')
                compute_difference(dic_path_all_ands['B11'], dic_path_all_ands['B8'], dic_path_all_ands['B2'],
                                   output_name)
                dic_path_tif['NDBI'] = output_name

    if indexes_to_compute is not None:
        for e in indexes_to_compute:
            if ('HARALICK' in e):
                print('HARALICK')
                output_name = join(path_sat_image, 'HARALICK_index.tif')
                _command = 'otbcli_HaralickTextureExtraction -in ' + str(dic_path_tif[
                                                                             'NDVI']) + ' -step 3 -parameters.xrad 5 -parameters.yrad 5 -parameters.min -1 -parameters.max 1 -texture simple -out ' + str(
                    output_name)
                os.system(_command)
                output_entropy = join(path_sat_image, 'haralick_endvi_index.tif')

                _ds_haralick = gdal.Open(output_name, GA_ReadOnly)
                _tab = _ds_haralick.GetRasterBand(2).ReadAsArray()
                _geo_transform = _ds_haralick.GetGeoTransform()

                driver = gdal.GetDriverByName("GTiff")
                dst_ds_endvi = driver.Create(output_entropy, _ds_haralick.RasterXSize, _ds_haralick.RasterYSize, 1,
                                             gdal.GDT_Float32)
                dst_ds_endvi.SetProjection(_ds_haralick.GetProjection())
                dst_ds_endvi.SetGeoTransform(_geo_transform)
                dst_ds_endvi.GetRasterBand(1).WriteArray(_tab)

                dic_path_tif['HARALICK'] = output_name

    dataset_raster_ref = gdal.Open(path_ground_reference, GA_ReadOnly)
    dataset_raster_img = None

    if sat_lvl == 'SPOT':
        dataset_raster_img = gdal.Open(dic_path_all_ands['2'], GA_ReadOnly)
    else:
        dataset_raster_img = gdal.Open(dic_path_all_ands['B8'], GA_ReadOnly)

    if (dataset_raster_ref.RasterXSize == dataset_raster_img.RasterXSize) and (
            dataset_raster_ref.RasterYSize == dataset_raster_img.RasterYSize):

        tab_ref = dataset_raster_ref.GetRasterBand(1).ReadAsArray()
        geo_transform = dataset_raster_img.GetGeoTransform()
        img_size_x = int(dataset_raster_img.RasterXSize)
        img_size_y = int(dataset_raster_img.RasterYSize)

        _tmp_bands = []

        for key in dic_path_tif:
            _dataset_raster_img_tmp = gdal.Open(dic_path_tif[key], GA_ReadOnly)
            _tab_raster_tmp = resize_bands(_dataset_raster_img_tmp)
            _tmp_bands.append(_tab_raster_tmp)
            _dataset_raster_img_tmp = None

        _nb_bands = len(_tmp_bands)

        driver = gdal.GetDriverByName("GTiff")
        dst_ds = driver.Create(output_image, img_size_x, img_size_y, _nb_bands, gdal.GDT_Float32)
        dst_ds.SetProjection(dataset_raster_img.GetProjection())
        dst_ds.SetGeoTransform(geo_transform)

        for i in range(1, _nb_bands + 1):
            if learning_method == 'RF':
                dst_ds.GetRasterBand(i).WriteArray(_tmp_bands[i - 1])
            elif learning_method == 'NN':
                if sat_lvl == 'SPOT':
                    dst_ds.GetRasterBand(i).WriteArray(_tmp_bands[i - 1] / 255)
                else:
                    dst_ds.GetRasterBand(i).WriteArray(normalize(_tmp_bands[i - 1]))

        dst_ds = None

        dst_ds_ref = driver.Create(output_ref, img_size_x, img_size_y, len(attributes) + 1, gdal.GDT_UInt16)
        # dst_ds_ref = driver.Create(output_ref, img_size, img_size, 1, gdal.GDT_UInt16)
        dst_ds_ref.SetProjection(dataset_raster_ref.GetProjection())
        dst_ds_ref.SetGeoTransform(geo_transform)

        for i in range(1, len(attributes) + 1):
            _tmp_band_gf = np.where(tab_ref == int(attributes[i - 1]), 1, 0)
            dst_ds_ref.GetRasterBand(i).WriteArray(_tmp_band_gf)

        _tmp_band_gf = np.where(tab_ref == 0, 1, 0)
        dst_ds_ref.GetRasterBand(len(attributes) + 1).WriteArray(_tmp_band_gf)
        dst_ds_ref = None


def get_sat_lvl_ref(sat_lvl):
    res = None

    if str(sat_lvl) == 'L2A':
        res = 'SRE'
    elif str(sat_lvl) == 'L3A':
        res = 'FRC'
    else:
        raise Exception('Sat lvl unknown')

    return res


def normalize(band):
    """
    Normalize all datas in a numpy array
    @params:
        band   - Required  : a band in a numpy format
    """
    mean = np.mean(band)
    std = np.std(band)
    x = (band - mean) / std
    return x


def merge_classes(input_raster, apply_morphology, dic_classes, output_filename):
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

    if apply_morphology:
        se = morphology.disk(3)
        tab_return = morphology.closing(tab_return, se)

    if -99 in tab_return:
        tab_return = ma.masked_array(tab_return, tab_return == -99)

        for shift in (-3, 3):
            for axis in (0, 1):
                a_shifted = np.roll(tab_return, shift=shift, axis=axis)
                idx = ~a_shifted.mask * tab_return.mask
                tab_return[idx] = a_shifted[idx]

    driver = gdal.GetDriverByName("GTiff")
    dst_ds = driver.Create(output_filename, ds_raster.RasterXSize, ds_raster.RasterYSize, 1, gdal.GDT_Byte)
    dst_ds.SetProjection(ds_raster.GetProjection())
    geotransform = ds_raster.GetGeoTransform()
    dst_ds.SetGeoTransform(geotransform)
    dst_ds.GetRasterBand(1).WriteArray(tab_return)
    dst_ds = None


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


def remove_all_content_in_folder(folder):
    for filename in os.listdir(folder):
        file_path = os.path.join(folder, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except Exception as e:
            print('Failed to delete %s. Reason: %s' % (file_path, e))


def check_files_directories_integrity(folder_raster_ref, path_shapefile_to_raster_ua, path_directory):
    if os.path.exists(path_shapefile_to_raster_ua):
        print('%s: Ground reference shapefile is ok:' % (datetime.datetime.now().isoformat()))
        print(path_shapefile_to_raster_ua)
    else:
        print('%s: Ground reference shapefile is NOT ok:' % (datetime.datetime.now().isoformat()))
        print(path_shapefile_to_raster_ua)
        sys.exit(1)

    if os.path.exists(folder_raster_ref):
        print('%s: Sattelite image is ok:' % (datetime.datetime.now().isoformat()))
        print(folder_raster_ref)
    else:
        print('%s: Sattelite image is ok is NOT ok:' % (datetime.datetime.now().isoformat()))
        print(folder_raster_ref)
        sys.exit(1)

    if os.path.exists(path_directory):
        remove_all_content_in_folder(path_directory)
        os.makedirs(os.path.join(path_directory, 'mask'))
        os.makedirs(os.path.join(path_directory, 'img'))
        print('%s: Directories created:' % (datetime.datetime.now().isoformat()))
        print(os.path.join(path_directory, 'mask') + ' - ' + os.path.join(path_directory, 'img'))
    else:
        print('%s: Directories NOT created:' % (datetime.datetime.now().isoformat()))
        print(os.path.join(path_directory, 'mask') + ' - ' + os.path.join(path_directory, 'img'))
        raise Exception('Error')


def extract_dic(tab_classes):
    dic_res = {}

    for e in tab_classes:
        _tmp = e.split(':')
        dic_res[int(_tmp[0])] = int(_tmp[1])

    return dic_res


def create_folders(output, step_to_execute):
    mask_img_path = os.path.join(output, 'mask_img')
    mask_path = os.path.join(mask_img_path, 'mask')
    img_path = os.path.join(mask_img_path, 'img')
    output_path_roi = os.path.join(output, 'shapes_cut')
    output_path_patches = os.path.join(output, 'patches')

    if ('all' in step_to_execute) or ('pretreatment' in step_to_execute):
        if os.path.exists(mask_img_path):
            remove_all_content_in_folder(output)

        os.makedirs(output_path_patches)
        os.makedirs(mask_img_path)
        os.makedirs(mask_path)
        os.makedirs(img_path)
        os.makedirs(output_path_roi)

    return mask_path, img_path, output_path_roi


def get_raster_reference(folder_sat, sat_lvl):
    path_res = None

    if sat_lvl == 'L2A':
        for f in os.listdir(folder_sat):
            if ('SRE' in f) and ('B2' in f):
                path_res = os.path.join(folder_sat, f)
                break
    elif sat_lvl == 'L3A':
        for f in os.listdir(folder_sat):
            if ('FRC' in f) and ('B2' in f):
                path_res = os.path.join(folder_sat, f)
                break
    else:
        print('Sat lvl should be L2A or L3A.\n')
        sys.exit(1)

    return path_res


def delete_files_with_keyword(path_sat, keyword):
    for f in os.listdir(path_sat):
        if keyword in f:
            os.remove(join(path_sat, f))


def compute_difference(band1, band2, band_ref, output):
    dataset_band1 = gdal.Open(band1, GA_ReadOnly)
    dataset_band2 = gdal.Open(band2, GA_ReadOnly)
    dataset_ref = gdal.Open(band_ref, GA_ReadOnly)

    tab_band1 = resize_bands(dataset_band1)
    tab_band2 = resize_bands(dataset_band2)

    numerator = tab_band1 - tab_band2
    denominator = tab_band1 + tab_band2

    denominator[denominator == 0] = 1

    res = np.divide(numerator, denominator)

    shape = res.shape

    geo_transform = dataset_ref.GetGeoTransform()
    driver = gdal.GetDriverByName("GTiff")
    dst_ds = driver.Create(output, shape[0], shape[1], 1, gdal.GDT_Float32)
    dst_ds.SetProjection(dataset_band1.GetProjection())
    dst_ds.SetGeoTransform(geo_transform)
    dst_ds.GetRasterBand(1).WriteArray(res)


def generate_datas(path_sat_image, path_ground_reference, path_subset_gf, path_subset_img, learning_method, attributes,
                   sat_lvl, bands_to_compute=['B02', 'B03', 'B04', 'B08'], indexes_to_compute=None,
                   convert_to_tif=True):
    """
    Call in a loop to create terminal progress bar
    @params:
        path_sat_image   		- Required  : current iteration (Int)
        path_ground_reference   - Required  : total iterations (Int)
        path_subset_gf      	- Optional  : prefix string (Str)
        path_subset_img      	- Optional  : suffix string (Str)

    """
    if (len(os.listdir(path_subset_gf)) != 0) and (len(os.listdir(path_subset_img)) != 0):
        remove_all_content_in_folder(path_subset_gf)
        remove_all_content_in_folder(path_subset_img)

    _tmp_subset_gf = join(path_subset_gf, 'gf.tif')
    _tmp_subset_img = join(path_subset_img, 'img.tif')

    random_selection_image(path_sat_image, path_ground_reference, _tmp_subset_gf, _tmp_subset_img, learning_method,
                           attributes, sat_lvl, bands_to_compute, indexes_to_compute, convert_to_tif)
