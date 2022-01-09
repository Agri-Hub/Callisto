import os
import re
import fnmatch
import sys
import csv
import time
import requests
import json
#import yaml
import socket
import psutil
from scipy.interpolate import interp1d
from glob import glob
import numpy as np
import logging
import xml.etree.ElementTree as ET
from osgeo import gdal, ogr, osr
from datetime import datetime, date, timedelta
from dateutil.relativedelta import relativedelta



def get_ip_address(def_apache_port='45200'):
    ''' Get the local IP address of the current machine'''

    s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    s.connect(("8.8.8.8", 80))
    return s.getsockname()[0]+':'+def_apache_port


class MyLogger(object):
    logging.basicConfig(format='%(asctime)s %(levelname)s \t %(message)s', level=logging.DEBUG)
    logger = logging.getLogger(__name__)


class Mapfile(MyLogger):
    ''' Create a mapfile for the publication of layers, used by UMN Mapserver'''

    def __init__(self, base_dir, extent):
        self.base_dir = base_dir
        self.extent = extent
        self.path_tif = self._get_path_tifs()
        self.path_vect = self._get_path_shp()

    def create(self):
        ''' Create a mapfile '''

        rasters = ''
        for ras_name, ras_path in self.path_tif.items():
            rasters += self._add_raster_layer(ras_name, ras_path)

        vectors = ''
        for vec_name, vec_path in self.path_vect.items():
            vectors += self._add_vector_layer(vec_name, vec_path)

        mapfile = self._add_map(rasters, vectors)

        with open(os.path.join(self.base_dir, 'wms.map'), 'w') as map:
            map.write(mapfile)

        self.logger.info('Mapfile just created')


    def _add_raster_layer(self, name, tif_path):
        ''' Add a raster layer block to the mapfile - returns string '''

        # OFFSITE 255 30 30
        # added by me
        # PROCESSING "AUTO"
        # STATUS DEFAULT
        layer_block = '''
            LAYER
             NAME {name}
             PROJECTION
              "init=epsg:3857"
             END
             TYPE RASTER
             DUMP TRUE
             STATUS DEFAULT
             DATA ./{tif_path}
             OFFSITE 70 74 66
             PROCESSING "SCALE"
             STATUS DEFAULT
               CLASS
                 NAME "red1"
                 EXPRESSION ([pixel] >=-1 AND [pixel] < -0.1)
                 STYLE
                   COLORRANGE 0 0 0  255 255 255
                   DATARANGE -1 -0.1
                 END
               END
               CLASS
                 NAME "red2"
                 EXPRESSION ([pixel] >= -0.1 AND [pixel] < 0.2)
                 STYLE
                   COLOR "#A5160D"
                 END
               END
               CLASS
                 NAME "red3"
                 EXPRESSION ([pixel] >= 0.2 AND [pixel] < 0.35)
                 STYLE
                   COLOR "#BA4F33"
                 END
               END
               CLASS
                 NAME "red4"
                 EXPRESSION ([pixel] >= 0.35 AND [pixel] < 0.45)
                 STYLE
                   COLOR "#D08859"
                 END
               END
               CLASS
                 NAME "red5"
                 EXPRESSION ([pixel] >= 0.45 AND [pixel] < 0.55)
                 STYLE
                   COLOR "#E6C17F"
                 END
               END
               CLASS
                 NAME "red6"
                 EXPRESSION ([pixel] >= 0.55 AND [pixel] < 0.65)
                 STYLE
                   COLOR "#FCFAA6"
                 END
               END
               CLASS
                 NAME "red7"
                 EXPRESSION ([pixel] >= 0.65 AND [pixel] < 0.75)
                 STYLE
                   COLOR "#D8E08E"
                 END
               END
               CLASS
                 NAME "red8"
                 EXPRESSION ([pixel] >= 0.75 AND [pixel] < 0.85)
                 STYLE
                   COLOR "#B4C777"
                 END
               END
               CLASS
                 NAME "red9"
                 EXPRESSION ([pixel] >= 0.85 AND [pixel] < 0.95)
                 STYLE
                   COLOR "#91AE60"
                 END
               END
                CLASS
                 NAME "red10"
                 EXPRESSION ([pixel] >= 0.95 AND [pixel] < 1.05)
                 STYLE
                   COLOR "#6D9548"
                 END
               END
                CLASS
                 NAME "red11"
                 EXPRESSION ([pixel] >= 1.05 AND [pixel] < 1.15)
                 STYLE
                   COLOR "#4A7C31"
                 END
               END
                CLASS
                 NAME "red12"
                 EXPRESSION ([pixel] >= 1.15 AND [pixel] < 1.25)
                 STYLE
                   COLOR "#26631A"
                 END
               END
                CLASS
                 NAME "red13"
                 EXPRESSION ([pixel] >= 1.25 AND [pixel] < 1.35)
                 STYLE
                   COLOR "#034A03"
                 END
               END
                 CLASS
                 NAME "red14"
                 EXPRESSION ([pixel] >= 1.35 AND [pixel] < 1.5)
                 STYLE
                   COLORRANGE "#055376" "#0D94B5"
                   DATARANGE 1.35 1.5
                 END
               END
             METADATA
              "ows_title"    "{name}"
              "ows_abstract"    "NOA, RECAP, Ortho"
              "ows_keywordlist"    "Sentinel 2"
              "wms_opaque"    "0"
              "wms_metadataurl_format" "text/html"
              "wms_dataurl_format" "image/png"
              "ows_metadataurl_type" "FGDC"
              "ows_metadataurl_href" "https://www.astro.noa.gr/"
              "wms_dataurl_href" "https://www.astro.noa.gr/"
             END
            END
           '''
        layer_data = {'name': name, 'tif_path': tif_path}

        self.logger.info('Layer %s just added at mapfile.' % name)

        return layer_block.format(**layer_data)


    def _add_vector_layer(self, name, vect_path):
        ''' Add a vector layer block to the mapfile - returns string '''

        layer_block = '''
         LAYER
          NAME {name}
          PROJECTION
           "init=epsg:3857"
          END
          TYPE POLYGON
          STATUS ON
          DATA ./{vect_path}
          CLASS
           NAME "Parcels of Classification."
           STYLE
            COLOR       232 232 232
            OUTLINECOLOR 32 32 32
           END
          END
          METADATA
           "ows_title"    "{name}"
           "ows_abstract"    "NOA, RECAP, Ortho"
           "ows_keywordlist"    "Sentinel 2"
           "wms_opaque"    "0"
           "wms_metadataurl_format" "text/html"
           "wms_dataurl_format" "image/png8"
           "ows_metadataurl_type" "FGDC"
           "ows_metadataurl_href" "https://www.astro.noa.gr/"
           "wms_dataurl_href" "https://www.astro.noa.gr/"
          END
         END
        '''
        layer_data = {'name': name, 'vect_path': vect_path}

        self.logger.info('Layer %s just added at mapfile.' %name)

        return layer_block.format(**layer_data)

    def _add_map(self, in_rasters, in_vectors):
        ''' Add a map block - returns string'''

        base_block = '''
        MAP
         NAME "WMS"
         IMAGECOLOR 10000 10000 10000
         SIZE 600 400
         STATUS ON

         PROJECTION
          "init=epsg:3857"
         END
         EXTENT {xmin} {ymin} {xmax} {ymax}

         OUTPUTFORMAT
          NAME png8
          DRIVER "GD/PNG8"
          MIMETYPE "image/png8"
          IMAGEMODE RGB
          TRANSPARENT ON
          EXTENSION "png8"
         END

         WEB
          IMAGEPATH "/var/www/ms_tmp/"
          IMAGEURL "/ms_tmp/"
          METADATA
           "ows_schemas_location" "http://schemas.opengis.net"
           "ows_title" "NOA RECAP WMS Service- Sentinel 2 ortho RGB images"
           "ows_abstract" "NOA RECAP WMS Service- Sentinel 2 ortho RGB images"
           "ows_keywordlist" "NOA,RECAP,MapServer,Ortho,RGB,Sentinel2"
           "ows_enable_request" "*"
           "ows_onlineresource" "{url}?map={path_mapfile}"
           "ows_fees" "none"
           "ows_accessconstraints" "none"
           "ows_srs" "EPSG:4326 EPSG:900913 EPSG:3857"
           "wms_feature_info_mime_type"  "text/html"
           "wms_getmap_formatlist" "image/png8,image/gif"
           END
         END

         QUERYMAP
          STATUS ON
          SIZE 200 200
          STYLE HILITE
          COLOR 255 255 0
         END

         LEGEND
          LABEL
           TYPE BITMAP
           SIZE MEDIUM
           COLOR 0 0 0
          END
         END

          {rasters_lyr}

          {vectors_lyr}

        END
        '''

        base_data = {'xmin': self.extent[0], 'ymin': self.extent[1], 'xmax': self.extent[2], 'ymax': self.extent[3], 'url': 'http://' + get_ip_address() +'/cgi-bin/mapserv', 'path_mapfile': os.path.join(self.base_dir, 'wms.map'), 'rasters_lyr': in_rasters, 'vectors_lyr': in_vectors}

        self.logger.info('Mapfile %s just added at mapfile.' %os.path.join(self.base_dir, 'wms.map'))

        return base_block.format(**base_data)

    def _get_path_tifs(self):
        ''' Get the relative paths of NDVI and RGB tif'''

        rasters = {}
        for root, dirs, files in os.walk(self.base_dir):
            for tif in fnmatch.filter(files, '*.tif'):
                if re.match("(.*)(NDVI|RGB|NDWI|SAVI)(.*)", tif):
                    rel_tif_path = os.path.join(os.path.split(root)[1], tif)
                    rasters[tif.split('.')[0]] = rel_tif_path

        self.logger.info('Raster\'s pathname of the directory %s have been successfully initialized', os.path.split(self.base_dir)[1])
        return rasters

    def _get_path_shp(self):
        ''' Get the relative paths of shapefiles'''

        vectors = {}
        for shp in glob(os.path.join(self.base_dir, '*.shp')):
            file_shp = os.path.split(shp)[1]
            vectors[file_shp.split('.')[0]]=file_shp

        self.logger.info('Vector\'s pathname of the directory %s have been successfully initialized', os.path.split(self.base_dir)[1])
        return vectors


class MapProxy(MyLogger):
    ''' Update (and overwrite) a mapproxy yaml file, with layers of an incoming WMS'''

    def __init__(self, wms_url, mapproxy_file, recap_fkey=None):
        self.wms_url = wms_url
        self.mapproxy_file = mapproxy_file
        self.recap_fkey = recap_fkey if recap_fkey else datetime.now().time().strftime('%H%M%S')
        #self.bbox = geom.extent

    def update_mapproxy_file(self):
        ''' Update MapProxy file (.yaml) '''

        with open(self.mapproxy_file, 'r') as map_proxy_file:
            map_proxy = yaml.load(map_proxy_file)
            MyLogger.logger.info("Successful open")
        layers = self._get_wms_layers()
        #layers = []
        #for item in layer_name:
            #layers.append((item,layer_bbox))
        MyLogger.logger.info(layers)
        for lyr_name, lyr_bbox in layers:
            MyLogger.logger.info(lyr_name)
            idx = '%s_%s' % (self.recap_fkey, lyr_name)
            src = self._add_src(self.wms_url, lyr_name, lyr_bbox)
            lyr = self._add_lyr(idx)
            map_proxy['sources'][idx] = src
            map_proxy['layers'].append(lyr)

        with open(self.mapproxy_file, 'w') as map_proxy_file:
            yaml.dump(map_proxy, map_proxy_file, default_flow_style=False)

        self.logger.info('MapProxy just completed.')

    def _get_wms_layers(self):
        ''' Returns a list of layer and bbox '''

        response = requests.get(self.wms_url + '&service=WMS&request=GetCapabilities&version=1.3.0')
        wms_xml = ET.fromstring(response.content)
        layer_name = [child.getchildren()[0].text for i, child in enumerate(wms_xml.iter('{http://www.opengis.net/wms}Layer')) if i > 0]
        logging.info(layer_name)
        layer_bbox = [ [float(child.getchildren()[0].text), float(child.getchildren()[2].text), float(child.getchildren()[1].text),float(child.getchildren()[3].text)] for i, child in enumerate(wms_xml.iter('{http://www.opengis.net/wms}EX_GeographicBoundingBox')) if i > 0]
        logging.info(layer_bbox)
        '''
        geoms = []
        for item in layer_name:
            x = layer_bbox[0]
            geoms.append(x)
        '''

        return zip(layer_name,layer_bbox)


    def _add_src(self, in_mapservice, in_layer, in_bbox):
        ''' Add a source to mapproxy file - dictionary '''

        src = {}
        src['type'] = 'wms'
        src['req'] = {
            'url': in_mapservice,
            'layers': in_layer,
            'srs': ['EPSG:4326', 'EPSG:3857', 'EPSG:900913']
            }
        src['coverage'] = {
            'bbox': in_bbox,
            'srs': 'EPSG:4326'
            }

        self.logger.info('WMS Source Layer (%s) has been successfully added at MapProxy file.' %in_layer )

        return src

    def _add_lyr(self, in_src):
        ''' Add a layer to mapproxy file - dictionary '''

        lyr = {}
        lyr['sources'] = [in_src]
        lyr['name'] = in_src
        lyr['title'] = in_src

        self.logger.info('Layer (%s) has been successfully added at MapProxy WMS service.' % in_src)

        return lyr


class Make_request(MyLogger):
    ''' Execute GET, POST, PUT, PATCH, DELETE requests'''

    def __init__(self, url, method='GET', data=None, headers=None):
        self.session = requests.Session()
        self.url = url
        self.method = method
        self.data = json.dumps(data)
        self.headers = {'Content-Type': 'application/json'}
        if headers:
            self.headers.update(headers)

    def execute_request(self):
        ''' Execute a request - returns dict {status, response}'''

        if self.method=='GET':
            response=self.session.get(url=self.url, headers=self.headers, verify=False)
            # response = requests.get(url=self.url, headers=self.headers, verify=False)

        elif self.method=='POST':
            response=self.session.post(url=self.url, headers=self.headers, data=self.data, verify=False)

        elif self.method=='PUT':
            response=self.session.put(url=self.url, headers=self.headers, data=self.data, verify=False)

        elif self.method=='PATCH':
            response=self.session.patch(url=self.url, headers=self.headers, data=self.data, verify=False)

        elif self.method=='DELETE':
            response=self.session.delete(url=self.url, headers=self.headers, data=self.data, verify=False)

        self.logger.info('URL: %s - STATUS: %d' %(self.url, response.status_code))
        return {'status': response.status_code , 'response': response.json()}

    @classmethod
    def create_token(cls, url, data):
        ''' Create a token - returns make_request '''

        token = cls(url=url, method='POST', data=data)

        return token.execute_request()

    def close_request(self):
        ''' Close the request '''

        self.session.close()


class Sentinel2_pre_process(MyLogger):
    ''' Pre-processing of Sentinel 2 products '''

    prefixes = [('B02_10m', 'blue'), ('B03_10m', 'green'), ('B04_10m', 'red'), ('B05_20m', 'vre1'), ('B06_20m', 'vre2'), ('B07_20m', 'vre3'), ('B08_10m', 'nir1'), ('B8A_20m', 'nir2'), ('B11_20m', 'swir1'), ('B12_20m', 'swir2'), ('SCL_20m', 'cloud'),('TCI_10m','tci')]


    def __init__(self, base_dir, product):
        self.base_dir = base_dir
        self.product = product
        self.ws = os.path.join(self.base_dir, self.product)
        self.bands = self._get_path_bands()
        logging.info('%s %s %s %s',base_dir,product,self.ws,self.bands)

    def converting_rasters(self):
        ''' batch resampling and convertion rasters for the specified bands - returns nothing '''
        # self.logger.info("Start to convert")
        # logging.info("Start to convert")
        for p, in_image, out_image, prj_image in self.bands.values():
            logging.info('%s',in_image)
            if '20m' in in_image:
                self.resample_raster(p, in_image, out_image)
            if '10m' in in_image:
                self.raster2geotiff(p, in_image, out_image)

        self.logger.info('Resampling and Transformation of 20_m and 10_m bands of %s, respectively, just finished.', self.product)

    def delete_jp2_files(self):
        ''' Delete jp2 files from SAFE directory - returns nothing '''

        for p, in_image, out_image, prj_image in self.bands.values():
            if os.path.exists(os.path.join(p, out_image)):
                os.remove(os.path.join(p, in_image))
                self.logger.info('Removing of %s.', in_image)

    @staticmethod
    def clip_raster(img_dir, in_ras_file, out_ras_file, **extent):
        ''' Clip a raster (one band) - returns nothing'''

        xmin = extent.get('xmin')
        ymin = extent.get('ymin')
        xmax = extent.get('xmax')
        ymax = extent.get('ymax')

        srs = osr.SpatialReference()
        srs.ImportFromEPSG(3857)

        in_ras = gdal.Open(os.path.join(img_dir, in_ras_file))
        # in_ras = gdal.Open(in_ras_file)
        gt = in_ras.GetGeoTransform()
        inv_gt = gdal.InvGeoTransform(gt)

        off_ulx, off_uly = map(int, gdal.ApplyGeoTransform(inv_gt, xmin, ymax))
        off_lrx, off_lry = map(int, gdal.ApplyGeoTransform(inv_gt, xmax, ymin))

        # changed for korea
        #rows, columns = (off_lry - off_uly) + 1, (off_lrx - off_ulx) + 1
        rows, columns = (off_lry - off_uly), (off_lrx - off_ulx)
        in_band = in_ras.GetRasterBand(1)

        driver = gdal.GetDriverByName('GTiff')
        out_ds = driver.Create(os.path.join(img_dir, out_ras_file), columns, rows, 1, in_band.DataType)
        out_ds.SetProjection(in_ras.GetProjection())
        ulx, uly = gdal.ApplyGeoTransform(gt, off_ulx, off_uly)
        out_gt = list(gt)
        out_gt[0], out_gt[3] = ulx, uly
        out_ds.SetGeoTransform(out_gt)

        out_ds.GetRasterBand(1).WriteArray(in_band.ReadAsArray(off_ulx, off_uly, columns, rows))

        del in_ras
        del in_band
        del out_ds

        Sentinel2_pre_process.logger.info('Clipping of %s has been successfully finished', out_ras_file)

    @staticmethod
    def reproject_raster(img_dir, in_ras_file, out_ras_file):
        ''' Reproject a raster to EPSG:3857 - returns nothing'''

        # TODO, the output raster to have the same resolution as the input
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(3857)

        in_ras = gdal.Open(os.path.join(img_dir, in_ras_file))
        projected_vrt = gdal.AutoCreateWarpedVRT(in_ras, None, srs.ExportToWkt(), gdal.GRA_Bilinear)
        out_ds = gdal.GetDriverByName('GTiff').CreateCopy(os.path.join(img_dir, out_ras_file),  projected_vrt)

        Sentinel2_pre_process.logger.info('Reprojection of %s has been successfully finished', in_ras_file)

        del in_ras
        del projected_vrt
        del out_ds

    @staticmethod
    def raster2geotiff(img_dir, in_ras_file, out_ras_file):
        '''Converts a raster file to GeoTiff - returns nothing'''

        # TODO, you can also use gdal.Wrap

        # in_ras = gdal.Open(os.path.join(img_dir, in_ras_file))

        x = os.path.join(img_dir,in_ras_file)
        y = os.path.join(img_dir,out_ras_file)

        # os.system("gdal_translate -of GTiff " + x + " " + y)

        if '_TCI' in in_ras_file:
            os.system("gdal_translate -of GTiff "+x+" "+y)
        else:
            in_ras = os.path.join(img_dir, in_ras_file)
            os.system("gdal_translate -of GTiff " + x + " " + y)
            # in_band = in_ras.GetRasterBand(1)
            # driver = gdal.GetDriverByName('GTiff')
            # out_ds = driver.Create(os.path.join(img_dir, out_ras_file), in_ras.RasterXSize, in_ras.RasterYSize, 1, in_band.DataType)
            # out_ds.SetProjection(in_ras.GetProjection())
            # out_ds.SetGeoTransform(in_ras.GetGeoTransform())
            # out_band = out_ds.GetRasterBand(1)
            # out_band.WriteArray(in_band.ReadAsArray())
            # out_band.FlushCache()
            # out_band.ComputeStatistics(False)
            #
            # del out_band
            # del out_ds
            # del in_ras
            # os.system('gdalwarp -t_srs EPSG:3857 -dstnodata -9999.0 -tr 10.0 10.0 -r near -of GTiff %s %s' %(in_ras,os.path.join(img_dir, out_ras_file)))
            # gdal.Warp(in_ras,os.path.join(img_dir, out_ras_file),dstSRS='EPSG:3857', xRes=10,yRes=10)
            # gdal.Warp(os.path.join(img_dir, out_ras_file),in_ras,dstSRS='EPSG:3857', xRes=10,yRes=10)

        Sentinel2_pre_process.logger.info('Conversion of %s has been successfully finished', out_ras_file)

    @staticmethod
    def resample_raster(img_dir, in_ras_file, out_ras_file):
        '''Resample a 20m raster file to 10m GeoTiff - returns nothing'''
        if not 'TCI' in in_ras_file:
            in_ras = os.path.join(img_dir, in_ras_file)
            # in_band = in_ras.GetRasterBand(1)
            # out_rows = in_band.YSize * 2
            # out_columns = in_band.XSize * 2
            #
            # tmp_ras = gdal.GetDriverByName('MEM').Create('', out_columns, out_rows, 1, in_band.DataType)
            # tmp_ras.SetProjection(in_ras.GetProjection())
            # geotransform = list(in_ras.GetGeoTransform())
            # geotransform[1] /= 2
            # geotransform[5] /= 2
            # tmp_ras.SetGeoTransform(geotransform)
            # data = in_band.ReadAsArray(buf_xsize=out_columns, buf_ysize=out_rows)
            # out_band = tmp_ras.GetRasterBand(1)
            # out_band.WriteArray(data)
            #
            # drv = gdal.GetDriverByName('GTiff')
            # drv.CreateCopy(os.path.join(img_dir, out_ras_file), tmp_ras)
            #
            # del drv
            # del tmp_ras
            # del in_ras
            # try vrt
            # os.system('gdalwarp -t_srs EPSG:3857 -dstnodata -9999.0 -tr 10.0 10.0 -r near -of vrt %s %s' %(in_ras,os.path.join(img_dir, out_ras_file.split('.')[0]+'.vr')))
            # old, working
            gdal.Warp(os.path.join(img_dir, out_ras_file),in_ras, xRes=10,yRes=10)

        Sentinel2_pre_process.logger.info('Resampling of %s has been successfully finished.', out_ras_file)

    def _get_path_bands(self):
        '''Get the path for each band - returns a dictionary'''

        try:
            metadta_xml = ET.parse(os.path.join(self.ws, 'MTD_MSIL2A.xml'))
            metadta_xml_root = metadta_xml.getroot()
            bands = {}

            cnt = 0
            for child in metadta_xml_root.iter('IMAGE_FILE_2A'):
                cnt+=1
            if cnt == 0:
                treename = 'IMAGE_FILE'
            else:
                treename = 'IMAGE_FILE_2A'

            for child in metadta_xml_root.iter(treename):
                for prf, b in Sentinel2_pre_process.prefixes:
                    if prf in child.text:
                        logging.info(b)
                        image_dir = os.path.dirname(os.path.join(self.ws, child.text))
                        in_image_name = child.text.split('/')[-1] + '.jp2'
                        if '20m' in child.text.split('/')[-1]:
                            out_image_name = child.text.split('/')[-1].split('_20m')[0] + '_10m.tif'
                        else:
                            out_image_name = child.text.split('/')[-1] + '.tif'

                        # changed for rodopi 2017
                        #out_proj_image_name = out_image_name.split('.tif')[0] + '_prj.' + out_image_name.split('.')[1]

                        out_proj_image_name = out_image_name.split('.')[0] + '_prj.' + out_image_name.split('.')[1]

                        bands[b] = [image_dir, in_image_name, out_image_name, out_proj_image_name]
            self.logger.info('Directory pathname for the specified bands of %s have been successfully initialized', self.product)
            return bands

        except:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            self.logger.error('Line: %d,\t Type: %s,\t Message: %s', exc_tb.tb_lineno, exc_type, exc_obj)


class Sentinel2_process(MyLogger):
    ''' Processing of Sentinel 2 products '''

    prefixes = {'B02': 'blue', 'B03': 'green', 'B04': 'red', 'B05': 'vre1', 'B06': 'vre2', 'B07': 'vre3', 'B08': 'nir', 'B8A': 'nir_vre', 'B11': 'swir1', 'B12': 'swir2', 'SCL': 'cloud','RGB': 'tci'}

    def __init__(self, ws):
        self.ws = ws
        self.bands = self._get_bands()


    def maskcloud_lth(self):

        months = {'01': 'January_', '02': 'February_', '03': 'March_', '04': 'April_', '05': 'May_', '06': 'June_',
                  '07': 'July_', '08': 'August_', '09': 'September_', '10': 'October_',
                  '11': 'November_', '12': 'December_'}

        cloud_ds = gdal.Open(self.bands['cloud'])

        x = cloud_ds.GetRasterBand(1)

        cloud_band = x.ReadAsArray()

        cloud = np.array(cloud_band, dtype=int)
        shape = cloud_band.shape

        del cloud_band

        # xsize, ysize = red_band.XSize, red_band.YSize
        # block_xsize, block_ysize = red_band.GetBlockSize()
        # check = np.logical_or(red > 0, nir > 0)
        # cloud = np.where(cloud<=7 & cloud>3 & cloud<3,0,1)
        cloud = np.where(np.logical_or(cloud > 7, cloud == 3), 1, 0)

        geo = cloud_ds.GetGeoTransform()
        proj = cloud_ds.GetProjection()
        driver = gdal.GetDriverByName("GTiff")
        dst_ds = driver.Create(
            os.path.join(self.ws, os.path.split(self.bands['cloud'])[1].split('_')[0] + '_CloudMask.tif'), shape[1],
            shape[0], 1, gdal.GDT_Float32)
        dst_ds.SetGeoTransform(geo)
        dst_ds.SetProjection(proj)
        dst_ds.GetRasterBand(1).WriteArray(cloud)
        del cloud
        del dst_ds
        del cloud_ds

        logging.info("Cloud mask has been generated")

    def maskcloud(self):

        months = {'01': 'January_', '02': 'February_', '03': 'March_', '04': 'April_', '05': 'May_',       '06': 'June_', '07': 'July_', '08': 'August_', '09': 'September_', '10': 'October_',
                  '11': 'November_', '12': 'December_'}

        cloud_ds = gdal.Open(self.bands['cloud'])
        cloud_band = cloud_ds.GetRasterBand(1)

        band_ds = gdal.Open(self.bands['cloud'])
        band_band = band_ds.GetRasterBand(1)

        xsize, ysize = band_band.XSize, band_band.YSize
        block_xsize, block_ysize = band_band.GetBlockSize()

        out_data = np.full((ysize, xsize), 0)

        for x in range(0, xsize, block_xsize):
            if x + block_xsize < xsize:
                cols = block_xsize
            else:
                cols = xsize - x
            for y in range(0, ysize, block_ysize):
                if y + block_ysize < ysize:
                    rows = block_ysize
                else:
                    rows = ysize - y

                cloud = cloud_band.ReadAsArray(x, y, cols, rows).astype(int)
                #cloud = np.where(cloud<=7 & cloud>3 & cloud<3,0,1)
                cloud = np.where(np.logical_or(cloud>7,cloud==3),1,0)

                #cloud1[cloud1>=7] = 0
                #cloud1[cloud1<7] = 1
                #band1 = band_band.ReadAsArray(x, y, cols, rows).astype(int)
                #ndvi = cloud1 * band1
                #ndvi = ndvi.filled(-9999)

                out_data[y:y + rows, x:cols] = cloud

            #out_psri = self.create_raster(gdal.Open(self.bands[item]), os.path.join(self.ws, os.path.split(
            #self.bands[item])[1].split('_')[0]  + '_' + os.path.split(
            #self.bands[item])[1].split('_')[1]), out_data, gdal.GDT_Float32, -9999)


        out_ndvi = self.create_raster(gdal.Open(self.bands['cloud']), os.path.join(self.ws,os.path.split(self.bands['cloud'])[1].split('_')[0] + '_CloudMask.tif'), out_data,gdal.GDT_Float32, -9999)

        del cloud_ds
        del band_ds
        del out_ndvi

    # cloud mask for phenology
    def maskcloudphenology(self):

        months = {'01': 'January_', '02': 'February_', '03': 'March_', '04': 'April_', '05': 'May_', '06': 'June_',
                  '07': 'July_', '08': 'August_', '09': 'September_', '10': 'October_',
                  '11': 'November_', '12': 'December_'}

        cloud_ds = gdal.Open(self.bands['cloud'])
        cloud_band = cloud_ds.GetRasterBand(1)

        band_ds = gdal.Open(self.bands['cloud'])
        band_band = band_ds.GetRasterBand(1)

        xsize, ysize = band_band.XSize, band_band.YSize
        block_xsize, block_ysize = band_band.GetBlockSize()

        out_data = np.full((ysize, xsize), 0)

        for x in range(0, xsize, block_xsize):
            if x + block_xsize < xsize:
                cols = block_xsize
            else:
                cols = xsize - x
            for y in range(0, ysize, block_ysize):
                if y + block_ysize < ysize:
                    rows = block_ysize
                else:
                    rows = ysize - y

                cloud = cloud_band.ReadAsArray(x, y, cols, rows).astype(int)
                # cloud = np.where(cloud<=7 & cloud>3 & cloud<3,0,1)


                # this is the rule based on the SLC raster
                cloud = np.where(np.logical_or(cloud>=7,np.logical_or(cloud==2,cloud==3)),1,0)

                # cloud1[cloud1>=7] = 0
                # cloud1[cloud1<7] = 1
                # band1 = band_band.ReadAsArray(x, y, cols, rows).astype(int)
                # ndvi = cloud1 * band1
                # ndvi = ndvi.filled(-9999)

                out_data[y:y + rows, x:cols] = cloud

                # out_psri = self.create_raster(gdal.Open(self.bands[item]), os.path.join(self.ws, os.path.split(
                # self.bands[item])[1].split('_')[0]  + '_' + os.path.split(
                # self.bands[item])[1].split('_')[1]), out_data, gdal.GDT_Float32, -9999)

        out_ndvi = self.create_raster(gdal.Open(self.bands['cloud']), os.path.join(self.ws,
                                                                                   os.path.split(self.bands['cloud'])[
                                                                                       1].split('_')[
                                                                                       0] + '_CloudMaskPhenology.tif'),
                                      out_data, gdal.GDT_Float32, -9999)

        del cloud_ds
        del band_ds
        del out_ndvi

        # cloud mask for phenology

    def maskupdated(self):

        months = {'01': 'January_', '02': 'February_', '03': 'March_', '04': 'April_', '05': 'May_', '06': 'June_',
                  '07': 'July_', '08': 'August_', '09': 'September_', '10': 'October_',
                  '11': 'November_', '12': 'December_'}

        cloud_ds = gdal.Open(self.bands['cloud'])
        cloud_band = cloud_ds.GetRasterBand(1)

        band_ds = gdal.Open(self.bands['cloud'])
        band_band = band_ds.GetRasterBand(1)

        xsize, ysize = band_band.XSize, band_band.YSize
        block_xsize, block_ysize = band_band.GetBlockSize()

        out_data = np.full((ysize, xsize), 0)

        for x in range(0, xsize, block_xsize):
            if x + block_xsize < xsize:
                cols = block_xsize
            else:
                cols = xsize - x
            for y in range(0, ysize, block_ysize):
                if y + block_ysize < ysize:
                    rows = block_ysize
                else:
                    rows = ysize - y

                cloud = cloud_band.ReadAsArray(x, y, cols, rows).astype(int)
                # cloud = np.where(cloud<=7 & cloud>3 & cloud<3,0,1)


                # this is the rule based on the SLC raster
                cloud = np.where(np.logical_and(cloud>=4,cloud<=6), 0, 1)
                # cloud1[cloud1>=7] = 0
                # cloud1[cloud1<7] = 1
                # band1 = band_band.ReadAsArray(x, y, cols, rows).astype(int)
                # ndvi = cloud1 * band1
                # ndvi = ndvi.filled(-9999)

                out_data[y:y + rows, x:cols] = cloud

                # out_psri = self.create_raster(gdal.Open(self.bands[item]), os.path.join(self.ws, os.path.split(
                # self.bands[item])[1].split('_')[0]  + '_' + os.path.split(
                # self.bands[item])[1].split('_')[1]), out_data, gdal.GDT_Float32, -9999)

        out_ndvi = self.create_raster(gdal.Open(self.bands['cloud']), os.path.join(self.ws,
                                                                                   os.path.split(self.bands['cloud'])[
                                                                                       1].split('_')[
                                                                                       0] + '_CloudMaskUpdated.tif'),
                                      out_data, gdal.GDT_Int16, -9999)

        del cloud_ds
        del band_ds
        del out_ndvi

    def maskcloudvegatation(self):

        months = {'01': 'January_', '02': 'February_', '03': 'March_', '04': 'April_', '05': 'May_', '06': 'June_',
                  '07': 'July_', '08': 'August_', '09': 'September_', '10': 'October_',
                  '11': 'November_', '12': 'December_'}



        band_ds = gdal.Open(self.bands['cloud'])
        band_band = band_ds.GetRasterBand(1)

        xsize, ysize = band_band.XSize, band_band.YSize
        block_xsize, block_ysize = band_band.GetBlockSize()

        out_data = np.full((ysize, xsize), 0)

        for x in range(0, xsize, block_xsize):
            if x + block_xsize < xsize:
                cols = block_xsize
            else:
                cols = xsize - x
            for y in range(0, ysize, block_ysize):
                if y + block_ysize < ysize:
                    rows = block_ysize
                else:
                    rows = ysize - y

                cloud = band_band.ReadAsArray(x, y, cols, rows).astype(int)
                # cloud = np.where(cloud<=7 & cloud>3 & cloud<3,0,1)
                cloud = np.where(cloud==4,0,1)

                # cloud1[cloud1>=7] = 0
                # cloud1[cloud1<7] = 1
                # band1 = band_band.ReadAsArray(x, y, cols, rows).astype(int)
                # ndvi = cloud1 * band1
                # ndvi = ndvi.filled(-9999)

                out_data[y:y + rows, x:cols] = cloud

                # out_psri = self.create_raster(gdal.Open(self.bands[item]), os.path.join(self.ws, os.path.split(
                # self.bands[item])[1].split('_')[0]  + '_' + os.path.split(
                # self.bands[item])[1].split('_')[1]), out_data, gdal.GDT_Float32, -9999)

        out_ndvi = self.create_raster(gdal.Open(self.bands['cloud']), os.path.join(self.ws,
                                                                                   os.path.split(self.bands['cloud'])[
                                                                                       1].split('_')[
                                                                                       0] + '_CloudMaskVegetation.tif'),
                                      out_data, gdal.GDT_Float32, -9999)

        del band_ds
        del out_ndvi


    def maskcloud_hard(self):

        months = {'01': 'January_', '02': 'February_', '03': 'March_', '04': 'April_', '05': 'May_',       '06': 'June_', '07': 'July_', '08': 'August_', '09': 'September_', '10': 'October_',
                  '11': 'November_', '12': 'December_'}

        cloud_ds = gdal.Open(self.bands['cloud'])
        cloud_band = cloud_ds.GetRasterBand(1)

        band_ds = gdal.Open(self.bands['cloud'])
        band_band = band_ds.GetRasterBand(1)

        xsize, ysize = band_band.XSize, band_band.YSize
        block_xsize, block_ysize = band_band.GetBlockSize()

        out_data = np.full((ysize, xsize), 0)

        for x in range(0, xsize, block_xsize):
            if x + block_xsize < xsize:
                cols = block_xsize
            else:
                cols = xsize - x
            for y in range(0, ysize, block_ysize):
                if y + block_ysize < ysize:
                    rows = block_ysize
                else:
                    rows = ysize - y

                cloud = cloud_band.ReadAsArray(x, y, cols, rows).astype(int)
                #cloud = np.where(cloud<=7 & cloud>3 & cloud<3,0,1)
                #cloud = np.where(np.logical_or(cloud>8,cloud==3),1,0)

                #
                cloud = np.where(cloud>8,1,0)

                #cloud1[cloud1>=7] = 0
                #cloud1[cloud1<7] = 1
                #band1 = band_band.ReadAsArray(x, y, cols, rows).astype(int)
                #ndvi = cloud1 * band1
                #ndvi = ndvi.filled(-9999)

                out_data[y:y + rows, x:cols] = cloud

            #out_psri = self.create_raster(gdal.Open(self.bands[item]), os.path.join(self.ws, os.path.split(
            #self.bands[item])[1].split('_')[0]  + '_' + os.path.split(
            #self.bands[item])[1].split('_')[1]), out_data, gdal.GDT_Float32, -9999)


        out_ndvi = self.create_raster(gdal.Open(self.bands['cloud']), os.path.join(self.ws,os.path.split(self.bands['cloud'])[1].split('_')[0] + '_CloudMaskHard.tif'), out_data,gdal.GDT_Float32, -9999)

        del cloud_ds
        del band_ds
        del out_ndvi

    def maskcloud_light(self):

        months = {'01': 'January_', '02': 'February_', '03': 'March_', '04': 'April_', '05': 'May_',       '06': 'June_', '07': 'July_', '08': 'August_', '09': 'September_', '10': 'October_',
                  '11': 'November_', '12': 'December_'}

        cloud_ds = gdal.Open(self.bands['cloud'])
        cloud_band = cloud_ds.GetRasterBand(1)

        band_ds = gdal.Open(self.bands['cloud'])
        band_band = band_ds.GetRasterBand(1)

        xsize, ysize = band_band.XSize, band_band.YSize
        block_xsize, block_ysize = band_band.GetBlockSize()

        out_data = np.full((ysize, xsize), 0)

        for x in range(0, xsize, block_xsize):
            if x + block_xsize < xsize:
                cols = block_xsize
            else:
                cols = xsize - x
            for y in range(0, ysize, block_ysize):
                if y + block_ysize < ysize:
                    rows = block_ysize
                else:
                    rows = ysize - y

                cloud = cloud_band.ReadAsArray(x, y, cols, rows).astype(int)
                #cloud = np.where(cloud<=7 & cloud>3 & cloud<3,0,1)
                #cloud = np.where(np.logical_or(cloud>8,cloud==3),1,0)
                cloud = np.where(cloud>7,1,0)

                #cloud1[cloud1>=7] = 0
                #cloud1[cloud1<7] = 1
                #band1 = band_band.ReadAsArray(x, y, cols, rows).astype(int)
                #ndvi = cloud1 * band1
                #ndvi = ndvi.filled(-9999)

                out_data[y:y + rows, x:cols] = cloud

            #out_psri = self.create_raster(gdal.Open(self.bands[item]), os.path.join(self.ws, os.path.split(
            #self.bands[item])[1].split('_')[0]  + '_' + os.path.split(
            #self.bands[item])[1].split('_')[1]), out_data, gdal.GDT_Float32, -9999)


        out_ndvi = self.create_raster(gdal.Open(self.bands['cloud']), os.path.join(self.ws,os.path.split(self.bands['cloud'])[1].split('_')[0] + '_CloudMaskLight.tif'), out_data,gdal.GDT_Float32, -9999)

        del cloud_ds
        del band_ds
        del out_ndvi

    def maskcloud_update(self):

        months = {'01': 'January_', '02': 'February_', '03': 'March_', '04': 'April_', '05': 'May_', '06': 'June_',
                  '07': 'July_', '08': 'August_', '09': 'September_', '10': 'October_',
                  '11': 'November_', '12': 'December_'}

        cloud_ds = gdal.Open(self.bands['cloud'])
        cloud_band = cloud_ds.GetRasterBand(1)

        band_ds = gdal.Open(self.bands['cloud'])
        band_band = band_ds.GetRasterBand(1)

        xsize, ysize = band_band.XSize, band_band.YSize
        block_xsize, block_ysize = band_band.GetBlockSize()

        out_data = np.full((ysize, xsize), 0)

        for x in range(0, xsize, block_xsize):
            if x + block_xsize < xsize:
                cols = block_xsize
            else:
                cols = xsize - x
            for y in range(0, ysize, block_ysize):
                if y + block_ysize < ysize:
                    rows = block_ysize
                else:
                    rows = ysize - y

                cloud = cloud_band.ReadAsArray(x, y, cols, rows).astype(int)
                # cloud = np.where(cloud<=7 & cloud>3 & cloud<3,0,1)
                # cloud = np.where(np.logical_or(cloud>8,cloud==3),1,0)
                cloud = np.where(np.logical_or(cloud < 3, cloud > 8), 1, 0)

                # cloud1[cloud1>=7] = 0
                # cloud1[cloud1<7] = 1
                # band1 = band_band.ReadAsArray(x, y, cols, rows).astype(int)
                # ndvi = cloud1 * band1
                # ndvi = ndvi.filled(-9999)

                out_data[y:y + rows, x:cols] = cloud

                # out_psri = self.create_raster(gdal.Open(self.bands[item]), os.path.join(self.ws, os.path.split(
                # self.bands[item])[1].split('_')[0]  + '_' + os.path.split(
                # self.bands[item])[1].split('_')[1]), out_data, gdal.GDT_Float32, -9999)

        out_ndvi = self.create_raster(gdal.Open(self.bands['cloud']), os.path.join(self.ws,
                                                                                   os.path.split(self.bands['cloud'])[
                                                                                       1].split('_')[
                                                                                       0] + '_CloudMaskUpdate.tif'),
                                      out_data, gdal.GDT_Float32, -9999)

        del cloud_ds
        del band_ds
        del out_ndvi

    def shadow(self):

        months = {'01': 'January_', '02': 'February_', '03': 'March_', '04': 'April_', '05': 'May_', '06': 'June_',
                  '07': 'July_', '08': 'August_', '09': 'September_', '10': 'October_',
                  '11': 'November_', '12': 'December_'}

        cloud_ds = gdal.Open(self.bands['cloud'])
        cloud_band = cloud_ds.GetRasterBand(1)

        band_ds = gdal.Open(self.bands['cloud'])
        band_band = band_ds.GetRasterBand(1)

        xsize, ysize = band_band.XSize, band_band.YSize
        block_xsize, block_ysize = band_band.GetBlockSize()

        out_data = np.full((ysize, xsize), 0)

        for x in range(0, xsize, block_xsize):
            if x + block_xsize < xsize:
                cols = block_xsize
            else:
                cols = xsize - x
            for y in range(0, ysize, block_ysize):
                if y + block_ysize < ysize:
                    rows = block_ysize
                else:
                    rows = ysize - y

                cloud = cloud_band.ReadAsArray(x, y, cols, rows).astype(int)
                # cloud = np.where(cloud<=7 & cloud>3 & cloud<3,0,1)
                cloud = np.where(np.logical_or(cloud<=3, np.logical_or(cloud==6,cloud==7)), 1, 0)

                # cloud1[cloud1>=7] = 0
                # cloud1[cloud1<7] = 1
                # band1 = band_band.ReadAsArray(x, y, cols, rows).astype(int)
                # ndvi = cloud1 * band1
                # ndvi = ndvi.filled(-9999)

                out_data[y:y + rows, x:cols] = cloud

                # out_psri = self.create_raster(gdal.Open(self.bands[item]), os.path.join(self.ws, os.path.split(
                # self.bands[item])[1].split('_')[0]  + '_' + os.path.split(
                # self.bands[item])[1].split('_')[1]), out_data, gdal.GDT_Float32, -9999)

        out_ndvi = self.create_raster(gdal.Open(self.bands['cloud']), os.path.join(self.ws, os.path.split(
            self.bands['cloud'])[1].split('_')[0] + '_Shadows.tif'), out_data, gdal.GDT_Float32, -9999)

        del cloud_ds
        del band_ds
        del out_ndvi

    def shadow_lth(self):
        months = {'01': 'January_', '02': 'February_', '03': 'March_', '04': 'April_', '05': 'May_', '06': 'June_',
                  '07': 'July_', '08': 'August_', '09': 'September_', '10': 'October_',
                  '11': 'November_', '12': 'December_'}

        red_ds = gdal.Open(self.bands['cloud'])

        x = red_ds.GetRasterBand(1)

        red_band = x.ReadAsArray()

        red = np.array(red_band, dtype=int)
        shape = red_band.shape

        del red_band

        # xsize, ysize = red_band.XSize, red_band.YSize
        # cloud = np.where(cloud<=7 & cloud>3 & cloud<3,0,1)
        cloud = np.where(np.logical_or(red <= 3, np.logical_or(red == 6, red == 7)), 1, 0)

        # ndvi = (nir - red) / (nir + red) * 1000
        # ndvi = ndvi.filled(-9999)
        del red
        geo = red_ds.GetGeoTransform()
        proj = red_ds.GetProjection()
        driver = gdal.GetDriverByName("GTiff")
        dst_ds = driver.Create(
            os.path.join(self.ws, os.path.split(self.bands['red'])[1].split('_')[0] + '_Shadows.tif'), shape[1],
            shape[0], 1, gdal.GDT_Int16)
        dst_ds.SetGeoTransform(geo)
        dst_ds.SetProjection(proj)
        dst_ds.GetRasterBand(1).WriteArray(cloud)

        del red_ds
        del dst_ds
        del cloud

    def create_ndvi(self):
        ''' Returns Normalized Difference Vegetation Index '''
        red_ds = gdal.Open(self.bands['red'])
        nir_ds = gdal.Open(self.bands['nir'])

        red_band = red_ds.GetRasterBand(1)
        nir_band = nir_ds.GetRasterBand(1)

        xsize, ysize = red_band.XSize, red_band.YSize
        block_xsize, block_ysize = red_band.GetBlockSize()

        out_data = np.full((ysize, xsize), 0)
        for x in range(0, xsize, block_xsize):
            if x + block_xsize < xsize:
                cols = block_xsize
            else:
                cols = xsize - x
            for y in range(0, ysize, block_ysize):
                if y + block_ysize < ysize:
                    rows = block_ysize
                else:
                    rows = ysize - y

                red = red_band.ReadAsArray(x, y, cols, rows).astype(np.float32)
                nir = nir_band.ReadAsArray(x, y, cols, rows).astype(np.float32)

                red = np.ma.masked_where(nir + red == 0, red)
                ndvi = ((nir - red) / (nir + red)) * 1000
                ndvi = ndvi.filled(-9999)

                out_data[y:y+rows, x:cols] = ndvi

        out_ndvi = self.create_raster(gdal.Open(self.bands['red']), os.path.join(self.ws, os.path.split(self.bands['red'])[1].split('_')[0] + '_NDVI_INDEX.tif'), out_data, gdal.GDT_Int16, -9999)

        self.logger.info('Normalized Difference Vegetation Index (NDVI) calculation has been successfully finished')

        del red_ds
        del nir_ds
        del out_ndvi

    def create_ndvi_lth(self):
        ''' Returns Normalized Difference Vegetation Index '''
        red_ds = gdal.Open(self.bands['red'])
        nir_ds = gdal.Open(self.bands['nir'])

        x = red_ds.GetRasterBand(1)
        y = nir_ds.GetRasterBand(1)

        red_band = x.ReadAsArray()
        nir_band = y.ReadAsArray()

        red = np.array(red_band, dtype="float32")
        nir = np.array(nir_band, dtype="float32")
        shape = red_band.shape

        del red_band
        del nir_band

        # xsize, ysize = red_band.XSize, red_band.YSize
        # block_xsize, block_ysize = red_band.GetBlockSize()
        check = np.logical_or(red > 0, nir > 0)
        #red = np.ma.masked_where(nir + red == 0, red)
        ndvi = np.where(check, (nir - red) / (nir + red)*1000, -9999)
        del check
        #ndvi = (nir - red) / (nir + red) * 1000
        #ndvi = ndvi.filled(-9999)
        del red
        del nir
        geo = red_ds.GetGeoTransform()
        proj = red_ds.GetProjection()
        driver = gdal.GetDriverByName("GTiff")
        dst_ds = driver.Create(
            os.path.join(self.ws, os.path.split(self.bands['red'])[1].split('_')[0] + '_NDVI_INDEX.tif'), shape[1],
            shape[0], 1, gdal.GDT_Float32)
        dst_ds.SetGeoTransform(geo)
        dst_ds.SetProjection(proj)
        dst_ds.GetRasterBand(1).WriteArray(ndvi)

        del red_ds
        del nir_ds
        del dst_ds
        del ndvi

        '''
        out_data = np.full((ysize, xsize), 0)
        for x in range(0, xsize, block_xsize):
            if x + block_xsize < xsize:
                cols = block_xsize
            else:
                cols = xsize - x
            for y in range(0, ysize, block_ysize):
                if y + block_ysize < ysize:
                    rows = block_ysize
                else:
                    rows = ysize - y

                red = red_band.ReadAsArray(x, y, cols, rows).astype(np.float)
                nir = nir_band.ReadAsArray(x, y, cols, rows).astype(np.float)

                red = np.ma.masked_where(nir + red == 0, red)
                ndvi = ((nir - red) / (nir + red)) * 1000
                ndvi = ndvi.filled(-9999)

                out_data[y:y+rows, x:cols] = ndvi
        '''
        # out_ndvi = self.create_raster(gdal.Open(self.bands['red']), os.path.join(self.ws, os.path.split(self.bands['red'])[1].split('_')[0] + '_NDVI_INDEX.tif'), out_data, gdal.GDT_Float32, -9999)

        self.logger.info('Normalized Difference Vegetation Index (NDVI) calculation has been successfully finished')

    def create_bsm(self):

        #ndvi calculation
        red_ds = gdal.Open(self.bands['red'])
        nir_ds = gdal.Open(self.bands['nir'])
        x = red_ds.GetRasterBand(1)
        y = nir_ds.GetRasterBand(1)
        red_band = x.ReadAsArray()
        del red_ds
        del nir_ds
        nir_band = y.ReadAsArray()
        del x
        del y
        red = np.array(red_band, dtype="float16")
        shape = red_band.shape
        del red_band
        nir = np.array(nir_band, dtype="float16")
        del nir_band
        check = np.logical_or(red > 0, nir > 0)
        ndvi = np.where(check, (nir - red) / (nir + red), -9999)
        del check
        del red
        del nir


        #ndwi calculation
        green_ds = gdal.Open(self.bands['green'])
        nir_ds = gdal.Open(self.bands['swir2'])
        x = green_ds.GetRasterBand(1)
        y = nir_ds.GetRasterBand(1)
        del green_ds
        del nir_ds
        green_band = x.ReadAsArray()
        nir_band = y.ReadAsArray()
        del x
        del y
        green = np.array(green_band, dtype="float16")
        del green_band
        nir = np.array(nir_band, dtype="float16")
        del nir_band
        check = np.logical_or(green > 0, nir > 0)
        check = np.logical_or(green > 0, nir > 0)
        ndwi = np.where(check, (green - nir) / (green + nir), -9999)
        del check
        del nir
        del green

        #lidx
        swir_ds = gdal.Open(self.bands['swir1'])
        swir_ds2 = gdal.Open(self.bands['swir2'])
        x = swir_ds.GetRasterBand(1)
        y = swir_ds2.GetRasterBand(1)
        del swir_ds2
        del swir_ds
        swir_band = x.ReadAsArray()
        swir2_band = y.ReadAsArray()
        del x
        del y
        swir = np.array(swir_band, dtype="float16")
        del swir_band
        swir2 = np.array(swir2_band, dtype="float16")
        del swir2_band
        check = swir2 > 0
        lidx = np.where(check, swir/swir2, -9999)
        del check
        del swir2
        del swir


        #nbr
        nir_ds = gdal.Open(self.bands['nir'])
        swir2_ds = gdal.Open(self.bands['swir2'])
        x = swir2_ds.GetRasterBand(1)
        y = nir_ds.GetRasterBand(1)
        del nir_ds
        geo = swir2_ds.GetGeoTransform()
        proj = swir2_ds.GetProjection()
        del swir2_ds
        swir2_band = x.ReadAsArray()
        nir_band = y.ReadAsArray()
        del x
        del y
        swir2 = np.array(swir2_band, dtype="float16")
        del swir2_band
        nir = np.array(nir_band, dtype="float16")
        del nir_band
        check = np.logical_or(swir2 > 0, nir > 0)
        nbr = np.where(check, (nir - swir2) / (nir + swir2), -9999)
        del check
        del swir2
        del nir


        #numexpr = '(1.3*nbr + 0.6*ndvi + 0.60*ndwi+ 0.50*lidx) < -0.32'
        #SDIThreshold = np.evaluate(numexpr).astype(np.float32)
        std1= np.std(nbr[nbr > -9999])
        nbr = nbr[nbr>-9999]/std1
        std2 = np.std(ndvi[ndvi > -9999])
        ndvi = ndvi[ndvi > -9999] / std2
        std3 = np.std(ndwi[ndwi > -9999])
        ndwi = ndwi[ndwi > -9999] / std3
        std4 = np.std(lidx[lidx > -9999])
        lidx = lidx[lidx > -9999] / std4

        check = (1.3*nbr + 0.6*ndvi + 0.60*ndwi+ 0.50*lidx)
        del nbr,ndvi,ndwi,lidx
        #SDIThreshold = np.where(check,1,0)


        driver = gdal.GetDriverByName("GTiff")
        dst_ds = driver.Create(os.path.join(self.ws, os.path.split(self.bands['red'])[1].split('_')[0] + '_bsm.tif'),shape[1],shape[0], 1, gdal.GDT_Float32)
        dst_ds.SetGeoTransform(geo)
        dst_ds.SetProjection(proj)
        dst_ds.GetRasterBand(1).WriteArray(check)


    def create_psri_lth(self):
        red_ds = gdal.Open(self.bands['red'])
        blue_ds = gdal.Open(self.bands['blue'])
        vre1_ds = gdal.Open(self.bands['vre2'])

        x = red_ds.GetRasterBand(1)
        y = blue_ds.GetRasterBand(1)
        z = vre1_ds.GetRasterBand(1)

        red_band = x.ReadAsArray()
        blue_band = y.ReadAsArray()
        vre1_band = z.ReadAsArray()

        red = np.array(red_band, dtype=float)
        blue = np.array(blue_band, dtype=float)
        vre1 = np.array(vre1_band, dtype=float)
        shape = red_band.shape

        del red_band
        del blue_band
        del vre1_band


        #red = np.ma.masked_where(vre1== 0, vre1)
        # ndvi = np.where(check, (nir - red) / (nir + red)*1000, -9999)
        #psri = ((red - blue) / vre1) * 1000
        check = vre1 > 0
        # red = np.ma.masked_where(nir + red == 0, red)
        psri = np.where(check, (red - blue) / vre1 * 1000, -9999)
        del check
        del red
        del blue
        del vre1

        geo = red_ds.GetGeoTransform()
        proj = red_ds.GetProjection()
        driver = gdal.GetDriverByName("GTiff")
        dst_ds = driver.Create(os.path.join(self.ws, os.path.split(self.bands['red'])[1].split('_')[0] + '_PSRI_INDEX.tif'), shape[1],
            shape[0], -9999, gdal.GDT_Float32)
        dst_ds.SetGeoTransform(geo)
        dst_ds.SetProjection(proj)
        dst_ds.GetRasterBand(1).WriteArray(psri)

        del red_ds
        del blue_ds
        del vre1_ds
        # del out_ndvi
        del dst_ds
        del psri

    def create_psri(self):
        ''' Returns Plant Senescence Reflectance Index '''

        red_ds = gdal.Open(self.bands['red'])
        blue_ds = gdal.Open(self.bands['blue'])
        vre1_ds = gdal.Open(self.bands['vre2'])

        red_band = red_ds.GetRasterBand(1)
        blue_band = blue_ds.GetRasterBand(1)
        vre1_band = vre1_ds.GetRasterBand(1)

        xsize, ysize = red_band.XSize, red_band.YSize
        block_xsize, block_ysize = red_band.GetBlockSize()

        out_data = np.full((ysize, xsize), 0)
        for x in range(0, xsize, block_xsize):
            if x + block_xsize < xsize:
                cols = block_xsize
            else:
                cols = xsize - x
            for y in range(0, ysize, block_ysize):
                if y + block_ysize < ysize:
                    rows = block_ysize
                else:
                    rows = ysize - y

                red = red_band.ReadAsArray(x, y, cols, rows).astype(np.float)
                blue = blue_band.ReadAsArray(x, y, cols, rows).astype(np.float)
                vre1 = vre1_band.ReadAsArray(x, y, cols, rows).astype(np.float)

                vre1 = np.ma.masked_where(vre1 == 0, vre1)
                psri = ((red - blue) / vre1) * 1000
                psri.filled(-9999)

                out_data[y:y+rows, x:cols] = psri

        out_psri = self.create_raster(gdal.Open(self.bands['red']), os.path.join(self.ws,  os.path.split(self.bands['red'])[1].split('_')[0] + '_PSRI_INDEX.tif'), out_data, gdal.GDT_Int32, -9999)

        self.logger.info('Plant Senescence Reflectance Index (PSRI) calculation has been successfully finished')

        del red_ds
        del blue_ds
        del vre1_ds
        del out_psri

    def create_bori(self):
        ''' Returns Plant Senescence Reflectance Index '''

        red_ds = gdal.Open(self.bands['red'])
        blue_ds = gdal.Open(self.bands['blue'])
        vre3_ds = gdal.Open(self.bands['vre3'])

        red_band = red_ds.GetRasterBand(1)
        blue_band = blue_ds.GetRasterBand(1)
        vre3_band = vre3_ds.GetRasterBand(1)

        xsize, ysize = red_band.XSize, red_band.YSize
        block_xsize, block_ysize = red_band.GetBlockSize()

        out_data = np.full((ysize, xsize), 0)
        for x in range(0, xsize, block_xsize):
            if x + block_xsize < xsize:
                cols = block_xsize
            else:
                cols = xsize - x
            for y in range(0, ysize, block_ysize):
                if y + block_ysize < ysize:
                    rows = block_ysize
                else:
                    rows = ysize - y

                red = red_band.ReadAsArray(x, y, cols, rows).astype(np.float)
                blue = blue_band.ReadAsArray(x, y, cols, rows).astype(np.float)
                vre3 = vre3_band.ReadAsArray(x, y, cols, rows).astype(np.float)


                bori = ((665-493)*(vre3-blue)-(783-493)*(red-blue))/float(2.0)

                out_data[y:y + rows, x:cols] = bori

        out_bori = self.create_raster(gdal.Open(self.bands['red']), os.path.join(self.ws,
                                                                                 os.path.split(self.bands['red'])[
                                                                                     1].split('_')[
                                                                                     0] + '_BORI_INDEX.tif'), out_data,
                                      gdal.GDT_Float32, -9999)

        self.logger.info('BORI Index calculation has been successfully finished')

        del red_ds
        del blue_ds
        del vre3_ds
        del out_bori

    def create_bari(self):
        ''' Returns BARI Index '''

        red_ds = gdal.Open(self.bands['red'])
        blue_ds = gdal.Open(self.bands['blue'])
        vre3_ds = gdal.Open(self.bands['vre3'])

        red_band = red_ds.GetRasterBand(1)
        blue_band = blue_ds.GetRasterBand(1)
        vre3_band = vre3_ds.GetRasterBand(1)

        xsize, ysize = red_band.XSize, red_band.YSize
        block_xsize, block_ysize = red_band.GetBlockSize()

        out_data = np.full((ysize, xsize), 0)
        for x in range(0, xsize, block_xsize):
            if x + block_xsize < xsize:
                cols = block_xsize
            else:
                cols = xsize - x
            for y in range(0, ysize, block_ysize):
                if y + block_ysize < ysize:
                    rows = block_ysize
                else:
                    rows = ysize - y

                red = red_band.ReadAsArray(x, y, cols, rows).astype(np.float)
                blue = blue_band.ReadAsArray(x, y, cols, rows).astype(np.float)
                vre3 = vre3_band.ReadAsArray(x, y, cols, rows).astype(np.float)


                bari = ((665-493)*(vre3-blue)-(783-493)*(red-blue))/(2.0*red)

                out_data[y:y + rows, x:cols] = bari

        out_bari = self.create_raster(gdal.Open(self.bands['red']), os.path.join(self.ws,
                                                                                 os.path.split(self.bands['red'])[
                                                                                     1].split('_')[
                                                                                     0] + '_BARI_INDEX.tif'), out_data,
                                      gdal.GDT_Float32, -9999)

        self.logger.info('BARI Index calculation has been successfully finished')

        del red_ds
        del blue_ds
        del vre3_ds
        del out_bari

    def create_ndwi_lth(self):
        green_ds = gdal.Open(self.bands['swir1'])
        nir_ds = gdal.Open(self.bands['nir'])

        x = green_ds.GetRasterBand(1)
        y = nir_ds.GetRasterBand(1)


        green_band = x.ReadAsArray()
        nir_band = y.ReadAsArray()
        shape = green_band.shape

        if  shape[1] < nir_band.shape[1]:
            nir_band = nir_band[:,:shape[1]]
        else:
            green_band = green_band[:,:nir_band.shape[1]]

        nir = np.array(nir_band, dtype=float)
        green = np.array(green_band, dtype=float)

        del green_band
        del nir_band

        check = np.logical_or(green > 0, nir > 0)

        #nir = np.ma.masked_where(green + nir == 0, nir)
        ndwi = np.where(check, (nir - green) / (green + nir) * 1000, -9999)

        del check

        #ndwi = (nir - green) / (green + nir) * 1000
        #ndwi.filled(-9999)
        del nir
        del green

        out_ndwi = self.create_raster(gdal.Open(self.bands['green']), os.path.join(self.ws, os.path.split(self.bands['green'])[1].split('_')[0] + '_NDWI_INDEX.tif'), ndwi, gdal.GDT_Int16, -9999)

        # geo = green_ds.GetGeoTransform()
        # proj = green_ds.GetProjection()
        # driver = gdal.GetDriverByName("GTiff")
        # dst_ds = driver.Create(
        #     os.path.join(self.ws, os.path.split(self.bands['green'])[1].split('_')[0] + '_NDWI_INDEX.tif'), shape[1],
        #     shape[0], 1, gdal.GDT_Float32)
        # dst_ds.SetGeoTransform(geo)
        # dst_ds.SetProjection(proj)
        # dst_ds.GetRasterBand(1).WriteArray(ndwi)
        del out_ndwi
        del green_ds
        del nir_ds
        # del dst_ds
        del ndwi

    def create_ndwi(self):
        ''' Returns Normalized Difference Water Index '''

        green_ds = gdal.Open(self.bands['swir1'])
        nir_ds = gdal.Open(self.bands['nir'])

        green_band = green_ds.GetRasterBand(1)
        nir_band = nir_ds.GetRasterBand(1)

        xsize, ysize = green_band.XSize, green_band.YSize
        block_xsize, block_ysize = green_band.GetBlockSize()

        out_data = np.full((ysize, xsize), 0)
        for x in range(0, xsize, block_xsize):
            if x + block_xsize < xsize:
                cols = block_xsize
            else:
                cols = xsize - x
            for y in range(0, ysize, block_ysize):
                if y + block_ysize < ysize:
                    rows = block_ysize
                else:
                    rows = ysize - y

                green = green_band.ReadAsArray(x, y, cols, rows).astype(np.float)
                nir = nir_band.ReadAsArray(x, y, cols, rows).astype(np.float)
                nir = np.ma.masked_where(green + nir == 0, nir)
                ndwi = (nir - green) / (green + nir) * 1000
                ndwi.filled(-9999)

                out_data[y:y+rows, x:cols] = ndwi

        out_ndwi = self.create_raster(gdal.Open(self.bands['green']), os.path.join(self.ws, os.path.split(self.bands['green'])[1].split('_')[0] + '_NDWI_INDEX.tif'), out_data, gdal.GDT_Int16, -9999)

        self.logger.info('Normalized Difference Water Index (NDWI) calculation has been successfully finished')

        del green_ds
        del nir_ds
        del out_ndwi


    def create_savi_lth(self):
        red_ds = gdal.Open(self.bands['red'])
        nir_ds = gdal.Open(self.bands['nir'])

        x = red_ds.GetRasterBand(1)
        y = nir_ds.GetRasterBand(1)


        red_band = x.ReadAsArray()
        nir_band = y.ReadAsArray()

        red = np.array(red_band, dtype=float)
        nir = np.array(nir_band, dtype=float)
        shape = red_band.shape

        del red_band
        del nir_band

        check = (nir + red + 0.5) != 0

        # nir = np.ma.masked_where(green + nir == 0, nir)
        savi = np.where(check, (1.5 * (nir - red)) / (nir + red + 0.5), -9999)

        #savi = (1.5 * (nir - red)) / (nir + red + 0.5)
        #savi = savi.filled(-9999)
        del check
        del red
        del nir

        geo = red_ds.GetGeoTransform()
        proj = red_ds.GetProjection()
        driver = gdal.GetDriverByName("GTiff")
        dst_ds = driver.Create(
            os.path.join(self.ws, os.path.split(self.bands['red'])[1].split('_')[0] + '_SAVI_INDEX.tif'), shape[1],
            shape[0], 1, gdal.GDT_Float32)
        dst_ds.SetGeoTransform(geo)
        dst_ds.SetProjection(proj)
        dst_ds.GetRasterBand(1).WriteArray(savi)

        del red_ds
        del nir_ds
        del dst_ds
        del savi


    def create_savi(self):
        ''' Returns Soil Adjusted Vegetation Index '''

        red_ds = gdal.Open(self.bands['red'])
        nir_ds = gdal.Open(self.bands['nir'])

        red_band = red_ds.GetRasterBand(1)
        nir_band = nir_ds.GetRasterBand(1)

        xsize, ysize = red_band.XSize, red_band.YSize
        block_xsize, block_ysize = red_band.GetBlockSize()

        out_data = np.full((ysize, xsize), 0)
        for x in range(0, xsize, block_xsize):
            if x + block_xsize < xsize:
                cols = block_xsize
            else:
                cols = xsize - x
            for y in range(0, ysize, block_ysize):
                if y + block_ysize < ysize:
                    rows = block_ysize
                else:
                    rows = ysize - y

                red = red_band.ReadAsArray(x, y, cols, rows).astype(np.float)
                nir = nir_band.ReadAsArray(x, y, cols, rows).astype(np.float)

                red = np.ma.masked_where(nir + red + 0.5 == 0, red)
                savi = ((1.5 * (nir - red)) / (nir + red + 0.5))*1000
                savi = savi.filled(-9999)

                out_data[y:y+rows, x:cols] = savi

        out_savi = self.create_raster(gdal.Open(self.bands['red']), os.path.join(self.ws, os.path.split(self.bands['red'])[1].split('_')[0] + '_SAVI_INDEX.tif'), out_data, gdal.GDT_Int16, -9999)

        self.logger.info('Soil Adjusted Vegetation Index (SAVI) calculation has been successfully finished')

        del red_ds
        del nir_ds
        del out_savi

    def create_nbr_lth(self):
        nir_ds = gdal.Open(self.bands['nir'])
        swir2_ds = gdal.Open(self.bands['swir2'])

        x = swir2_ds.GetRasterBand(1)
        y = nir_ds.GetRasterBand(1)


        swir2_band = x.ReadAsArray()
        nir_band = y.ReadAsArray()

        swir2 = np.array(swir2_band, dtype=float)
        nir = np.array(nir_band, dtype=float)
        shape = swir2_band.shape

        del swir2_band
        del nir_band

        check = np.logical_or(swir2 > 0, nir > 0)

        # nir = np.ma.masked_where(green + nir == 0, nir)
        nbr = np.where(check, (nir - swir2) / (nir + swir2) * 1000, -9999)
        del check
        del swir2
        del nir
        geo = nir_ds.GetGeoTransform()
        proj = nir_ds.GetProjection()
        driver = gdal.GetDriverByName("GTiff")
        dst_ds = driver.Create(
            os.path.join(self.ws, os.path.split(self.bands['nir'])[1].split('_')[0] + '_NBR_INDEX.tif'), shape[1],
            shape[0], 1, gdal.GDT_Float32)
        dst_ds.SetGeoTransform(geo)
        dst_ds.SetProjection(proj)
        dst_ds.GetRasterBand(1).WriteArray(nbr)

        del swir2_ds
        del nir_ds
        del dst_ds
        del nbr

    def create_nbr(self):
        ''' Returns Normalized Burn Ratio Index '''

        nir_ds = gdal.Open(self.bands['nir'])
        swir2_ds = gdal.Open(self.bands['swir2'])

        nir_band = nir_ds.GetRasterBand(1)
        swir2_band = swir2_ds.GetRasterBand(1)

        xsize, ysize = nir_band.XSize, nir_band.YSize
        block_xsize, block_ysize = nir_band.GetBlockSize()

        out_data = np.full((ysize, xsize), 0)
        for x in range(0, xsize, block_xsize):
            if x + block_xsize < xsize:
                cols = block_xsize
            else:
                cols = xsize - x
            for y in range(0, ysize, block_ysize):
                if y + block_ysize < ysize:
                    rows = block_ysize
                else:
                    rows = ysize - y

                nir = nir_band.ReadAsArray(x, y, cols, rows).astype(np.float)
                swir2 = swir2_band.ReadAsArray(x, y, cols, rows).astype(np.float)

                nir = np.ma.masked_where(nir + swir2 ==0, nir)
                nbr = (nir - swir2)/(nir + swir2)
                nbr = nbr.filled(-9999)

                out_data[y:y+rows, x:cols] = nbr

        out_nbr = self.create_raster(gdal.Open(self.bands['nir']), os.path.join(self.ws,  os.path.split(self.bands['nir'])[1].split('_')[0] + '_NBR_INDEX.tif'), out_data, gdal.GDT_Float32, -9999)

        self.logger.info('Normalized Burn Ratio (NBR) calculation has been successfully finished')

        del nir_ds
        del swir2_ds
        del out_nbr

    def create_rgb_8bit(self):
        ''' Returns a 8bit True Color (RGB) raster '''
        red_ds = self.convert_to_8bit(self.bands['red'])
        red_band = red_ds.GetRasterBand(1)
        green_ds = self.convert_to_8bit(self.bands['green'])
        blue_ds = self.convert_to_8bit(self.bands['blue'])

        out_ds = red_ds.GetDriver().Create(os.path.join(self.ws, os.path.split(self.bands['red'])[1].split('_')[0] + '_RGB.tif'), red_ds.RasterXSize, red_ds.RasterYSize, 3, red_band.DataType, options = ['COMPRESS=LZW', 'TILED=YES', 'BIGTIFF=YES'] )
        out_ds.SetProjection(red_ds.GetProjection())
        out_ds.SetGeoTransform(red_ds.GetGeoTransform())
        out_ds.GetRasterBand(1).WriteArray(red_band.ReadAsArray())
        out_ds.GetRasterBand(2).WriteArray(green_ds.ReadAsArray())
        out_ds.GetRasterBand(3).WriteArray(blue_ds.ReadAsArray())
        out_ds.FlushCache()

        for i in range(1, 4):
            out_ds.GetRasterBand(i).ComputeStatistics(False)
        out_ds.BuildOverviews('average', [2, 4, 8, 16, 32])

        self.logger.info('Production of RGB raster of %s has been successfully fnished', os.path.split(self.bands['red'])[1].split('_')[0] + '_RGB.tif')

        del red_ds
        del green_ds
        del blue_ds
        del out_ds

    @staticmethod
    def convert_to_8bit(in_ras_file):
        ''' Convert a 16 bit raster to 8bit - returns raster'''

        ras_ds = gdal.Open(in_ras_file)
        ras_band = ras_ds.GetRasterBand(1)
        min_bl, max_bl = ras_band.ComputeRasterMinMax()
        mean_bl, sd_bl = ras_band.ComputeBandStats()

        band_val = [0, max(min_bl, mean_bl - 2*sd_bl), min(max_bl, mean_bl + 2*sd_bl), 65536]
        transf_val = [1, 2, 254, 255]
        transf_Func = interp1d(band_val, transf_val)

        out_ds = ras_ds.GetDriver().Create(os.path.join(os.path.split(in_ras_file)[0], os.path.split(in_ras_file)[1].split('.')[0] + '_8bit.tif'), ras_ds.RasterXSize, ras_ds.RasterYSize, 1, gdal.gdalconst.GDT_Byte)
        out_ds.SetGeoTransform(ras_ds.GetGeoTransform())
        out_ds.SetProjection(ras_ds.GetProjection())

        blue_np = ras_band.ReadAsArray()
        b_zero = np.zeros((ras_ds.RasterYSize, ras_ds.RasterXSize), dtype=int)
        for i in range(ras_ds.RasterYSize):
            b_zero[i] = transf_Func(blue_np[i]).astype(int)
        out_ds.GetRasterBand(1).WriteArray(b_zero)
        out_ds.FlushCache()

        Sentinel2_process.logger.info('Convertion to 8bit of %s has been successfully finished.', in_ras_file)

        del ras_ds
        del b_zero

        return out_ds

    @staticmethod
    def create_raster(in_ras, ras_pathname, np_data, data_type, nodata=None):
        '''Create a one-band GeoTIFF '''

        driver = gdal.GetDriverByName('GTiff')
        out_ds = driver.Create(ras_pathname, in_ras.RasterXSize, in_ras.RasterYSize, 1, data_type)
        out_ds.SetProjection(in_ras.GetProjection())
        out_ds.SetGeoTransform(in_ras.GetGeoTransform())
        out_band = out_ds.GetRasterBand(1)
        if nodata is not None:
            out_band.SetNoDataValue(nodata)
        out_band.WriteArray(np_data)
        out_band.FlushCache()
        out_band.ComputeStatistics(False)

        Sentinel2_process.logger.info('Creation of %s finished', ras_pathname)

        return out_ds

    def _get_bands(self):
        ''' Get the path for each band - returns a dictionary'''

        bands={}
        for b in glob(os.path.join(self.ws, '*.tif')):
            try:
                bnd = Sentinel2_process.prefixes[os.path.split(b)[-1].split('_')[1].split('.')[0]]
                # bnd = Sentinel2_process.prefixes[b.split('.')[0].split('_')[-2]]
                bands[bnd]=os.path.join(self.ws, b)
            except:
                print('No such band')
            #print bands[bnd]
        self.logger.info('Directory pathname for the specified bands of %s have been successfully initialized', self.ws)
        return bands


def get_first_day(in_day):
    '''Get the first day of the in_day's belonging month - returns string '''

    return date(in_day.year, in_day.month, 1).strftime('%Y-%m-%d')


def get_last_day(in_day):
    '''Get the last day of the in_day's belonging month - returns string '''

    return (date(in_day.year, (in_day + relativedelta(months=1)).month, 1) - timedelta(days=1)).strftime('%Y-%m-%d')


def get_s2_new_naming_convention(in_name):
    '''Get the NEW format naming convention of Sentinel-2 L1C products - returns list'''

    mission = in_name[0:3]
    product_level = in_name[4:10]
    sensing_timestamp = in_name[11:26]
    processing_baseline_number = in_name[27:32]
    orbit_number = in_name[33:37]
    tile_number_field = in_name[39:44]
    product_discriminator = in_name[45:59]

    return [mission, product_level, sensing_timestamp, processing_baseline_number, orbit_number, tile_number_field,
            product_discriminator]


def scihub_credentials():
    '''Get Scihub user, pass, url - returns tuple'''

    # get SciHub url instance, username and password
    sci_hub = Scihub_instance.objects.get(instance_desc__icontains='sci')
    sci_users = sci_hub.scihub_users_set.all()
    # TODO, each user for each country, make dict
    sci_user = sci_users[1].sci_user.encode('utf-8')
    sci_passwd = sci_users[1].sci_passwd.encode('utf-8')
    sci_url = sci_hub.instance_url.encode('utf-8')

    return (sci_user, sci_passwd, sci_url)


def create_vrt(out_path_vrt, *rasters):
    ''' Create a VRT file - returns nothing '''

    if len(rasters)==1 and hasattr(rasters[0], '__iter__'):
        rasters = rasters[0]

        in_ras = gdal.Open(rasters[0])
        columns, rows = in_ras.RasterXSize, in_ras.RasterYSize

        ras_vrt = gdal.GetDriverByName('VRT').Create(out_path_vrt, columns, rows, 0)
        ras_vrt.SetProjection(in_ras.GetProjection())
        ras_vrt.SetGeoTransform(in_ras.GetGeoTransform())

        del in_ras

        for idx, ras in enumerate(rasters):

            in_ras = gdal.Open(ras)
            in_ras_band = in_ras.GetRasterBand(1)

            ras_vrt.AddBand(in_ras_band.DataType)
            vrt_band = ras_vrt.GetRasterBand(idx + 1)

            band_meta = {'Band_%d' % (idx + 1): os.path.split(ras)[1]}
            vrt_band.SetMetadata(band_meta)

            xml = '''
                <SimpleSource>
                    <SourceFilename relativeToVRT="1">{ras}</SourceFilename>
                    <SourceBand>1</SourceBand>
                    <SrcRect xOff="0" yOff="0" xSize="{cols}" ySize="{rows}" />
                    <DstRect xOff="0" yOff="0" xSize="{cols}" ySize="{rows}" />
                </SimpleSource>
                '''
            data = {'ras': ras, 'cols': columns, 'rows': rows}
            meta = {'source_0': xml.format(**data)}
            vrt_band.SetMetadata(meta, 'vrt_sources')
            #print data
            MyLogger.logger.info('Raster %s just added', os.path.split(ras)[1])

        del in_ras
        del ras_vrt



def create_feature_space(in_aoi_id, in_path_name_vrt, out_path_name_csv):
    ''' Create feature space - returns nothing'''

    ras = gdal.Open(in_path_name_vrt)
    gt = ras.GetGeoTransform()
    inv_gt = gdal.InvGeoTransform(gt)

    with open(out_path_name_csv, 'wb') as csv_space:
        # writer = csv.writer(csv_space, quotechar='"', quoting=csv.QUOTE_ALL)
        writer = csv.writer(csv_space)
        csv_columns = ['id']

        for nband in range(ras.RasterCount):
            k_meta = 'Band_%d' % (nband + 1)
            band_name = ras.GetRasterBand(nband + 1).GetMetadata_Dict()[k_meta]
            band_name = band_name.split('.')[0]
            csv_columns.append(band_name)
        csv_columns.append('crop_code')
        csv_columns.append('crop_family')
        csv_columns.append('crop_season')

        writer.writerow(csv_columns)
        #cnt = 1
        for parcel in Parcels.objects.filter(aoi_id=in_aoi_id):
            parcel_data = [parcel.id]
            parcel_data = unicode(parcel_data, "utf-8")
            for nband in range(ras.RasterCount):
                try:
                    in_band = ras.GetRasterBand(nband + 1)
                    mesos_oros = parcel.zonal_stats(in_band, ras.GetProjection(), gt, inv_gt)
                    parcel_data.append(mesos_oros)
                except:
                    exc_type, exc_obj, exc_tb = sys.exc_info()
                    MyLogger.logger.error('Line: %d,\t Type: %s,\t Message: %s', exc_tb.tb_lineno, exc_type, exc_obj)
                    continue
            parcel_data.append(parcel.crop_code.id)
            #parcel_data.append(parcel.crop_code.crop_code)
            #print parcel.crop_code.crop_code
            parcel_data.append(parcel.crop_code.family_code)
            parcel_data.append(parcel.crop_code.seasonal_code)
            #parcel_data.append(parcel.famtype)
            #parcel_data.append(parcel.seastype)


            writer.writerow(parcel_data)



def create_feature_space_cloud(in_aoi_id, in_path_name_vrt, out_path_name_csv):
    ''' Create feature space - returns nothing'''

    ras = gdal.Open(in_path_name_vrt)
    gt = ras.GetGeoTransform()
    inv_gt = gdal.InvGeoTransform(gt)

    with open(out_path_name_csv, 'wb') as csv_space:
        # writer = csv.writer(csv_space, quotechar='"', quoting=csv.QUOTE_ALL)
        writer = csv.writer(csv_space)
        csv_columns = ['id']

        for nband in range(ras.RasterCount):
            k_meta = 'Band_%d' % (nband + 1)
            band_name = ras.GetRasterBand(nband + 1).GetMetadata_Dict()[k_meta]
            band_name = band_name.split('.')[0]
            csv_columns.append(band_name)
        #csv_columns.append('crop_code')
        #csv_columns.append('crop_family')
        #csv_columns.append('crop_season')

        writer.writerow(csv_columns)
        for parcel in Parcels.objects.filter(aoi_id=in_aoi_id):
            parcel_data = [parcel.id]

            for nband in range(ras.RasterCount):
                try:
                    in_band = ras.GetRasterBand(nband + 1)
                    res = parcel.zonal_stats(in_band, ras.GetProjection(), gt, inv_gt,stat_alg='majority')
                    athroisma = float(sum(res[0]))
                    if res[0]:
                        if len(res[0])>1:
                            zeros = res[0][0]
                            aces = res[0][1]
                            #print zeros,aces
                            #print "Perc:",aces/(zeros+aces)
                            if aces/(zeros+aces) > 0.5:
                                decision = 1
                            else:
                                decision = 0
                        else:
                            #print res[1]
                            decision = res[1]
                    else:
                        decision = 0
                    parcel_data.append(decision)
                except:
                    exc_type, exc_obj, exc_tb = sys.exc_info()
                    MyLogger.logger.error('Line: %d,\t Type: %s,\t Message: %s', exc_tb.tb_lineno, exc_type, exc_obj)
                    continue
            parcel_data.append(parcel.crop_code.id)
            #parcel_data.append(parcel.crop_code.crop_code)
            #print parcel.crop_code.crop_code
            #parcel_data.append(parcel.crop_code.family_code)
            #parcel_data.append(parcel.crop_code.seasonal_code)
            #parcel_data.append(parcel.famtype)
            #parcel_data.append(parcel.seastype)


            writer.writerow(parcel_data)


def create_feature_space_BAK(in_path_name_vrt, out_path_name_csv,aoi_id):
    ''' Create feature space - returns nothing'''

    #parcels = in_polygons_geoqueryset # Parcels.objects.filter(aoi=549)
    in_aoi_id = aoi_id
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(3857)

    ras = gdal.Open(in_path_name_vrt)
    gt = ras.GetGeoTransform()
    inv_gt = gdal.InvGeoTransform(gt)

    with open(out_path_name_csv, 'wb') as csv_space:
        # writer = csv.writer(csv_space, quotechar='"', quoting=csv.QUOTE_ALL)
        writer = csv.writer(csv_space)
        csv_columns = ['id']

        for nband in range(ras.RasterCount):
            k_meta = 'Band_%d' % (nband + 1)
            band_name = ras.GetRasterBand(nband + 1).GetMetadata_Dict()[k_meta]
            band_name = band_name.split('.')[0]
            csv_columns.append(band_name)
        csv_columns.append('crop_code')
        csv_columns.append('crop_family')
        csv_columns.append('crop_season')

        writer.writerow(csv_columns)

        for parcel in Parcels.objects.filter(aoi_id=in_aoi_id):
            parcel_data = [parcel.id]
            #thanos
            exists =1
            for nband in range(ras.RasterCount):
                try:
                    in_band = ras.GetRasterBand(nband + 1)

                    vect_tmp_drv = ogr.GetDriverByName('MEMORY')
                    vect_tmp_src = vect_tmp_drv.CreateDataSource('')
                    vect_tmp_lyr = vect_tmp_src.CreateLayer('', srs, ogr.wkbPolygon)
                    vect_tmp_lyr.CreateField(ogr.FieldDefn("id", ogr.OFTInteger))

                    feat = ogr.Feature(vect_tmp_lyr.GetLayerDefn())
                    feat.SetField("id", parcel.id)
                    feat_geom = ogr.CreateGeometryFromWkt(parcel.geom.wkt)
                    feat.SetGeometry(feat_geom)
                    vect_tmp_lyr.CreateFeature(feat)

                    xmin, xmax, ymin, ymax = feat_geom.GetEnvelope()

                    off_ulx, off_uly = map(int, gdal.ApplyGeoTransform(inv_gt, xmin, ymax))
                    off_lrx, off_lry = map(int, gdal.ApplyGeoTransform(inv_gt, xmax, ymin))
                    rows, columns = (off_lry - off_uly) + 1, (off_lrx - off_ulx) + 1
                    ras_tmp = gdal.GetDriverByName('MEM').Create('', columns, rows, 1, gdal.GDT_Byte)
                    ras_tmp.SetProjection(ras.GetProjection())
                    ras_gt = list(gt)
                    ras_gt[0], ras_gt[3] = gdal.ApplyGeoTransform(gt, off_ulx, off_uly)
                    ras_tmp.SetGeoTransform(ras_gt)

                    gdal.RasterizeLayer(ras_tmp, [1], vect_tmp_lyr, burn_values=[1])
                    data_band = in_band.ReadAsArray(off_ulx, off_uly, columns, rows).astype(np.float)
                    mask = ras_tmp.GetRasterBand(1).ReadAsArray()
                    zone_ras = np.ma.masked_array(data_band, np.logical_not(mask), fill_value=-9999)
                    zone_ras_list = zone_ras.compressed().tolist()
                    parcel_data.append(np.mean(zone_ras_list))

                except:
                    # thanos
                    exists = 0
                    exc_type, exc_obj, exc_tb = sys.exc_info()
                    MyLogger.logger.error('Line: %d,\t Type: %s,\t Message: %s', exc_tb.tb_lineno, exc_type, exc_obj)
                    continue

            #if exists==1:
            parcel_data.append(parcel.crop_code.id)
            parcel_data.append(parcel.crop_code.family_code)
            parcel_data.append(parcel.crop_code.seasonal_code)
            writer.writerow(parcel_data)


def create_shp(in_aoi_id, out_ws, out_shp_name='parcels.shp'):
    ''' Create a shapefile, consisted of land parcel polygons - return nothing'''

    srs = osr.SpatialReference()
    srs.ImportFromEPSG(3857)

    out_drv = ogr.GetDriverByName("ESRI Shapefile")
    out_ds = out_drv.CreateDataSource(os.path.join(out_ws, out_shp_name))
    out_lyr = out_ds.CreateLayer('', srs, ogr.wkbPolygon)
    out_lyr.CreateField(ogr.FieldDefn("id", ogr.OFTInteger))

    for parcel in Parcels.objects.filter(aoi_id=in_aoi_id):
        feat = ogr.Feature(out_lyr.GetLayerDefn())
        feat.SetField("id", parcel.id)
        feat.SetGeometry(ogr.CreateGeometryFromWkt(parcel.geom.wkt))

        out_lyr.CreateFeature(feat)

    out_ds.Destroy()


def check_bsm_running():
    ''' check if bsm is running '''

    for proc in psutil.process_iter():
        cmdline = proc.cmdline()
        if len(cmdline) > 1 and cmdline[1] == r'/home/userdev/PycharmProjects/BSM_Sentinel/scripts/s2_fire.py':
            return True


# den xrhsimopoiheitai kapou (25/7/2017)
def get_extent(in_ras_file):
    ''' Returns min_x, max_y, max_x, min_y'''

    in_ras = gdal.Open(in_ras_file)
    gt = in_ras.GetGeoTransform()
    return (gt[0], gt[3], gt[0] + gt[1] * in_ras.RasterXSize, gt[3] + gt[5] * in_ras.RasterYSize)

def get_month(string):

    m = {
        'jan': 1,
        'feb': 2,
        'mar': 3,
        'apr':4,
         'may':5,
         'jun':6,
         'jul':7,
         'aug':8,
         'sep':9,
         'oct':10,
         'nov':11,
         'dec':12
        }
    s = string.strip()[:3].lower()

    try:
        out = m[s]
        return out
    except:
        raise ValueError('Not a month')


def create_fs(in_path_name_vrt,aoi_id,out_csv):
    stat_algorithm = {
        'mean': np.mean,
        'max': np.max,
        'min': np.min,
        'majority': np.bincount
    }

    in_aoi_id = aoi_id
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(3857)

    ras = gdal.Open(in_path_name_vrt)
    gt = ras.GetGeoTransform()
    inv_gt = gdal.InvGeoTransform(gt)

    total = []
    first = []
    i = 0

    for band in range(ras.RasterCount):
        values = ras.GetRasterBand(band + 1).ReadAsArray().astype(np.float)
        i += 1
        k_meta = 'Band_%d' % (band + 1)
        band_name = ras.GetRasterBand(band + 1).GetMetadata_Dict()[k_meta]
        band_name = band_name.split('.')[0]
        classes = []
        cnt = 1
        col_list = []
        col_list.append(band_name)
        first.append("id")
        classes.append("Crop Type")
        family = []
        season = []
        family.append("Crop Family")
        season.append("Crop Season")
        for parcel in Parcels.objects.filter(aoi_id=in_aoi_id):
            id = parcel.id
            if i == 1:
                first.append(id)
            cropid = parcel.crop_code.id
            #parcel_data.append(parcel.crop_code.id)
            #parcel_data.append(parcel.crop_code.family_code)
            #parcel_data.append(parcel.crop_code.seasonal_code)
            classes.append(cropid)
            family.append(parcel.crop_code.family_code)
            season.append(parcel.crop_code.seasonal_code)
            parcel_data = [id]
            try:
                vect_tmp_drv = ogr.GetDriverByName('MEMORY')
                vect_tmp_src = vect_tmp_drv.CreateDataSource('')
                vect_tmp_lyr = vect_tmp_src.CreateLayer('', srs, ogr.wkbPolygon)
                vect_tmp_lyr.CreateField(ogr.FieldDefn("id", ogr.OFTInteger))

                feat = ogr.Feature(vect_tmp_lyr.GetLayerDefn())
                feat.SetField("id", parcel.id)
                feat_geom = ogr.CreateGeometryFromWkt(parcel.geom.wkt)
                feat.SetGeometry(feat_geom)
                vect_tmp_lyr.CreateFeature(feat)

                xmin, xmax, ymin, ymax = feat_geom.GetEnvelope()

                off_ulx, off_uly = map(int, gdal.ApplyGeoTransform(inv_gt, xmin, ymax))
                off_lrx, off_lry = map(int, gdal.ApplyGeoTransform(inv_gt, xmax, ymin))
                rows, columns = (off_lry - off_uly) + 1, (off_lrx - off_ulx) + 1
                ras_tmp = gdal.GetDriverByName('MEM').Create('', columns, rows, 1, gdal.GDT_Byte)
                ras_tmp.SetProjection(ras.GetProjection())
                ras_gt = list(gt)
                ras_gt[0], ras_gt[3] = gdal.ApplyGeoTransform(gt, off_ulx, off_uly)
                ras_tmp.SetGeoTransform(ras_gt)

                gdal.RasterizeLayer(ras_tmp, [1], vect_tmp_lyr, burn_values=[1])
                mask = ras_tmp.GetRasterBand(1).ReadAsArray()
                aa = off_uly
                bb = off_lry + 1
                cc = off_ulx
                dd = off_lrx + 1
                zone_ras = np.ma.masked_array(values[aa:bb, cc:dd], np.logical_not(mask), fill_value=-9999)
                zone_ras_list = zone_ras.compressed().tolist()
                # parcel_data.append(np.mean(zone_ras_list))
                #parcel_data.append(cropid)
                # writer.writerow(parcel_data)
                # f[cnt][i] = str(np.mean(zone_ras_list))
                # cnt+=1
                col_list.append(np.mean(zone_ras_list))
                #maj = stat_algorithm['mean'](zone_ras_list)
                #res = maj.tolist(), maj.argmax()
            except:
                col_list.append("Out")
                cnt += 1
                exc_type, exc_obj, exc_tb = sys.exc_info()
                # logger.info('Line: %d,\t Type: %s,\t Message: %s', exc_tb.tb_lineno, exc_type, exc_obj)
                continue


                # if exists==1:
                # parcel_data.append(family)
                # parcel_data.append(season)
        if i == 1:
            total.append(first)
        total.append(col_list)

        del values


    total.append(classes)
    total.append(family)
    total.append(season)
    rows = zip(*total)
    with open(out_csv, 'w') as f:
        writer = csv.writer(f);
        for row in rows:
            writer.writerow(row)


def create_cloud(in_path_name_vrt,aoi_id,out_csv):
    stat_algorithm = {
        'mean': np.mean,
        'max': np.max,
        'min': np.min,
        'majority': np.bincount
    }

    in_aoi_id = aoi_id
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(3857)

    ras = gdal.Open(in_path_name_vrt)
    gt = ras.GetGeoTransform()
    inv_gt = gdal.InvGeoTransform(gt)

    total = []
    first = []
    i = 0

    for band in range(ras.RasterCount):
        values = ras.GetRasterBand(band + 1).ReadAsArray().astype(np.float)
        k_meta = 'Band_%d' % (band + 1)
        band_name = ras.GetRasterBand(band + 1).GetMetadata_Dict()[k_meta]
        band_name = band_name.split('.')[0]
        classes = []
        cnt = 1
        col_list = []
        col_list.append(band_name)
        first.append("id")
        classes.append("Crop Type")
        family = []
        season = []
        family.append("Crop Family")
        season.append("Crop Season")
        for parcel in Parcels.objects.filter(aoi_id=in_aoi_id):
            id = parcel.id
            if i == 1:
                first.append(id)
            cropid = parcel.crop_code.id
            classes.append(cropid)
            family.append(parcel.crop_code.family_code)
            season.append(parcel.crop_code.seasonal_code)
            parcel_data = [id]
            try:
                vect_tmp_drv = ogr.GetDriverByName('MEMORY')
                vect_tmp_src = vect_tmp_drv.CreateDataSource('')
                vect_tmp_lyr = vect_tmp_src.CreateLayer('', srs, ogr.wkbPolygon)
                vect_tmp_lyr.CreateField(ogr.FieldDefn("id", ogr.OFTInteger))

                feat = ogr.Feature(vect_tmp_lyr.GetLayerDefn())
                feat.SetField("id", parcel.id)
                feat_geom = ogr.CreateGeometryFromWkt(parcel.geom.wkt)
                feat.SetGeometry(feat_geom)
                vect_tmp_lyr.CreateFeature(feat)

                xmin, xmax, ymin, ymax = feat_geom.GetEnvelope()

                off_ulx, off_uly = map(int, gdal.ApplyGeoTransform(inv_gt, xmin, ymax))
                off_lrx, off_lry = map(int, gdal.ApplyGeoTransform(inv_gt, xmax, ymin))
                rows, columns = (off_lry - off_uly) + 1, (off_lrx - off_ulx) + 1
                ras_tmp = gdal.GetDriverByName('MEM').Create('', columns, rows, 1, gdal.GDT_Byte)
                ras_tmp.SetProjection(ras.GetProjection())
                ras_gt = list(gt)
                ras_gt[0], ras_gt[3] = gdal.ApplyGeoTransform(gt, off_ulx, off_uly)
                ras_tmp.SetGeoTransform(ras_gt)

                gdal.RasterizeLayer(ras_tmp, [1], vect_tmp_lyr, burn_values=[1])
                mask = ras_tmp.GetRasterBand(1).ReadAsArray()
                aa = off_uly
                bb = off_lry + 1
                cc = off_ulx
                dd = off_lrx + 1
                zone_ras = np.ma.masked_array(values[aa:bb, cc:dd], np.logical_not(mask), fill_value=-9999)
                zone_ras_list = zone_ras.compressed().tolist()
                # parcel_data.append(np.mean(zone_ras_list))
                #parcel_data.append(cropid)
                # writer.writerow(parcel_data)
                # f[cnt][i] = str(np.mean(zone_ras_list))
                # cnt+=1
                maj = stat_algorithm['majority'](zone_ras_list)
                res = maj.tolist(), maj.argmax()
                athroisma = float(sum(res[0]))
                if res[0]:
                    if len(res[0]) > 1:
                        zeros = res[0][0]
                        aces = res[0][1]
                        if aces / (zeros + aces) > 0.5:
                            decision = 1
                        else:
                            decision = 0
                    else:
                        decision = res[1]
                else:
                    decision = 0
                col_list.append(decision)
            except:
                col_list.append(-9999)
                cnt += 1
                exc_type, exc_obj, exc_tb = sys.exc_info()
        if i == 1:
            total.append(first)
            for k in range(0, 14):
                total.append(col_list)
        else:
            for k in range(0,14):
                total.append(col_list)
        del values


    total.append(classes)
    total.append(family)
    total.append(season)
    rows = zip(*total)
    with open(out_csv, 'w') as f:
        writer = csv.writer(f);
        for row in rows:
            writer.writerow(row)





def get_token(in_url='https://app.recap-h2020.eu/backend/api/authentication/login/', in_data={"username": "alelentzis", "email": "alelentzis@noa.gr", "password": "r3c@pNoa"}):
    ''' Get token from RECAP system '''

    while True:
        req_token = Make_request.create_token(url=in_url, data=in_data)
        if req_token['status']==200:
            token = req_token['response']['token'].encode('utf-8')
            token_header = {'Authorization': 'Bearer %s' %token,}
            MyLogger.logger.info('Json Web Token (%s) has been succeesfully acquired', token)

            return token_header

        else:
            MyLogger.logger.warning('Bad response for token request: %s', req_token['response'])
            time.sleep(300)
