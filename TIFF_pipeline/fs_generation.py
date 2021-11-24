
# Feature Space Generation

from osgeo import gdal
from osgeo import ogr
from osgeo import osr
from settings import SETTINGS
from utils import *
import numpy as np
import zipfile
import fnmatch
import csv
import pdb
import sys
import os
import re

# This section needs to be rewritten
# tile = SETTINGS['tiles'][0]
outname = SETTINGS['fs_filename']
id_colname = "id"

# Doesn't make sense when its getcwd but it may change (most probably)
# ws_tmp = os.getcwd()
workdir = SETTINGS['working_directory']
# os.chdir(ws_tmp)

#### Shapefile should be in 3857 projection!! #####
shapefile = SETTINGS['shapefile']

order = {
    'BLUE':  0,
    'GREEN': 1,
    'RED':   2,
    'NIR':   3,
}
number_of_bands = len(order)

stat_algorithm = {
    'mean': np.mean,
    'max': np.max,
    'min': np.min,
    'majority': np.bincount
}

bands = []
orderedbands = []
band_cloudmasks = {}
acquisitions = {}
for root, dirs, files in sorted(os.walk(workdir)):
    for tif in fnmatch.filter(files, '*.tif'):
        if any(x in tif for x in ['Jan', 'Feb']):
            continue
        if re.match('(.*)(BLUE|GREEN|RED|NIR)(.*)(intersection)(.*)', tif):
            bands.append(os.path.join(root, tif))
            orderedbands.append('')
            print ('Matched: {}'.format(tif))
            # Split different folders to process each acquisition separately
            acq = '_'.join(tif.split('_')[:5])
            if acq not in acquisitions: # create the acquisition key if its the first band we encounter
                acquisitions[acq] = {}
                acquisitions[acq]['bands'] = []
                acquisitions[acq]['ordered_bands'] = []
            acquisitions[acq]['bands'].append(os.path.join(root,tif)) # and append it
            acquisitions[acq]['ordered_bands'].append('')
        # Now check if it is a TRIPLESAT acquisition. We will use this info
        # to get the proper cloudmask path
        if re.match('TRIPLESAT(.*)',tif):
            print('  It is a TRIPLESAT acquisition')
            # We know where to find the cloudmask in this case
            cloud_dir = '/'.join(root.split('/')[:-1]) + '/Cloudmask'
            cloud_file = 'cloudmask.tif'
            cloudmask = os.path.join(cloud_dir,cloud_file)
            # Map band to cloudmask
            band_cloudmasks[tif] = cloudmask
            print("  Corresponding cloudmask directory: {0}".format(cloud_dir))
            print("  Filename: {0}".format(cloud_file))

# Let's try to perform different runs for different acquisition tiles
print("Bands:")
for acq in acquisitions.items():
    acq_key = acq[0]
    acq_dict = acq[1]
    for item in acq_dict['bands']:
        print(item)

iteration = 0
for acq in acquisitions.items():
    iteration += 1
    acq_key = acq[0]
    acq_dict = acq[1]
    k = 1
    x = 0
    for item in acquisitions[acq_key]['bands']:
        # For the Sentinel scripts this referred to the preprocessed bands as they appear on the
        # preprocessed filname, ie. B02, B03, B04, B05, B06, B07, B08, B8A, B11, B12
        # key = item.split('.')[0].split('/')[-1].split('_')[1]

        # We just keep the filename without the extension, and then extract the (BLUE|GREEN|RED|NIR)
        # We use [-3] in the last split because that is the 3rd from the end (as the filename
        # ends -- curently -- with eg. RED_3857_intersection)
        key = item.split('.')[0].split('/')[-1].split('_')[-3]
        acquisitions[acq_key]['ordered_bands'][order[key]+x*number_of_bands] = item
        if k%number_of_bands==0:
            x+=1
        k+=1


    ########################
    # VRT Creation attempt #
    ########################
    # based on the same function for the S2 scripts
    ''' Create a VRT file - returns nothing '''
    print("Creating vrt file of VHR satellite bands ...")
    create_vrt(os.path.join(workdir,'stack_{}.vrt'.format(iteration)), acquisitions[acq_key]['ordered_bands'])
    in_path_name_vrt = os.path.join(workdir,'stack_{}.vrt'.format(iteration))

    # pdb.set_trace()

    # cloudbands = list(band_cloudmasks.values())
    cloudband = band_cloudmasks[acquisitions[acq_key]['bands'][0].split('/')[-1]]
    # # This loop should happen in the Cloudmask folder
    # for root, dirs, files in sorted(os.walk(os.path.join(workdir,'../Cloudmask/'))):
    #     for tif in fnmatch.filter(files, '*.tif'):
    #         if any(x in tif for x in ['Jan', 'Feb']):
    #             continue
    #         if re.match("(.*)([Cc]loudmask)(.*)", tif):
    #             for i in range(0,number_of_bands):
    #                 bands.append(os.path.join(root, tif))

    print("Creating vrt file for cloudmasks ...")
    create_vrt(os.path.join(workdir,'cloudstack_{}.vrt'.format(iteration)), [cloudband])
    in_path_name_vrt2 = os.path.join(workdir,'cloudstack_{}.vrt'.format(iteration))

    # pdb.set_trace()

    srs = osr.SpatialReference()
    srs.ImportFromEPSG(3857)


    ras = gdal.Open(in_path_name_vrt)
    gt = ras.GetGeoTransform()
    inv_gt = gdal.InvGeoTransform(gt)


    ras2 = gdal.Open(in_path_name_vrt2)
    gt2 = ras2.GetGeoTransform()
    inv_gt2 = gdal.InvGeoTransform(gt2)

    # This was previously in the band loop, but I don't get why
    clouded = ras2.GetRasterBand(1).ReadAsArray().astype(np.float32)

    total = []
    first = []
    second = []
    i = 0


    out_csv = os.path.join(workdir, outname)
    # TEMP
    out_csv = os.path.join(workdir, 'VHR_fs_dandrimont_{}.csv'.format(iteration))

    for band in range(ras.RasterCount):
        k_meta = 'Band_%d' % (band + 1)
        band_full_name = ras.GetRasterBand(band + 1).GetMetadata_Dict()[k_meta]
        band_name = band_full_name.split('_')[4] + '_' + band_full_name.split('_')[9]
        print('Processing band: {0}'.format(band_name))

        values = ras.GetRasterBand(band + 1).ReadAsArray().astype(np.float32)
        result = np.where(clouded==0,values,np.nan)

        i += 1

        # band_name2 = ras2.GetRasterBand(band + 1).GetMetadata_Dict()[k_meta]
        # # print (band_name2)
        # band_name2 = band_name2.split('.')[0]

        cnt = 1
        col_list = []
        col_list.append(band_name)
        first.append("id")
        #second.append("code")
        driver = ogr.GetDriverByName("ESRI Shapefile")
        dataSource = driver.Open(shapefile, 0)
        layer = dataSource.GetLayer()

        print('Starting feature loop')
        for feature in layer:
            geom = feature.GetGeometryRef()
            geom = geom.ExportToWkt()

            id = feature.GetField(id_colname)
            if i == 1:
                first.append(id)
            # parcel_data = [id]
            try:
                vect_tmp_drv = ogr.GetDriverByName('MEMORY')
                vect_tmp_src = vect_tmp_drv.CreateDataSource('')
                vect_tmp_lyr = vect_tmp_src.CreateLayer('', srs, ogr.wkbPolygon)
                vect_tmp_lyr.CreateField(ogr.FieldDefn("id", ogr.OFTInteger))

                feat = ogr.Feature(vect_tmp_lyr.GetLayerDefn())
                feat.SetField("id", id)
                feat_geom = ogr.CreateGeometryFromWkt(geom)
                feat.SetGeometry(feat_geom)
                vect_tmp_lyr.CreateFeature(feat)

                xmin, xmax, ymin, ymax = feat_geom.GetEnvelope()

                off_ulx, off_uly = map(int, gdal.ApplyGeoTransform(inv_gt2, xmin, ymax))
                off_lrx, off_lry = map(int, gdal.ApplyGeoTransform(inv_gt2, xmax, ymin))
                rows, columns = off_uly - off_lry + 1, off_lrx - off_ulx + 1

                ras_tmp = gdal.GetDriverByName('MEM').Create('', columns, rows, 1, gdal.GDT_Byte)
                ras_tmp.SetProjection(ras2.GetProjection())
                ras_gt = list(gt2)
                ras_gt[0], ras_gt[3] = gdal.ApplyGeoTransform(gt, off_ulx, off_uly)
                ras_tmp.SetGeoTransform(ras_gt)

                gdal.RasterizeLayer(ras_tmp, [1], vect_tmp_lyr, burn_values=[1])
                mask = ras_tmp.GetRasterBand(1).ReadAsArray()
                # aa = off_uly
                # bb = off_lry + 1
                # cc = off_ulx
                # dd = off_lrx + 1
                ### GC Change
                aa = off_ulx
                bb = off_lrx + 1
                cc = off_lry
                dd = off_uly + 1

                zone_ras = np.ma.masked_array(clouded[aa:bb, cc:dd], np.logical_not(mask), fill_value=np.nan)
                zone_ras_list = zone_ras.compressed().tolist()

                maj = stat_algorithm['majority'](zone_ras_list)
                res = maj.tolist(), 2

                athroisma = float(sum(res[0]))
                #maj_value = res[1]
                if len(res[0])>1:
                    perc = float(res[0][1]/athroisma)
                else:
                    perc = 0

                vect_tmp_drv2 = ogr.GetDriverByName('MEMORY')
                vect_tmp_src2 = vect_tmp_drv2.CreateDataSource('')
                vect_tmp_lyr2 = vect_tmp_src2.CreateLayer('', srs, ogr.wkbPolygon)
                vect_tmp_lyr2.CreateField(ogr.FieldDefn("id", ogr.OFTInteger))

                feat2 = ogr.Feature(vect_tmp_lyr2.GetLayerDefn())
                feat2.SetField("id", id)
                feat_geom2 = ogr.CreateGeometryFromWkt(geom)
                feat2.SetGeometry(feat_geom2)
                vect_tmp_lyr2.CreateFeature(feat2)

                xmin, xmax, ymin, ymax = feat_geom2.GetEnvelope()

                # Why weren't these updated below?
                off_ulx2, off_uly2 = map(int, gdal.ApplyGeoTransform(inv_gt2, xmin, ymax))
                off_lrx2, off_lry2 = map(int, gdal.ApplyGeoTransform(inv_gt2, xmax, ymin))
                # rows2, columns2 = off_lry2 - off_uly2 + 1, off_lrx2 - off_ulx2 + 1
                rows2, columns2 = off_uly2 - off_lry2 + 1, off_lrx2 - off_ulx2 + 1
                ras_tmp = gdal.GetDriverByName('MEM').Create('', columns2, rows2, 1, gdal.GDT_Byte)
                ras_tmp.SetProjection(ras.GetProjection())
                ras_gt = list(gt)
                ras_gt[0], ras_gt[3] = gdal.ApplyGeoTransform(gt, off_ulx2, off_uly2)
                ras_tmp.SetGeoTransform(ras_gt)

                gdal.RasterizeLayer(ras_tmp, [1], vect_tmp_lyr, burn_values=[1])
                mask = ras_tmp.GetRasterBand(1).ReadAsArray()
                # aa = off_uly
                # bb = off_lry + 1
                # cc = off_ulx
                # dd = off_lrx + 1
                ### GC Change
                aa = off_ulx2
                bb = off_lrx2 + 1
                cc = off_lry2
                dd = off_uly2 + 1

                zone_ras = np.ma.masked_array(result[aa:bb, cc:dd], np.logical_not(mask), fill_value=np.nan)
                zone_ras_init = np.ma.masked_array(values[aa:bb, cc:dd], np.logical_not(mask), fill_value=np.nan)

                zone_ras_list = zone_ras_init.compressed().tolist()
                zone_ras_list2 = zone_ras.compressed()#.tolist()

                if perc > 0.8:
                    col_list.append(np.nan)
                else:
                    if np.all((zone_ras_list2 == np.nan)):
                        x = np.nan
                    else:
                        x = np.nanmean(zone_ras_list2)
                        # if ('PSRI' in band_name and (x<-5000 or x>5000)) or np.isnan(x):
                        #     x = -999999
                        if ('PSRI' in band_name and x==1):
                            x = np.nan
                        elif (not ('INDEX' in band_name)) and x==0:
                            x = np.nan

                    col_list.append(x)
                del zone_ras_list
                del zone_ras_list2
                del zone_ras
            except:
                # print('exception occurred for id {}'.format(id))
                col_list.append(np.nan)
                cnt += 1
                exc_type, exc_obj, exc_tb = sys.exc_info()
                continue

        if i == 1:
            total.append(first)
        total.append(col_list)

        del values
        # del clouded
        del result

    rows = zip(*total)
    with open(out_csv, 'w') as f:
        writer = csv.writer(f)
        for row in rows:
            # print row
            writer.writerow(row)
