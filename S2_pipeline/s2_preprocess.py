import zipfile
import shutil
from osgeo import gdal
from utils import *
from settings import SETTINGS, CREDENTIALS


tiles     = SETTINGS['tiles']
sudo_pass = CREDENTIALS['sudo_pass_bsc_vm']
cnt       = len(tiles) # number of tiles

def basic_preprocess(ws_tmp):

    months = SETTINGS['months']
    bands  = SETTINGS['bands']

    clip = 0

    #Lithuania
    # xmin,ymin,xmax,ymax =

    '''
    if cnt > 1:
        logger.info('The specified AOI (id=%d) is contained by more than one tile' % aoi.id)

        ossim_cmd = 'ossim-orthoigen -P %)s --ossim-logfile %s --hist-auto-minmax --combiner-type ossimFeatherMosaic --srs EPSG:3857 %s %s.tif'

    else:
        logger.info('It is not required the mosaicing process for the specified AOI (id=%d)' % aoi.id)

        ossim_cmd = 'ossim-orthoigen -P %s --ossim-logfile %s --srs EPSG:3857 %s %s.tif'
    '''
    for mon in os.listdir(ws_tmp):
        if '.' in mon:
            continue
        print("Processing month:",mon)
        mon_path = os.path.join(ws_tmp, mon)
        # transforming (converting and resampling) each band (10m and 20m) of SAFE product to geotiff
        for safe in glob(os.path.join(mon_path, '*.SAFE')):
            p = Sentinel2_pre_process(mon_path, safe)
            p.converting_rasters()
            p.delete_jp2_files()
        logging.info('Transforming (converting and resampling of basic rasters) process has been successfully finished')
        # if cnt == 1:
        #     for b in bands:
        #         for safe in glob(os.path.join(mon_path, '*.SAFE')):
        #             rasterdate = safe[-14:-12] + safe[-20:-16]
        #             fulldate = safe[-20:-12]
        #             for root, dirs, files in os.walk(mon_path):
        #                 if fulldate in root:
        #                     for tif in fnmatch.filter(files, '*' + b):
        #                         destination = os.path.join(mon_path, months[mon].split('_')[0] + rasterdate + '_all' + '_' +b.split('_')[0] + '.tif')
        #                         shutil.move(os.path.join(root, tif),destination)
        if cnt == 1:
            for b in bands:
                print(b)
                print((os.path.join(mon_path, '*.SAFE')))
                for safe in glob(os.path.join(mon_path, '*.SAFE')):
                    #   OLD      day of month           year
                    old_rasterdate = safe[-14:-12] + safe[-20:-16]
                    #   NEW      day of month           year
                    rasterdate = safe[-48:-46] + safe[-54:-50]
                    #   OLD
                    old_fulldate = safe[-20:-12]
                    #   NEW
                    fulldate = safe[-54:-46]
                    for root, dirs, files in os.walk(mon_path):
                        if fulldate in root:
                            for tif in fnmatch.filter(files, '*' + b):
                                # print os.path.join(root,tif)
                                gdal.Warp(os.path.join(mon_path, months[mon].split('_')[0] + rasterdate + '_all' + '_' +b.split('_')[0] + '.tif'), os.path.join(root, tif),
                                          dstSRS='EPSG:3857',xRes=10,yRes=10)
            logging.info('Reproject with GDAL and SRC files of %s files has been successfully generated',
                        months[mon].split('_')[0])

        # delete L2A files
        for safe in glob(os.path.join(mon_path, '*.SAFE')):
            shutil.rmtree(safe)
        logging.info('L2A files removing process has been successfully finished')

        if clip==1:
            # clip mosaic images
            for tif in glob(os.path.join(mon_path, '*.tif')):
                if '_TCI' in tif:
                    x = tif.split("_all_")[0]
                    xx = x + "_RGB.tif"
                    cmd = "gdalwarp -te %f %f %f %f %s %s" % (xmin, ymin, xmax, ymax, tif, xx)
                    os.system(cmd)
                    os.remove(tif)
                elif not "_RGB" in tif:
                    out_tif = os.path.split(tif)[1].split('_all')[0] + os.path.split(tif)[1].split('_all')[1]
                    # print out_tif
                    Sentinel2_pre_process.clip_raster(os.path.dirname(tif), tif, out_tif, xmin=xmin, ymin=ymin,
                                                      xmax=xmax,
                                                      ymax=ymax)
                    os.remove(tif)
            logging.info('Clipping process has been successfully finished')

        else:
            #without clip process
            for tif in glob(os.path.join(mon_path, '*.tif')):
                if '_TCI' in tif:
                    x = tif.split("_all_TCI")[0]
                    xx = x + "_RGB.tif"
                    os.rename(tif, xx)
                    # print x, xx
                elif 'all' in tif:
                    # print tif
                    x = tif.split('_all')
                    xx = x[0] + x[1]
                    # print xx
                    os.rename(tif, xx)

        dates = []
        for tif in glob(os.path.join(mon_path, '*.tif')):
            rasterdate = tif.split('.')[0].split('/')[-1].split('_')[0][-6:-4]
            if rasterdate not in dates:
                dates.append(rasterdate)

        for date in dates:

            os.mkdir(os.path.join(mon_path, date))
            for tif in glob(os.path.join(mon_path, '*.tif')):
                rasterdate = tif.split('.')[0].split('/')[-1].split('_')[0][-6:-4]
                if rasterdate == date:
                    cmd = 'echo %s | sudo -S mv %s %s' % (sudo_pass, tif, os.path.join(mon_path, date))
                    os.system(cmd)
        for date in dates:
            mon_path2 = os.path.join(mon_path, date)
            s2 = Sentinel2_process(mon_path2)
            s2.maskupdated()
            s2.create_ndvi()
            s2.create_psri()
            s2.create_ndwi_lth()
            s2.create_savi()
            s2.create_nbr_lth()
            # s2.create_bori()
            # s2.create_bari()
            for tif_8bit in glob(os.path.join(mon_path2, '*8bit.tif')):
                    os.remove(tif_8bit)
            logging.info('Removing process of intermediate files (8bit rasters) has been successfully finished')

for tile in tiles:
    ws_tmp = '{0}/2017/{1}'.format(SETTINGS['outdir'],tile)
    for month in os.listdir(ws_tmp):
        if "." in month:
            continue
        print("PROCESSING MONTH {0}".format(month))
        month_path = os.path.join(ws_tmp,month)
        os.chdir(month_path)
        for item in os.listdir(month_path):  # loop through items in dir
            if item.endswith('.zip'):  # check for ".zip" extension
                file_name = os.path.abspath(item)  # get full path of files
                with zipfile.ZipFile(file_name, 'r') as zip_ref:
                    zip_ref.extractall(month_path)  # extract file to dir
                    zip_ref.close()  # close file
                    os.remove(file_name)  # delete zipped file

        # os.system("unzip \*.zip")
        # for item in os.listdir(month_path):
        #     if item.endswith('.zip'):
        #         os.remove(os.path.join(month_path,item))

    basic_preprocess(ws_tmp)

