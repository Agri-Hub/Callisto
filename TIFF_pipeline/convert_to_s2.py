
from osgeo import gdal, osr
from utils import *
from settings import SETTINGS
import numpy as np
import os
import re
import pdb

# Search for files to process in the given directory
workdir = SETTINGS['working_directory']

# Will hold tuples of (acquisition, cloudmask)
# If cloudmask doesn't exist we can set it to the empty string ''
FILES_TO_PROCESS = []

# Walk the current working directory to find acquisitions and cloudmasks
for root, dirs, files in os.walk(workdir):
    for file in files:
        # Searching for TRIPLESAT acquisitions (Dutch Orthophotos)
        if re.match(r'^TRIPLESAT.*MS.*[0-9]+\.tif$',file):
            print("Found TRIPLESAT multispectral acquisition: {0}".format(file))
            print(" Directory: {0}".format(root))
            # We know where to find the cloudmask in this case
            cloud_dir = '/'.join(root.split('/')[:-1]) + '/Cloudmask'
            cloud_file = ''.join(file.split('MS_'))
            cloud_file = '.'.join(cloud_file.split('.')[:-1]) + '.shp'
            cloudmask = os.path.join(cloud_dir,cloud_file)
            print("  Corresponding cloudmask directory: {0}".format(cloud_dir))
            print("  Filename: {0}".format(cloud_file))
            FILES_TO_PROCESS.append((os.path.join(root,file),cloudmask))

for input, cloudmask in FILES_TO_PROCESS:
    print('Processing: {0}'.format(input.split('/')[-1]))

    # Convert 1 multiband tiff to N single-band tiffs
    band_names = ['BLUE', 'GREEN', 'RED', 'NIR']
    VHR_native_band_filenames = {}
    VHR_repr_filenames = {}
    VI_filenames = {}

    file_prefix = input.split('.')[0]
    file_suffix = input.split('.')[-1]
    # Projection to use
    s2_proj = 'EPSG:3857'

    data = gdal.Open(input)
    print('Number of bands: {}'.format(data.RasterCount))
    # Extract the bands to single image tiffs
    for band in range(data.RasterCount):
        print('Extracting band {} ...'.format(band+1))
        out_file = '{0}_{1}.{2}'.format(
                                    file_prefix,
                                    band_names[band],
                                    file_suffix)
        repr_out_file = '{0}_{1}_{2}.{3}'.format(
                                    file_prefix,
                                    band_names[band],
                                    s2_proj.split(':')[1],
                                    file_suffix)
        os.system(
            'gdal_translate {0} -b {1} {2}'
                .format(input, band+1, out_file)
        )

        # This is a dictionary that is being built which has the band names as keys
        # and the single-band filenames as values. Thus, it will look like this in the end:
        # VHR_native_band_filenames = {
        #     # before reprojecttion
        #     'BLUE':  'TRIPLESAT_3_MS_L1_20170923094819_001237VI_008_0220170831001001_029_BLUE.tif',
        #     'GREEN': 'TRIPLESAT_3_MS_L1_20170923094819_001237VI_008_0220170831001001_029_GREEN.tif',
        #     'RED':   'TRIPLESAT_3_MS_L1_20170923094819_001237VI_008_0220170831001001_029_RED.tif',
        #     'NIR':   'TRIPLESAT_3_MS_L1_20170923094819_001237VI_008_0220170831001001_029_NIR.tif',
        # }
        VHR_native_band_filenames[band_names[band]] = out_file

        # This is the dictionary that will hold the single-band filenames after reprojection
        # It will look like this in the end:
        # VHR_repr_filenames = {
        #     # after reprojecttion
        #     'BLUE':  'TRIPLESAT_3_MS_L1_20170923094819_001237VI_008_0220170831001001_029_BLUE_3857.tif',
        #     'GREEN': 'TRIPLESAT_3_MS_L1_20170923094819_001237VI_008_0220170831001001_029_GREEN_3857.tif',
        #     'RED':   'TRIPLESAT_3_MS_L1_20170923094819_001237VI_008_0220170831001001_029_RED_3857.tif',
        #     'NIR':   'TRIPLESAT_3_MS_L1_20170923094819_001237VI_008_0220170831001001_029_NIR_3857.tif',
        # }
        VHR_repr_filenames[band_names[band]] = repr_out_file

    for vi in ['NDVI', 'SAVI']:
        vi_out_file = '{0}_{1}_{2}.{3}'.format(
                        file_prefix,
                        vi,
                        s2_proj.split(':')[1],
                        file_suffix)
        # And this is the dictionary that will hold the filenames with info on specific vegetation indices
        # Once populated it should look like this:
        # VI_filenames = {
        #     # Vegetation Indices (after reprojection)
        #     'NDVI': 'TRIPLESAT_3_MS_L1_20170923094819_001237VI_008_0220170831001001_029_NDVI_3857.tif',
        #     'SAVI': 'TRIPLESAT_3_MS_L1_20170923094819_001237VI_008_0220170831001001_029_SAVI_3857.tif'
        # }
        VI_filenames[vi] = vi_out_file


    ###############################
    # Reproject single band tiffs #
    ###############################

    srs = osr.SpatialReference()
    srs.ImportFromEPSG(3857)

    for band, file in VHR_native_band_filenames.items():
        print('Reprojection of {0} ...'.format(file))

        file_prefix = file.split('.')[0]
        file_suffix = file.split('.')[-1]
        proj_name_part = s2_proj.split(':')[1]

        # Running through a GDAL system call
        # os.system('gdalwarp -t_srs {0} -r bilinear -tr 10 10 {1} {2}_{3}.{4}'.format(s2_proj, f_input, file_prefix, proj_name_part, file_suffix))

        # Running through the GDAL python library
        gdal.Warp(
            # Output file will be the same as the previous but with the projection at the end
            # '{0}_{1}.{2}'.format(file_prefix, proj_name_part, file_suffix),
            VHR_repr_filenames[band],
            file,
            xRes = 10,
            yRes = 10,
            dstSRS=s2_proj,
            resampleAlg='bilinear'
        )

    ###################
    # NDVI Generation #
    ###################
    ''' Returns Normalized Difference Vegetation Index '''

    print("NDVI Calculation initiated ...")

    red_ds = gdal.Open(VHR_repr_filenames['RED'])
    nir_ds = gdal.Open(VHR_repr_filenames['NIR'])

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

    out_ndvi = create_raster(gdal.Open(VHR_repr_filenames['RED']), os.path.join(os.getcwd(), VI_filenames['NDVI']), out_data, gdal.GDT_Int16, -9999)
    del out_data
    print("Finished NDVI calculation - TIFF file created")

    ###################
    # SAVI Generation #
    ###################
    # Most of the below can be removed and performed jointly on the section above
    # I probably only need to add the savi calculation
    # Check it later
    ''' Returns Soil Adjusted Vegetation Index '''

    print("SAVI Calculation initiated ...")

    # We already have these 2
    # red_ds = gdal.Open(self.bands['red'])
    # nir_ds = gdal.Open(self.bands['nir'])

    # and these
    # red_band = red_ds.GetRasterBand(1)
    # nir_band = nir_ds.GetRasterBand(1)

    # xsize, ysize = red_band.XSize, red_band.YSize
    # block_xsize, block_ysize = red_band.GetBlockSize()

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

            red = np.ma.masked_where(nir + red + 0.5 == 0, red)
            savi = ((1.5 * (nir - red)) / (nir + red + 0.5))*1000
            savi = savi.filled(-9999)

            out_data[y:y+rows, x:cols] = savi

        out_savi = create_raster(gdal.Open(VHR_repr_filenames['RED']), os.path.join(os.getcwd(), VI_filenames['SAVI']), out_data, gdal.GDT_Int16, -9999)
        print("Finished SAVI calculation - TIFF file created")


        ############################################
        # Reproject cloudmask and convert to tiff #
        ############################################
        # Cloudmask is currently an shp file. We'll convert it to a tiff and reproject it
        # to the EPSG:3857 projection that is used in the Scripts.

        # Let's initially convert the cloudmask to a tif file
        # Get the file paths (this needs to be automated)
        # cloudmask_file_shp_reprojected = os.path.join(os.getcwd(), '../Cloudmask/cloudmask.shp')
        cloudmask_file_shp_reprojected = '/' + '/'.join(cloudmask.split('/')[:-1]) + '/cloudmask.shp'
        cloudmask_file_out_tif = '/' + '/'.join(cloudmask.split('/')[:-1]) + '/cloudmask.tif'

        # Reproject the cloud mask to EPSG:3857
        os.system('ogr2ogr -t_srs EPSG:3857 {0} {1}'.format(cloudmask_file_shp_reprojected, cloudmask))

        # Convert the cloudmask shapefile to raster
        convert2Raster(cloudmask_file_shp_reprojected, cloudmask_file_out_tif)

        print("Cloudmask converted to tiff")


    ############################################
    # Clip rasters generated against cloudmask #
    ############################################
    # Since there are slight differences in the size of the single-band tiffs and the cloudmask
    # provided we will clip them against on another. In fact, this was found to be the case with
    # the Dutch Orthophotos, it will not necessarily be a problem with the ESA CSC DAP data that
    # we got. Also it is not only a problem of the single-band reprojected tiffs that were generated
    # as a result of the processing to now, but it was also a problem with the original tif.

    # Iterate through the tiffs we have created (both the bands and the vegetation indices)
    # and clip them against the cloudmask.

    # Get extent of one of the reprojected band tiffs
    # Get the extent of the cloudmask tif
    cl_tif = gdal.Open(os.path.join(os.getcwd(), cloudmask_file_out_tif))
    # We will need the GeoTransform and its inverse
    cl_gt = cl_tif.GetGeoTransform()
    cl_inv_gt = gdal.InvGeoTransform(cl_gt)
    # Get its extent
    cl_xmin = cl_gt[0]
    cl_xmax = cl_xmin + (cl_gt[1] * cl_tif.RasterXSize)
    cl_ymin = cl_gt[3] + (cl_gt[5] * cl_tif.RasterYSize)
    cl_ymax = cl_gt[3]

    for f in list(VHR_repr_filenames.values()) + list(VI_filenames.values()):
        ###
        # 1st -- Clip image based on the extent of the cloudmask
        ###
        # clip_raster(os.path.join(os.getcwd(),f), os.path.join(os.getcwd(),'../Cloudmask/cloudmask.tif'))
        ### Trying warp
        os.system('gdalwarp -te {0} {1} {2} {3} -tr 10 10 -r bilinear {4} {5}'
            .format(
                cl_xmin,
                cl_ymin,
                cl_xmax,
                cl_ymax,
                os.path.join(os.getcwd(),f),
                f.split('.')[0] + '_intersection.tif'
            )
        )
        print('Intersected file created: {0}'.format(f.split('.')[0] + '_intersection.tif'))

    del red_ds
    del nir_ds
    del out_savi
    del out_ndvi

# Does this even work now?
