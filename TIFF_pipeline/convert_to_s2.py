
from osgeo import gdal, osr
from utils import *
import numpy as np
import os
import pdb


#####################################
# Utility function -- Create Raster #
#####################################
'''Create a one-band GeoTIFF '''
def create_raster(template_raster, template_path, np_data, data_type, nodata=None):
    driver = gdal.GetDriverByName('GTiff')
    out_ds = driver.Create(template_path, template_raster.RasterXSize, template_raster.RasterYSize, 1, data_type)
    out_ds.SetProjection(template_raster.GetProjection())
    out_ds.SetGeoTransform(template_raster.GetGeoTransform())
    out_band = out_ds.GetRasterBand(1)
    if nodata is not None:
        out_band.SetNoDataValue(nodata)
    out_band.WriteArray(np_data)
    out_band.FlushCache()
    out_band.ComputeStatistics(False)

    return


#############################################
# Pending - Split Multiband to Single Bands #
#############################################
# Code that converts 1 multiband tiff to N single-band tiffs
input_file = 'TRIPLESAT_3_MS_L1_20170923094819_001237VI_008_0220170831001001_029.tif'
band_names = ['BLUE', 'GREEN', 'RED', 'NIR']

data = gdal.Open(input_file)
print('Number of bands: {}'.format(data.RasterCount))
# Extract the bands to single image tiffs
for band in range(data.RasterCount):
    print('Extracting band {} ...'.format(band+1))
    file_prefix = input_file.split('.')[0]
    file_suffix = input_file.split('.')[-1]
    os.system(
        'gdal_translate {0} -b {1} {2}_{3}.{4}'
            .format(input_file, band+1, file_prefix, band_names[band], file_suffix)
    )

# Checkpoint
pdb.set_trace()

###########################################
# Single-band filenames that were created #
###########################################
VHR_native_band_filenames = {
    # before reprojecttion
    'BLUE':  'TRIPLESAT_3_MS_L1_20170923094819_001237VI_008_0220170831001001_029_BLUE.tif',
    'GREEN': 'TRIPLESAT_3_MS_L1_20170923094819_001237VI_008_0220170831001001_029_GREEN.tif',
    'RED':   'TRIPLESAT_3_MS_L1_20170923094819_001237VI_008_0220170831001001_029_RED.tif',
    'NIR':   'TRIPLESAT_3_MS_L1_20170923094819_001237VI_008_0220170831001001_029_NIR.tif',
}

##########################################################
# Single-band reprojected filenames that will be created #
##########################################################
VHR_repr_filenames = {
    # after reprojecttion
    'BLUE':  'TRIPLESAT_3_MS_L1_20170923094819_001237VI_008_0220170831001001_029_BLUE_3857.tif',
    'GREEN': 'TRIPLESAT_3_MS_L1_20170923094819_001237VI_008_0220170831001001_029_GREEN_3857.tif',
    'RED':   'TRIPLESAT_3_MS_L1_20170923094819_001237VI_008_0220170831001001_029_RED_3857.tif',
    'NIR':   'TRIPLESAT_3_MS_L1_20170923094819_001237VI_008_0220170831001001_029_NIR_3857.tif',
}

VI_filenames = {
    # Vegetation Indices
    'NDVI': 'TRIPLESAT_3_MS_L1_20170923094819_001237VI_008_0220170831001001_029_NDVI_3857.tif',
    'SAVI': 'TRIPLESAT_3_MS_L1_20170923094819_001237VI_008_0220170831001001_029_SAVI_3857.tif'
}

###############################
# Reproject single band tiffs #
###############################
s2_proj = 'EPSG:3857'

srs = osr.SpatialReference()
srs.ImportFromEPSG(3857)

for f_input in VHR_native_band_filenames.values():
    print('Reprojection of {0} ...'.format(f_input))

    # in_ras = gdal.Open(f_input)
    # projected_vrt = gdal.AutoCreateWarpedVRT(in_ras, None, srs.ExportToWkt(), gdal.GRA_Bilinear)
    # out_ds = gdal.GetDriverByName('GTiff').CreateCopy(
    #     f_input.split('.')[0]+ '_' + s2_proj.split(':')[0] + '_' + s2_proj.split(':')[1] + '.tif',
    #     projected_vrt
    # )

    # del in_ras
    # del projected_vrt
    # del out_ds

    file_prefix = f_input.split('.')[0]
    file_suffix = f_input.split('.')[-1]
    proj_name_part = s2_proj.split(':')[1]

    # Running through a GDAL system call
    # os.system('gdalwarp -t_srs {0} -r bilinear -tr 10 10 {1} {2}_{3}.{4}'.format(s2_proj, f_input, file_prefix, proj_name_part, file_suffix))

    # Running through the GDAL python library
    gdal.Warp(
        # Output file will be the same as the previous but with the projection at the end
        '{0}_{1}.{2}'.format(file_prefix, proj_name_part, file_suffix),
        f_input,
        xRes = 10,
        yRes = 10,
        dstSRS=s2_proj,
        resampleAlg='bilinear'
    )

#########################
# NDVI Creation attempt #
#########################
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

#########################
# SAVI Creation attempt #
#########################
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
# Clip rasters generated against cloudmask #
############################################
# Since there are slight differences in the size of the single-band tiffs and the cloudmask
# provided we will clip them against on another. In fact, this was found to be the case with
# the Dutch Orthophotos, it will not necessarily be a problem with the ESA CSC DAP data that
# we got. Also it is not only a problem of the single-band reprojected tiffs that were generated
# as a result of the processing to now, but it was also a problem with the original tif.

# Let's initially convert the cloudmask to a tif file
# Get the file paths (this needs to be automated)
cloudmask_file_in = os.path.join(os.getcwd(),'../Cloudmask/TRIPLESAT_3_L1_20170923094819_001237VI_008_0220170831001001_029.shp')
cloudmask_file_shp_reprojected = os.path.join(os.getcwd(), '../Cloudmask/cloudmask.shp')
cloudmask_file_out_tif = os.path.join(os.getcwd(), '../Cloudmask/cloudmask.tif')

# Reproject the cloud mask to EPSG:3857
os.system('ogr2ogr -t_srs EPSG:3857 {0} {1}'.format(cloudmask_file_shp_reprojected, cloudmask_file_in))

# Convert the cloudmask shapefile to raster
convert2Raster(cloudmask_file_shp_reprojected, cloudmask_file_out_tif)

print("Converted to tiff")

####


# Iterate through the tiffs we have created (both the bands and the vegetation indices)
# for f in VHR_repr_filenames + VI_filenames:




del red_ds
del nir_ds
del out_savi
del out_ndvi
