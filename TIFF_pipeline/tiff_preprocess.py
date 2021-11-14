# Attempt to preprocess tiff files in general (eg. coming from VHR acquisitions like PlanetScope)
import rasterio as rio

tiff_file = '/Users/gchoumos/Downloads/callisto_temp/Dutch_VHR/as_on_my_notes_aug_4/20170923_094819_Tri_12bit/MS/TRIPLESAT_3_MS_L1_20170923094819_001237VI_008_0220170831001001_029.tif'
raster = rio.open(tiff_file)

print('Raster file read')
print('Shape: {}'.format(raster.shape))
print('Number of bands: {}'.format(raster.count))

# Read the bands
band_array = raster.read()

# It seems that it is common for all Planetscope and other satellite acquisitions to store the
# bands in the raster using the following order
# 1. Blue
# 2. Green
# 3. Red
# 4. NIR

