
from osgeo import gdal
from osgeo import ogr
from osgeo import osr
from settings import SETTINGS
import os
import pdb

# Taken from S2 scripts
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

        del in_ras
        del ras_vrt


# Convert shapefile to raster
def convert2Raster(shapefile, output_raster):

    # Retreve the VHR data source
    vhr_src = SETTINGS['VHR_SOURCE']

    # Open the shapefile and extract the single layer (shapefiles have 1 layer AFAIK)
    shp = ogr.Open(shapefile)
    shp_layer = shp.GetLayer()

    # We want to match the resolution to the Sentinel optical bands resolution, thus 10m
    # If more scenarios come up in the future, we can convert this to a configuration item
    # in the SETTINGS dictionary
    pixel_size = 10

    # Get the layer extent and calculate the final image resolution
    x_min, x_max, y_min, y_max = shp_layer.GetExtent()
    x_res = int((x_max - x_min) / pixel_size)
    y_res = int((y_max - y_min) / pixel_size)

    # Prepare for the new GeoTiff creation
    driver = gdal.GetDriverByName('GTiff')
    new_raster = driver.Create(output_raster, x_res, y_res, 1, gdal.GDT_Byte)
    new_raster.SetGeoTransform((x_min, pixel_size, 0, y_min, 0, pixel_size))
    band = new_raster.GetRasterBand(1)
    no_data_value = -9999 # How does the -9999 occur? Is it sentinel-specific?
    band.SetNoDataValue(no_data_value)
    band.FlushCache()

    # Eventually rasterize the layer we got from the shapefile
    gdal.RasterizeLayer(
            new_raster,
            [1],
            shp_layer,
            options = ['ATTRIBUTE={}'.format(SETTINGS['cloudmask_attr_name'][vhr_src])] # from settings
    )

    # The projection of the new tiff will be EPSG:3857. In case we may need more
    # in the future, we can make this configurable as well.
    new_rasterSRS = osr.SpatialReference()
    new_rasterSRS.ImportFromEPSG(3857)
    new_raster.SetProjection(new_rasterSRS.ExportToWkt())

    # return
    return gdal.Open(output_raster)

