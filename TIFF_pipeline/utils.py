
from osgeo import gdal
from osgeo import ogr
from osgeo import osr
from settings import SETTINGS
import os
import pdb
import json
import matplotlib.pyplot as plt
from descartes import PolygonPatch

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


def get_extent(in_ras_file):
    ''' Returns min_x, max_y, max_x, min_y'''

    in_ras = gdal.Open(in_ras_file)
    gt = in_ras.GetGeoTransform()
    return (gt[0], gt[3], gt[0] + gt[1] * in_ras.RasterXSize, gt[3] + gt[5] * in_ras.RasterYSize)


# Utility function to clip a raster based on the extent of another
def clip_raster(input_tif, clip_tif):
    ''' Clip a raster (one band) - returns nothing'''

    # Get the extent of the input tif
    in_tif = gdal.Open(input_tif)
    # We will need the GeoTransform and its inverse
    in_gt = in_tif.GetGeoTransform()
    in_inv_gt = gdal.InvGeoTransform(in_gt)
    # Get its extent
    in_xmin = in_gt[0]
    in_xmax = in_xmin + (in_gt[1] * in_tif.RasterXSize)
    in_ymin = in_gt[3] + (in_gt[5] * in_tif.RasterYSize)
    in_ymax = in_gt[3]

    # Same for the clip tiff
    cl_tif = gdal.Open(clip_tif)
    cl_gt = cl_tif.GetGeoTransform()
    cl_inv_gt = gdal.InvGeoTransform(cl_gt)
    cl_xmin = cl_gt[0]
    cl_xmax = cl_xmin + (cl_gt[1] * cl_tif.RasterXSize)
    cl_ymin = cl_gt[3] + (cl_gt[5] * cl_tif.RasterYSize)
    cl_ymax = cl_gt[3]

    srs = osr.SpatialReference()
    srs.ImportFromEPSG(3857)

    off_ulx, off_uly = map(int, gdal.ApplyGeoTransform(in_inv_gt, cl_xmin, cl_ymax))
    off_lrx, off_lry = map(int, gdal.ApplyGeoTransform(in_inv_gt, cl_xmax, cl_ymin))

    rows, columns = abs(off_lry - off_uly), abs(off_lrx - off_ulx)
    in_band = in_tif.GetRasterBand(1)

    out_prefix = input_tif.split('/')[-1].split('.')[0]
    out_suffix = input_tif.split('.')[-1]
    out_tif = '{0}_clipped.{1}'.format(out_prefix,out_suffix)
    driver = gdal.GetDriverByName('GTiff')
    out_ds = driver.Create(out_tif, columns, rows, 1, in_band.DataType)
    pdb.set_trace()
    out_ds.SetProjection(in_tif.GetProjection())
    ulx, uly = gdal.ApplyGeoTransform(in_gt, off_ulx, off_uly)
    out_gt = list(in_gt)
    out_gt[0], out_gt[3] = ulx, uly
    out_ds.SetGeoTransform(out_gt)

    out_ds.GetRasterBand(1).WriteArray(in_band.ReadAsArray(off_ulx, off_uly, columns, rows))

    del in_tif
    del in_band
    del out_ds

    print('Clipping of {0} based on the extent ({1},{2},{3},{4}) has been successfully performed'
            .format(input_tif.split('/')[-1],in_xmin,in_xmax,in_ymin,in_ymax))

    return


def plot_geometry(geometry):
    geom_string = geometry.ExportToJson()
    geom_dict = json.loads(geom_string)
    BLUE = '#6699cc'
    fig = plt.figure()
    ax = fig.gca()
    ax.add_patch(PolygonPatch(geom_dict, fc=BLUE, ec=BLUE, alpha=0.5, zorder=2 ))
    ax.axis('scaled')
    plt.show()
