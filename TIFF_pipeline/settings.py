# The details here should be changed accordingly if anything gets updated.
# It is recommended, especially we eventually plan to publish the code, that we set the
# username and password either through the credentials.py file (which is gitignored) or
# we can prompt the users to enter them on execution
#
SETTINGS = {
    'working_directory': '/Users/george/Downloads/callisto_temp/Dutch_VHR/dandrimont_timeframe/2017',
    'shapefile': '/Users/george/Downloads/callisto_temp/shapefiles/lpis_dandrimond_area/lpis_dandrimond_3857.shp',
    'feature_space_filename': 'VHR_fs_dandrimont.csv',
    'bands': [
        'B02_10m.tif', 'B03_10m.tif', 'B04_10m.tif', 'B05_10m.tif',
        'B06_10m.tif', 'B07_10m.tif', 'B08_10m.tif', 'B8A_10m.tif',
        'B11_10m.tif', 'B12_10m.tif', 'SCL_10m.tif', 'TCI_10m.tif',
    ],
    # The VHR Source. This is useful for various reasons, like the expected file
    # structure, the cloumask attribute etc.
    # Possible types up to now:
    # - DUTCH_ORTHO
    #     The orthophotos that we got for the Netherlands from the site: satellietdataportaal.nl/
    #     This source was pointed by Marcel from RVO, a Dutch PA.
    'VHR_SOURCE': 'DUTCH_ORTHO',
    # The attribute name that will be used in the cloudmask
    'cloudmask_attr_name': {
        'DUTCH_ORTHO': 'DN',
    },
}
