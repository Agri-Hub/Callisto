import geopandas
import pandas as pd
from shapely.geometry import Point, Polygon
import matplotlib.pyplot as plt
from deep_translator import (GoogleTranslator, PonsTranslator, LingueeTranslator, MyMemoryTranslator, YandexTranslator, 
                            DeepL, QCRI, single_detection, batch_detection)
from datetime import datetime, timezone, timedelta
import json, time, math

# functions that calculates a new point based on distance and orientation from another
def new_lat(row, d = 0.02, slide = 45):
    lat = math.radians(row.lat)
    lon = math.radians(row.lon)
    theta = math.radians((row.compass_angle + slide)% 360)
    sinf = math.sin(lat)
    cosd = math.cos(d/6371)
    cosf = math.cos(lat) 
    sind = math.sin(d/6371)
    costheta = math.cos(theta)
    cosf = math.cos(lat) 
    return math.degrees(math.asin(sinf*cosd + cosf*sind*costheta))
def new_lon(row, d = 0.02, slide = 45):
    lat = math.radians(row.lat)
    lon = math.radians(row.lon)
    theta = math.radians((row.compass_angle + slide)% 360)
    new_lat = math.radians(row.new_lat) 
    sintheta = math.sin(theta)
    cosd = math.cos(d/6371)
    sind = math.sin(d/6371)
    cosf1 = math.cos(lat) 
    sinf1 = math.sin(lat) 
    sinf2 = math.sin(new_lat) 
    return row.lon + math.degrees(math.atan2(sintheta*sind*cosf1, cosd - sinf1*sinf2))

# Read the shapefile with the area of interest
gdp = geopandas.read_file(r'data\shapefiles\utrecht_4326.shp')
gdp.columns = ['CAT_GEWASCATEGORIE', 'GWS_GEWASCODE','GEOMETRIE_Length','GEOMETRIE_Area','GWS_GEWAS','id','geometry']

# Translate the description to english
english = {x:GoogleTranslator(source='nl', target='en').translate(x) for x in gdp.GWS_GEWAS.unique()}
gdp['GWS_GEWAS_eng'] = gdp['GWS_GEWAS'].map(english)
gdp.GWS_GEWASCODE = gdp.GWS_GEWASCODE.astype(int)

# Merge some basic categories
gdp["GWS_GEWAS_MERGED"] = 'Other'
gdp.loc[gdp.GWS_GEWASCODE.isin([265,266,331,336,383]),'GWS_GEWAS_MERGED'] = 'Grassland'
gdp.loc[gdp.GWS_GEWASCODE.isin([259,316,317]),'GWS_GEWAS_MERGED'] = 'Maize'
gdp.loc[gdp.GWS_GEWASCODE.isin([2014, 2015, 2016, 2017]),'GWS_GEWAS_MERGED'] = 'Potatos'
gdp.loc[gdp.GWS_GEWASCODE.isin([233]),'GWS_GEWAS_MERGED'] = 'Winter Wheat'
gdp.loc[gdp.GWS_GEWASCODE.isin([256]),'GWS_GEWAS_MERGED'] = 'Sugar Beet'
gdp.loc[gdp.GWS_GEWASCODE.isin([262, 1931]),'GWS_GEWAS_MERGED'] = 'Onions'
gdp.loc[gdp.GWS_GEWASCODE.isin([236]),'GWS_GEWAS_MERGED'] = 'Sumer Barley'
gdp.loc[gdp.GWS_GEWASCODE.isin([1004, 1002]),'GWS_GEWAS_MERGED'] = 'Flower'

# Load the metadata we got from Mapillary
f = open(r'data\images.geojson')
seq = json.load(f)
seq2 = []
for s in seq['features']:
    seq2.append({**s['properties'], **s['geometry']})
df = pd.DataFrame(seq2)

del seq, seq2

# Create a DataFrame with the images' metadata
df.captured_at = df.captured_at.map(lambda x: datetime.utcfromtimestamp(x//1000))
df = df[(df.captured_at > '2017-04-01') & (df.captured_at < '2017-11-01')]
df["lon"] = df.coordinates.map(lambda x: x[0])
df["lat"] = df.coordinates.map(lambda x: x[1])


# Calculate minimum and maximum coordinates of the parcels and then filter the Mapillary-images DataFrame 
minx = gdp.bounds.min().minx
miny = gdp.bounds.min().miny
maxx = gdp.bounds.max().maxx
maxy = gdp.bounds.max().maxy
df = df.loc[(df.lon >= minx) & (df.lat >= miny) & (df.lon <= maxx) & (df.lat <= maxy)]

# Calulcate new coords for right side of images
df["new_lat"] = df.apply(lambda x: new_lat(x, d = 0.01, slide = 45), axis = 1)
df["new_lon"] = df.apply(lambda x: new_lon(x, d = 0.01, slide = 45), axis = 1)

# Convert  to GeoDataFrame
new_points_gdf = geopandas.GeoDataFrame(df.copy(), geometry=geopandas.points_from_xy(df.new_lon, df.new_lat))
new_points_gdf["parcel_id"] = -1
new_points_gdf["cropType"] = None
start_time = time.time()

# For each parcel in gdp, find all images that illustrate it on the right side.
for i, row in gdp.iterrows():
    subset_points = new_points_gdf.loc[(new_points_gdf.new_lon > row.geometry.bounds[0]) & (new_points_gdf.new_lat > row.geometry.bounds[1]) &
                   (new_points_gdf.new_lon < row.geometry.bounds[2]) & (new_points_gdf.new_lat < row.geometry.bounds[3])].copy()
    subset_points = subset_points[subset_points.geometry.map(lambda x: x.within(row.geometry))].copy().index
    if len(subset_points) > 0:
        new_points_gdf.loc[subset_points, 'parcel_id'] = row.id
        new_points_gdf.loc[subset_points, 'cropType'] = row.GWS_GEWAS
print("Time elapsed: {}".format(time.time() - start_time))

new_points_gdf.to_file(r'data\str_level_right.json', driver="GeoJSON")

# Calulcate new coords for left side of images
df["new_lat"] = df.apply(lambda x: new_lat(x, d = 0.01, slide = -45), axis = 1)
df["new_lon"] = df.apply(lambda x: new_lon(x, d = 0.01, slide = -45), axis = 1)

# Same as before, convert to GeoDataFrame
new_points_gdf = geopandas.GeoDataFrame(df.copy(), geometry=geopandas.points_from_xy(df.new_lon, df.new_lat))
new_points_gdf["parcel_id"] = -1
new_points_gdf["cropType"] = None
start_time = time.time()

# For each parcel in gdp, find all images that illustrate it on the left side.
for i, row in gdp.iterrows():
    subset_points = new_points_gdf.loc[(new_points_gdf.new_lon > row.geometry.bounds[0]) & (new_points_gdf.new_lat > row.geometry.bounds[1]) &
                   (new_points_gdf.new_lon < row.geometry.bounds[2]) & (new_points_gdf.new_lat < row.geometry.bounds[3])].copy()
    subset_points = subset_points[subset_points.geometry.map(lambda x: x.within(row.geometry))].copy().index
    if len(subset_points) > 0:
        new_points_gdf.loc[subset_points, 'parcel_id'] = row.id
        new_points_gdf.loc[subset_points, 'cropType'] = row.GWS_GEWAS
print("Time elapsed: {}".format(time.time() - start_time))

new_points_gdf.to_file(r'data\str_level_left.json', driver="GeoJSON")

# Merge left and right images in one GeoDataFrame and store it.
new_points_gdf_right = geopandas.read_file('data/str_level_right.geojson')
new_points_gdf_left = geopandas.read_file('data/str_level_left.geojson')
new_points_gdf = new_points_gdf_right.iloc[:,:9]
new_points_gdf["parcel_left"] = new_points_gdf_left["parcel_id"]
new_points_gdf["parcel_right"] = new_points_gdf_right["parcel_id"]
new_points_gdf["cropType_left"] = new_points_gdf_left["cropType"]
new_points_gdf["cropType_right"] = new_points_gdf_right["cropType"]
new_points_gdf["geometry"] = new_points_gdf_right.apply(lambda x: Point(x.lon, x.lat), axis = 1)
new_points_gdf = geopandas.GeoDataFrame(new_points_gdf)
new_points_gdf.to_file(r'data\str_level.json', driver="GeoJSON")