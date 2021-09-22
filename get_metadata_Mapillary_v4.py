import mercantile, requests, json
from vt2geojson.tools import vt_bytes_to_geojson

# define an empty geojson as output
output= { "type": "FeatureCollection", "features": [] }

# vector tile endpoints -- change this in the API request to reference the correct endpoint
tile_coverage = 'mly1_public'

# tile layer depends which vector tile endpoints:
# 1. if map features or traffic signs, it will be "point" always
# 2. if looking for coverage, it will be "image" for points, "sequence" for lines, or "overview" for far zoom
tile_layer = "image"

# Mapillary access token -- user should provide their own
access_token = 'MLY|XXXX'

# a bounding box in [east_lng,_south_lat,west_lng,north_lat] format
west, south, east, north = [4.7, 51.5, 5.7, 52.5]

# get the list of tiles with x and y coordinates which intersect our bounding box
# MUST be at zoom level 14 where the data is available, other zooms currently not supported
tiles = list(mercantile.tiles(west, south, east, north, 14))

start_date='2017-04-01'
end_date='2017-11-01'

# loop through list of tiles to get tile z/x/y to plug in to Mapillary endpoints and make request
for tile in tiles:
    tile_url = 'https://tiles.mapillary.com/maps/vtp/{}/2/{}/{}/{}?access_token={}&start_time={}&end_time={}'.format(
        tile_coverage,tile.z,tile.x,tile.y,access_token, start_date, end_date)

    response = requests.get(tile_url)
    data = vt_bytes_to_geojson(response.content, tile.x, tile.y, tile.z,layer=tile_layer)

    # push to output geojson object if yes
    for feature in data['features']:
        output['features'].append(feature)

# save a local geojson with the filtered data
with open('images.geojson', 'w') as f:
    json.dump(output, f)