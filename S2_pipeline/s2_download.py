import logging
import requests
import os
import time
import shutil
from pathlib import Path
from settings import SETTINGS, CREDENTIALS

tiles        = SETTINGS['tiles']
outdir       = SETTINGS['download_outdir']
start_date   = SETTINGS['start_date']
end_date     = SETTINGS['end_date']
download_url = SETTINGS['download_url']
token_url    = SETTINGS['token_url']

username = CREDENTIALS['username']
password = CREDENTIALS['password']

for tile_name in tiles:

    print("Tile: {}".format(tile_name))

    # This needs to be dynamically generated, and all the parameters passed should be configurable through the settings
    query_url = 'https://finder.creodias.eu/resto/api/collections/Sentinel2/search.json?maxRecords=10&rocessingLevel=LEVEL2A&productIdentifier=%25' + tile_name + '%25&startDate=' + start_date + 'T00%3A00%3A00Z&completionDate=' + end_date + 'T23%3A59%3A59Z&sortParam=startDate&sortOrder=descending&status=0%7C34%7C37&dataset=ESA-DATASET&cloudCover=%5B0%2C90%5D&'

    initial_outdir = outdir

    def _get_next_page(links):
        for link in links:
            if link['rel'] == 'next':
                return link['href']
        return False


    query_response = {}
    while query_url:
        response = requests.get(query_url)
        response.raise_for_status()
        data = response.json()
        for feature in data['features']:
            query_response[feature['id']] = feature
        query_url = _get_next_page(data['properties']['links'])


    def _get_token(username, password):
        token_data = {
            'client_id': 'CLOUDFERRO_PUBLIC',
            'username': username,
            'password': password,
            'grant_type': 'password'
        }
        response = requests.post(token_url, data=token_data).json()
        try:
            return response['access_token']
        except KeyError:
            raise RuntimeError('Unable to get token. Response was {0}'.format(response))


    ids = [result['id'] for result in query_response.values()]

    print("Found the following ids:")
    for curid in ids:
        print(curid)

    i = 0
    isUpdated = 0

    logging.info("Found %d new products for tile %s..." % (len(ids),tile_name))
    d = {}
    for id in ids:
        i += 1
        time.sleep(0.1)
        outdir = initial_outdir
        identifier = query_response[id]['properties']['productIdentifier'].split('/')[-1]
        if 'MSIL1C' in identifier:
            continue
        sensing_date = identifier.split('.')[0].split('_')[2]
        print(sensing_date)
        sensing_day = sensing_date[-2:]
        sensing_year = sensing_date[:4]
        sensing_month = sensing_date[4:6]
        size = query_response[id]['properties']['services']['download']['size'] / 1024.0 / 1024.0

        token = _get_token(username, password)
        url = '{}/{}?token={}'.format(download_url, id, token)
        # pbar =  tqdm(total = len(ids), unit='files')
        outdir = outdir + '/' + sensing_year
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        outdir = outdir + '/' + tile_name
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        outdir = outdir + '/' + sensing_month
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        print(outdir)
        outfile = Path(outdir) / '{}.zip'.format(identifier)
        if os.path.exists(outfile):
            continue
        outfile_temp = str(outfile) + '.incomplete'
        try:
            downloaded_bytes = 0
            print(url)
            req = requests.get(url, stream=True, timeout=100)
            chunk_size = 2 ** 20  # download in 1 MB chunks
            with open(outfile_temp, 'wb') as fout:
                for chunk in req.iter_content(chunk_size=chunk_size):
                    if chunk:  # filter out keep-alive new chunks
                        fout.write(chunk)
                        downloaded_bytes += len(chunk)
            shutil.move(outfile_temp, str(outfile))
            isUpdated = 1
            msg = "Download for product "+id+" has been completed!"
            logging.info("%s" % msg)
        finally:
            try:
                Path(outfile_temp).unlink()
            except OSError:
                pass
