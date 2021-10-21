
# Settings dict with the details required for the scipts to work
#
# These scripts are:
# --> s2_download.py
# --> s2_preprocess.py
# --> s2_feature_space.py
#
# The details here should be changed accordingly if anything gets updated.
# It is recommended, especially we eventually plan to publish the code, that we set the
# username and password either through the credentials.py file (which is gitignored) or
# we can prompt the users to enter them on execution
#
SETTINGS = {
    'tiles': ['31UFT'],
    'outdir': '/home/callisto-noa/callisto/data/s2',
    'start_date': '2017-03-01',
    'end_date': '2017-11-30',
    'download_url': 'https://zipper.creodias.eu/download',
    'token_url': 'https://auth.creodias.eu/auth/realms/DIAS/protocol/openid-connect/token',
    'months': {
        '01': 'January_all',
        '02': 'February_all',
        '03': 'March_all',
        '04': 'April_all',
        '05': 'May_all',
        '06': 'June_all',
        '07': 'July_all',
        '08': 'August_all',
        '09': 'September_all',
        '10': 'October_all',
        '11': 'November_all',
        '12': 'December_all',
    },
    'bands': [
        'B02_10m.tif', 'B03_10m.tif', 'B04_10m.tif', 'B05_10m.tif',
        'B06_10m.tif', 'B07_10m.tif', 'B08_10m.tif', 'B8A_10m.tif',
        'B11_10m.tif', 'B12_10m.tif', 'SCL_10m.tif', 'TCI_10m.tif',
    ],

}

# Either update the credentials here, or create a credentials.py file defining the dict there
CREDENTIALS = {
    'username': 'placeholder_username',
    'password': 'placeholder_password',
    'sudo_pass_bsc_vm': 'placeholder_sudo_pass_bsc_vm',
}

# Try to import the redacted credentials, if they exist
# Gracely handle the case it doesn't
try:
    from credentials import CREDENTIALS
except ImportError:
    print("No valid redacted credentials found. Will use the ones in the settings.py file")
