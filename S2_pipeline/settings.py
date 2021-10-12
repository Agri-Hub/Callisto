
# Settings dict with the details required for the scipts to work
#
# The details here should be changed accordingly if anything gets updated.
# It is recommended, especially we eventually plan to publish the code, that we set the
# username and password either through the credentials.py file (which is gitignored) or
# we can prompt the users to enter them on execution
#
SETTINGS = {
    'tiles': ['31UFT'],
    'outdir': '/home-ext/callisto/akoukos/data/s2',
    'start_date': '2017-03-01',
    'end_date': '2017-11-30',
    'download_url': 'https://zipper.creodias.eu/download',
    'token_url': 'https://auth.creodias.eu/auth/realms/DIAS/protocol/openid-connect/token',
}

# Either update the credentials here, or create a credentials.py file defining the dict there
CREDENTIALS = {
    'username': 'placeholder_username',
    'password': 'placeholder_password'
}

# Try to import the redacted credentials, if they exist
# Gracely handle the case it doesn't
try:
    from credentials import CREDENTIALS
except ImportError:
    print("No valid redacted credentials found. Will use the ones in the settings.py file")
