try:
    from urllib import urlretrieve
except ImportError:
    from urllib.request import urlretrieve
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-a', '--all', action='store_true',
    help='Get all packages, not just the ones need for releasing SnapPy')

args = parser.parse_args()

url_base = 'https://download.microsoft.com/download/'

# Download links from
#
# https://support.microsoft.com/en-us/help/2977003/the-latest-supported-visual-c-downloads
#
# and
#
# https://stackoverflow.com/questions/12206314/detect-if-visual-c-redistributable-for-visual-studio-2012-is-installed

all_url_details = [
    ('2008', '5/D/8/5D8C65CB-C849-4025-8E95-C3966CAFD8AE/', 'vcredist_x86.exe'),
    ('2008', '5/D/8/5D8C65CB-C849-4025-8E95-C3966CAFD8AE/', 'vcredist_x64.exe'),
    ('2010', '1/6/5/165255E7-1014-4D0A-B094-B6A430A6BFFC/', 'vcredist_x86.exe'),
    ('2010', '1/6/5/165255E7-1014-4D0A-B094-B6A430A6BFFC/', 'vcredist_x64.exe'),
    ('2015', '6/A/A/6AA4EDFF-645B-48C5-81CC-ED5963AEAD48/', 'vc_redist.x86.exe'),
    ('2015', '6/A/A/6AA4EDFF-645B-48C5-81CC-ED5963AEAD48/', 'vc_redist.x64.exe'),
]

url_details = [
    ('2008', '5/D/8/5D8C65CB-C849-4025-8E95-C3966CAFD8AE/', 'vcredist_x86.exe'),
    ('2015', '6/A/A/6AA4EDFF-645B-48C5-81CC-ED5963AEAD48/', 'vc_redist.x86.exe'),
]

if args.all:
    url_details = all_url_details

for year, url, filename in url_details:
    arch = filename[-7:-4]
    target = 'vc_redist_%s_%s.exe' % (year, arch)
    print('Downloading %s...' % target)
    urlretrieve(url_base + url + filename, target)
print('Downloads complete!')
