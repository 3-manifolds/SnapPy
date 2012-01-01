try:
    from httplib import HTTPConnection
except ImportError:
    from http.client import HTTPConnection
    
from snappy.version import version as old_version

def get_current():
    try:
        connection = HTTPConnection('www.math.uic.edu', timeout=1)
        connection.request('GET','/t3m/SnapPy-nest/current.txt')
        response = connection.getresponse()
        result = response.read()
        connection.close()
    except:
        return None
    return result.strip()

def needs_updating():
    new_version = get_current()
    if new_version and new_version != old_version:
        return (new_version, old_version)
    return None

if __name__ == '__main__':
    print(needs_updating())
