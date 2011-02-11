from httplib import HTTPConnection
import re

def get_current():
    try:
        connection = HTTPConnection('www.math.uic.edu', timeout=1)
        connection.request('GET','/t3m/SnapPy/index.html')
        response = connection.getresponse()
        result = response.read()
        connection.close()
        version = re.search("SnapPy v([0-9.]*) documentation", result).group(1)
    except:
        return None
    return version 

if __name__ == '__main__':
    current = get_current()
    if current:
        print current
    else:
        print 'Connection failed.'
