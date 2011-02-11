from httplib import HTTPConnection

def get_current():
    try:
        connection = HTTPConnection('www.math.uic.edu')
        connection.request('GET','/t3m/SnapPy-nest/current.txt')
        response = connection.getresponse()
        result = response.read()
        connection.close()
    except:
        return None
    return result.strip()

if __name__ == '__main__':
current = get_current()
if current:
    print current
else:
    print 'Connection failed.'
