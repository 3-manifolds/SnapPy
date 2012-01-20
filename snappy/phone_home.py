try:
    from httplib import HTTPConnection
except ImportError: # Python 3
    from http.client import HTTPConnection
from threading import Thread    
from snappy.version import version as old_version

class Phoner(Thread):
    def __init__(self):
        Thread.__init__(self)
        self.new_version = None
    def run(self):
        try:
            connection = HTTPConnection('www.math.uic.edu', timeout=1)
            connection.request('GET','/t3m/SnapPy-nest/current.txt')
            response = connection.getresponse()
            self.new_version = response.read()
            connection.close()
        except:
            pass
        if isinstance(self.new_version, bytes):
            self.new_version = self.new_version.decode()

def needs_updating():
    ET = Phoner()
    ET.start()
    ET.join(0.5)
    if ET.new_version and ET.new_version != old_version:
        return (ET.new_version, old_version)
    return None

if __name__ == '__main__':
    print(needs_updating())
