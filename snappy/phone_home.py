try:
    from httplib import HTTPConnection
except ImportError: # Python 3
    from http.client import HTTPConnection
from threading import Thread    
from snappy.version import version as old_version

class Phoner(Thread):
    def __init__(self):
        Thread.__init__(self)
        self.answer = None
    def run(self):
        try:
            connection = HTTPConnection('www.math.uic.edu', timeout=1)
            connection.request('GET','/t3m/SnapPy-nest/current.txt')
            response = connection.getresponse()
            self.new_version = response.read().strip()
            connection.close()
        except:
            pass
        if isinstance(self.new_version, bytes):
            self.new_version = self.new_version.decode()
        if self.new_version and self.new_version != old_version:
            self.answer = (self.new_version, old_version)

if __name__ == '__main__':
    ET = Phoner()
    ET.start()
    ET.join()
    print(ET.answer)
