try:
    from httplib import HTTPConnection
except ImportError: # Python 3
    from http.client import HTTPConnection
from threading import Thread    
from snappy.version import version as this_version

class Phoner(Thread):
    def __init__(self):
        Thread.__init__(self)
        self.answer = None
    def run(self):
        newest_version = None
        try:
            connection = HTTPConnection('www.math.uic.edu', timeout=1)
            connection.request('GET','/t3m/SnapPy-nest/current.txt')
            response = connection.getresponse()
            newest_version = response.read().strip()
            connection.close()
        except:
            pass
        if isinstance(newest_version, bytes):
            newest_version = newest_version.decode()
        if newest_version and newest_version != this_version:
            self.answer = (newest_version, this_version)

def update_needed():
    ET = Phoner()
    ET.start()
    ET.join(1.0) # Wait up to 1.0 seconds for an answer
    if ET.is_alive():
        return ''
    elif ET.answer:
        return ("**Please upgrade to %s from %s via "
                "http://snappy.computop.org**\n" % ET.answer)
    else:
         return ''
        
if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        this_version = '0.0.0'
    else:
        print('To provoke an update response add an argument.')
    import time
    start = time.time()
    print(update_needed())
    print('Elapsed time: %ss'%(time.time() - start))
