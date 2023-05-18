import os
import webbrowser
from threading import Thread
from http.server import ThreadingHTTPServer, SimpleHTTPRequestHandler
from . import __path__

class SnapPyDocHandler(SimpleHTTPRequestHandler):
    def __init__(self, *args, **kwargs):
        doc_dir = os.path.join(__path__[0], 'doc')
        super().__init__(*args, directory=doc_dir, **kwargs)

    def log_message(self, *args, **kwargs):
        pass

class SnapPyDocServer(ThreadingHTTPServer):
    def __init__(self, *args, **kwargs):
        super().__init__(('127.0.0.1', 0), SnapPyDocHandler)
        self.URL = 'http://%s:%s'%(self.server_address)
        self.thread = Thread(target=self.serve_forever, daemon=True)
        self.thread.start()
    

