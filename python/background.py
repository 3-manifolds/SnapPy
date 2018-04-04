from __future__ import print_function

"""
This module defines objects which compute information about a
manifold in a background process.  Examples of appropriate tasks for
such an object include its census identification, its length spectrum
up to a certain number of geodesics, or its covering spaces of a
specified degree.

A background helper object should only work on one manifold at a
time. It should start a worker process when given a manifold.  If it
already has a worker process running then it should kill that process
before starting work on the new manifold.  For example, a browser
window can use a BackgroundIdentifier to get a list of names to
put into its aka window.  If the user does a Dehn filling before
the BackgroundIdentifier has finished identifying the current
manifold, then the browser could send the filled manifold to the
BackgroundIdentifier.  It would stop trying to identify the old
manifold and work on the filled manifold instead.
"""

from multiprocessing import Process, Queue
from multiprocessing.queues import Empty
from pickle import dumps, loads

class Identifier(object):
    """
    Manages a background process which identifies manifolds.

    The get method of this object returns a dict with keys 'weak' and 'strong'.
    The values are lists of names of manifolds.  A strong identification means
    that there is an isometry which preserves meridians.
    """

    def __init__(self):
        self.process = self.queues = None
        
    @property
    def state(self):
        if self.process and self.process.is_alive():
            return 'working'
        elif self.queues is None:
            return 'ready'
        else:
            return 'finished'
        
    def identify(self, manifold):
        if self.state == 'working':
            self.process.terminate()
        self.queues = (Queue(100), Queue(100))
        self.process = Process(target=self._task, args=(dumps(manifold),))
        self.process.start()

    def _task(self, pickle):
        M = loads(pickle)
        weak_queue, strong_queue = self.queues
        for N in M.identify():
            name = N.name()
            weak_queue.put(name, False)
        for N in M.identify(True):
            name = N.name()
            strong_queue.put(name, False)
        
    def get(self):
        if not self.state == 'finished':
            raise RuntimeError('Identifier state is: "%s"; should be "finished".'%self.state)
        weak, strong = [], []
        weak_queue, strong_queue = self.queues
        while not weak_queue.empty():
            weak.append(weak_queue.get(False))
        while not strong_queue.empty():
            strong.append(strong_queue.get(False))
        self.queues = None
        print(weak, strong)
        return {'weak': weak, 'strong': strong}

if __name__ == '__main__':
    from time import sleep, time
    from snappy import Manifold
    M = Manifold('L14a31811')
    I = Identifier()
    print('starting identification process ...')
    I.identify(M)
    start = time()
    count = 5
    while True:
        sleep(0.5)
        elapsed = time() - start
        print('elapsed: %.3f'%(elapsed))
        if I.state == 'finished':
            print('done.')
            print(I.get())
            break
        if elapsed > 1.0 and count > 0:
            print('restarting')
            I.identify(M)
            start = time()
            count -= 1
