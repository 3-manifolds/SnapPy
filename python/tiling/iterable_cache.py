from itertools import count

__all__ = ['IterableCache']

class IterableCache:
    def __init__(self, iterable):
        self._iterable = iterable
        self._cache = []

    def __iter__(self):
        for i in count():
            if i == len(self._cache):
                try:
                    n = next(self._iterable)
                except StopIteration:
                    break
                self._cache.append(n)
            yield self._cache[i]


