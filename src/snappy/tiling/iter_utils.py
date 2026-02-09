from itertools import count

class IteratorCache:
    def __init__(self, iterator):
        """
        Takes an iterator and produces an iterable such that an
        iterator created from the iterable produces the same iterated
        values.

        The iterates values, however, are lazily cached so that
        iterating over another iterator again is faster.

        (Silly) Example:

            >>> squares = IterableCache(iter(i ** 2 for i in range(1000)))
            >>> iterator0 = iter(squares)
            >>> next(iterator0)
            0
            >>> next(iterator0)
            1
            >>> next(iterator0)
            4
            >>> iterator1 = iter(squares)
            >>> next(iterator1) # Cached from previous run
            0
            >>> next(iterator1)
            1
            >>> next(iterator1)
            4
            >>> next(iterator1) # Not cached from previous run
            9

        """

        self._iterator = iterator
        self._cache = []

    def __iter__(self):
        for i in count():
            if i == len(self._cache):
                try:
                    n = next(self._iterator)
                except StopIteration:
                    break
                self._cache.append(n)
            yield self._cache[i]

def merge_iterators(iterators):
    """
    Merges iterators into one iterable picking the smallest element
    (using < operator) from the list of next elements for each iterator.
    """

    elements  = [ [ next(iter), i ]
                  for i, iter in enumerate(iterators) ]

    while True:
        element, i = min(elements)
        yield element

        elements[i][0] = next(iterators[i])

def merge_iterables(iterables):
    """
    Merges iterables into one iterable picking the smallest element
    (using < operator) from the list of next elements for each iterator.

    Example:

        >>> squares_and_cubes = iter(merge_iterables([
        ...          (i ** 2 for i in range(10000)),
        ...          (i ** 3 for i in range(10000)) ]))
        >>> list(islice(squares_and_cubes, 12))
        [0, 0, 1, 1, 4, 8, 9, 16, 25, 27, 36, 49]

    """

    return merge_iterators([ iter(i) for i in iterables ])
