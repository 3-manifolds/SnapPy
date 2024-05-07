class DictBasedSet:
    """
    This class can be constructed from a dictionary with keys in A.

    The resulting object is a set of A.

    >>> fruit_set = DictBasedSet(dict())
    >>> fruit_set.add('Apple')
    True
    >>> fruit_set.add('Banana')
    True
    >>> fruit_set.add('Apple') # Already added earlier
    False

    """

    def __init__(self, dictionary):
        self._dictionary = dictionary

    def add(self, key) -> bool:
        """
        Adds key to set of A.

        Returns true if key was not yet in the set.
        """
        visited = self._dictionary.setdefault(key, [False])
        if visited[0]:
            return False
        visited[0] = True
        return True

class DictBasedProductSet:
    """
    This class can be constructed from a dictionary with keys in A.

    The resulting object is a set of AxB.

    Can be combined with RealHashDict or CanonicalKeyDict.

    >>> fruit_number_set = DictBasedProductSet(dict())
    >>> fruit_number_set.add('Apple', 1)
    True
    >>> fruit_number_set.add('Banana', 1)
    True
    >>> fruit_number_set.add('Apple', 2)
    True
    >>> fruit_number_set.add('Apple', 1) # Already added earlier
    False

    """

    def __init__(self, dictionary):
        self._dictionary = dictionary

    def add(self, key, element) -> bool:
        """
        Adds (key, element) in AxB to set.

        Returns true if (key, element) was not yet in the set.
        """
        
        elements = self._get_elements_for_key(key)
        if element in elements:
            return False
        else:
            elements.add(element)
            return True

    def _get_elements_for_key(self, key):
        """
        Could be used for potential future optimization:
        When going to a neighboring tetrahedron paired through
        a trivial face-pairing matrix (so it is a tetrahedron within
        the same copy of the fundamental polyhedron), we can keep
        a reference to the set returned by _get_elements_for_key
        and just add the neighboring tetrahedron to this set rather
        than doing the look-up from scratch.
        """
        return self._dictionary.setdefault(key, set())
