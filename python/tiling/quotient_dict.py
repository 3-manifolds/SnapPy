class QuotientDict:
    """
    Applies a function to canonize key before looking it up
    or setting it in the supplied dictionary.

    Note that the function can return a set of multiple
    canidadates for the canonical key. The look-ups in the
    dictionary will work correctly, as long as the true canonical
    key is among the candidates.

    # Modulo 5, canonical representative is 0, 1, ... 4.

    >>> d = QuotientDict({}, canonical_representative_function = lambda x : [x % 5])
    >>> d.setdefault(1, []).append('A')
    >>> d.setdefault(2, []).append('B')
    >>> d.setdefault(6, []).append('C')
    >>> d.setdefault(1, [])
    ['A', 'C']

    # Test with a fuzzy "Modulo 5" that always includes
    # the canonical representative in 0, 1, ... 4 but sometimes
    # an additional representative.

    >>> d = QuotientDict({}, canonical_representative_function = lambda x : x)
    >>> d.setdefault([1], []).append('A')
    >>> d.setdefault([2, 7], []).append('B')
    >>> d.setdefault([1,6], [])
    ['A']
    >>> d.setdefault([1,6], []).append('C')
    >>> d.setdefault([1,6], [])
    ['A', 'C']
    >>> d.setdefault([1], [])
    ['A', 'C']
    >>> d.setdefault([1, 11], [])
    ['A', 'C']
    >>> d.setdefault([2], [])
    ['B']
    """

    def __init__(self, dictionary, canonical_representative_function):
        self._dictionary = dictionary
        self._canonical_representative_function = canonical_representative_function

    # Not implemented because it was not needed so far:
    # get
    # __getitem__
    # __setitem__

    # Note that __setitem__ requires updating all entries in
    # self._dictionary.
    # The easiest way to do this is to wrap the value into
    # wrapped_value = [value], and assign to wrapped_value[0].
    # This would also allow self.setdefault(key, None).

    def setdefault(self, key, default):
        if default is None:
            raise Exception("Implementation reserved default = None for "
                            "internal purposes.")

        computed_keys = []
        for computed_key in self._canonical_representative_function(key):
            value = self._dictionary.get(computed_key)
            if not value is None:
                return value
            computed_keys.append(computed_key)

        for computed_key in computed_keys:
            self._dictionary[computed_key] = default

        return default
