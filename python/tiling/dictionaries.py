

class MultiRealDict:
    """
    Dictionary where key is a real number.

    The real number can be represented as floating point number or interval.

    Returns all numbers within epsilon or for overlapping intervals.

    Implemented by interval tree or computing integer hash by dividing by epsilon.
    """

    def __init__(self, epsilon):
        """
        epsilon used depending on implementation.
        """
    
    # Key is a real number (represented as floating point number or interval)
    def add(self, key, value):
        """
        Adds key, value.
        """

    def get(self, key):
        """
        Gets all [ value' : (key', value') where key and key' are close ].
        """

class HashDict:
    """
    Keys together with hash and predicate to decide whether two
    keys are equal.
    """

    def __init__(self, multi_dict, equality_predicate, hash_function):
        pass

    # Ordinary dictionary API: get, [], setdefault

class MultiKeyDict:
    """
    If storing a key, stores value in each key returned by
    canonical_keys.
    When reading, reads for each key returned by canonical_keys.
    """

    def __init__(self, dictionary, canonical_keys):
        pass
        
    # Ordinary dictionary API: get, [], setdefault

class ProductSetViaDictionary:
    """
    Stores key -> set.
    """

    def __init__(self, dictionary):
        pass
    
    # in operator takes (key, element)
    # add takes (key, element)

    # Optimization: get set for key - when going to a tetrahedron within the same copy of the fundamental polyhedron, we can reuse that set.

class LiftedTetrahedronSet:
    def __init__(self, dictionary, base_point):
        pass

    # Ordinary set API
    # Takes (lifted_tetrahedron.matrix * base_point, lifted_tetrahedron.tetrahedron.Index) for ProductSetViaDictionary.

