from ..SnapPy import reduce_list_word, inverse_list_word

from typing import List

def simplify_geodesic_word(word : List[int]) -> List[int]:
    """
    Simplifies given word. It can change the word by a conjugate or
    invert it.

    More precisely, it cancels pairs of generator and inverse. First if
    they are next to each other and then if they are at opposite ends of
    the word. It then cyclically rotates the word and its inverse to take
    the lexicographically smallest when ordering the generators as
    1, -1, 2, -2, 3, -3, ...

    >>> simplify_geodesic_word([])
    []
    >>> simplify_geodesic_word([1])
    [1]
    >>> simplify_geodesic_word([1, 2])
    [1, 2]
    >>> simplify_geodesic_word([-1, 2, -2, 1])
    []
    >>> simplify_geodesic_word([-1, 3, 2, -2, 1])
    [3]
    >>> simplify_geodesic_word([-1, 3, 4, 2, -2, 1])
    [3, 4]
    >>> simplify_geodesic_word([-1, 4, 3, 2, -2, 1])
    [3, 4]
    >>> simplify_geodesic_word([-1, -3, -4, 5, 2, -2, 1])
    [3, -5, 4]
    """

    return (
        _rotate_and_optionally_invert(
            _cancel_conjugation(
                reduce_list_word(word))))

def _cancel_conjugation(word : List[int]) -> List[int]:
    """
    Cancels pairs of generator and inverse at opposite ends of
    given word.

    >>> _cancel_conjugation([])
    []
    >>> _cancel_conjugation([1])
    [1]
    >>> _cancel_conjugation([1, -1])
    []
    >>> _cancel_conjugation([1, 2, -1])
    [2]
    >>> _cancel_conjugation([1, 2, 3, -1])
    [2, 3]
    """
    n = len(word)
    for i in range(n // 2):
        k = n - i - 1
        if word[i] != -word[k]:
            return word[i:k+1]
    return word[n // 2 : (n + 1) // 2]

def _rotate_and_optionally_invert(word : List[int]) -> List[int]:
    """
    Rotate word and its inverse to pick lexicographically smallest when
    ordering the generators as 1, -1, 2, -2, 3, -3, ...

    >>> _rotate_and_optionally_invert([3, 4, 5])
    [3, 4, 5]
    >>> _rotate_and_optionally_invert([4, 5, 3])
    [3, 4, 5]
    >>> _rotate_and_optionally_invert([-3])
    [3]
    >>> _rotate_and_optionally_invert([5, 3, 4, -3])
    [3, 4, -3, 5]
    """
    n = len(word)

    if n == 0:
        return word

    return min(
        (candidates[i:] + candidates[:i]
         for candidates in [ word, inverse_list_word(word) ]
         for i in range(len(word))),
        key=lambda w: [ 2 * l if l > 0 else 2 * -l + 1
                        for l in w ])
