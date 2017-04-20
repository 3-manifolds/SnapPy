
# Fundamental Groups

Alphabet = '$abcdefghijklmnopqrstuvwxyzZYXWVUTSRQPONMLKJIHGFEDCBA'

# Helper functions for manipulating fg. words

def inverse_word(word):
    parts = list(word.swapcase())
    parts.reverse()
    return ''.join(parts)

reduce_word_regexp = re.compile('|'.join([x + x.swapcase()
                                          for x in string.ascii_letters]))

def reduce_word(word):
    """
    Cancels inverse generators.
    """
    ans, progress = word, 1

    while progress:
        ans, progress =  reduce_word_regexp.subn('', ans)
    return ans

def format_word(word, verbose_form):
    if not verbose_form:
        return word
    if re.search('\d', word):
        letters = re.findall('([xX]\d+)', word)
    else:
        letters = list(word)
    return '*'.join([a if a[0].islower() else a.lower() + '^-1' for a in letters])

cdef class CFundamentalGroup(object):
    cdef c_GroupPresentation *c_group_presentation
    cdef c_Triangulation *c_triangulation
    cdef readonly num_cusps

    cdef int_to_gen_string(self, int g):
        if self.num_generators() <=26:
            return Alphabet[g]
        else:
            ans = 'x' if g > 0 else 'X'
            return ans + '%d' % abs(g)

    cdef c_word_as_int_list(self, int *word):
        cdef int n = 0
        word_list = []
        while word[n] != 0:
            word_list.append(word[n])
            n += 1
        return word_list

    cdef c_word_as_string(self, int *word):
        cdef int n = 0
        word_list = []
        while word[n] != 0:
            word_list.append(self.int_to_gen_string(word[n]))
            n += 1
        return ''.join(word_list)

    cdef int *c_word_from_list(self, word_list):
        cdef int *c_word
        cdef int length, size, n
        length = <int>len(word_list)
        size = sizeof(int)*(1+length)
        c_word = <int *>malloc(size)
        for n from 0 <= n < length:
            c_word[n] = word_list[n]
        c_word[length] = 0
        return c_word

    def __cinit__(self, Triangulation triangulation,
                  simplify_presentation = True,
                  fillings_may_affect_generators = True,
                  minimize_number_of_generators = True,
                  try_hard_to_shorten_relators = True):
        if triangulation.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')
        copy_triangulation(triangulation.c_triangulation,
                           &self.c_triangulation)
        self.c_group_presentation = fundamental_group(
            self.c_triangulation,
            simplify_presentation,
            fillings_may_affect_generators,
            minimize_number_of_generators,
            try_hard_to_shorten_relators)
        self.num_cusps = triangulation.num_cusps()

    def __dealloc__(self):
        free_triangulation(self.c_triangulation)
        free_group_presentation(self.c_group_presentation)

    def __repr__(self):
        return 'Generators:\n   %s\nRelators:\n   %s'%(
            ','.join(self.generators()),
            '\n   '.join(self.relators()))

    def _word_as_list(self, word):
        if not isinstance(word, basestring):
            raise TypeError('Words must be represented '
                            'as Python strings.')
        word_list = []
        generators = self.generators()
        if len(generators) > 26:
            word = re.findall('[xX]\d+', word)
        for letter in word:
            try:
                if letter[0].islower():
                    word_list.append(1 + generators.index(letter))
                else:
                    word_list.append(-1 - generators.index(letter.lower()))
            except ValueError:
                raise RuntimeError('The word contains a non-generator.')
        return word_list

    def num_generators(self):
        """
        Return the number of generators for the presentation.
        """
        return fg_get_num_generators(self.c_group_presentation)

    def num_relators(self):
        """
        Return the number of generators for the presentation.
        """
        return fg_get_num_relations(self.c_group_presentation)
                            
    def num_original_generators(self):
        """
        Return the number of geometric generators (before simplification).
        """
        return fg_get_num_orig_gens(self.c_group_presentation)

    def original_generators(self, verbose_form=False):
        """
        Return the original geometric generators (before
        simplification) in terms of the current generators.
        """
        cdef int n
        cdef int *gen
        orig_gens = []
        num_orig_gens = fg_get_num_orig_gens(self.c_group_presentation)
        for n from 0 <= n < num_orig_gens:
            gen = fg_get_original_generator(self.c_group_presentation, n)
            word = format_word(self.c_word_as_string(gen), verbose_form)
            orig_gens.append(word)
            fg_free_relation(gen)
        return orig_gens

    def generators_in_originals(self, verbose_form=False, raw_form =False):
        """
        Return the current generators in terms of the original
        geometric generators (before simplification).

        If the flag "raw_form" is set to True, it returns a sequence of
        instructions for expressing the current generators in terms of
        the orignal ones.  This is sometimes much more concise, though
        the format is somewhat obscure.  See the source code of this
        function in SnapPy.pyx for details. 
        """
        moves = self._word_moves()
        if raw_form:
            return moves

        n = self.num_original_generators()
        if n > 26:
            raise ValueError('Too many generators.')

        words = [None] + list(string.ascii_letters[:n])

        while len(moves) > 0:
            a = moves.pop(0)
            if a >= len(words): # new generator added
                n = moves.index(a)  # end symbol location
                # word is the expression of the new generator in terms
                # of the old ones
                word, moves = moves[:n], moves[n+1:]
                words.append( reduce_word(''.join(
                    [words[g] if g > 0 else inverse_word(words[-g])
                     for g in word]
                    )))
            else:
                b = moves.pop(0)
                if a == b:  # generator removed
                    words[a] = words[-1]
                    words = words[:-1]
                elif a == -b: # invert generator
                    words[a] = inverse_word(words[a])
                else: #handle slide
                    A, B = words[abs(a)], words[abs(b)]
                    if a*b < 0:
                        B = inverse_word(B)
                    words[abs(a)] = reduce_word(  A+B if a > 0 else B+A ) 

        return [format_word( w, verbose_form) for w in words[1:]]

    def _word_moves(self):
        cdef int *c_moves
        c_moves = fg_get_word_moves(self.c_group_presentation)
        moves = self.c_word_as_int_list(c_moves)
        fg_free_relation(c_moves)
        return moves
        
    def generators(self):
        """
        Return the letters representing the generators in the presentation.
        """
        n = self.num_generators()
        return [ self.int_to_gen_string(i) for i in range(1, 1+n) ]

    def relators(self, verbose_form = False, as_int_list = False):
        """
        Return a list of words representing the relators in the presentation.

        If the optional argument verbose_form is True, then the
        relator is returned in the form "a*b*a^-1*b^-1" instead of "abAB".  
        """
        cdef int n
        cdef int *relation
        relation_list = []
        num_relations = fg_get_num_relations(self.c_group_presentation)
        for n from 0 <= n < num_relations:
            relation = fg_get_relation(self.c_group_presentation, n)
            if as_int_list:
                word = self.c_word_as_int_list(relation)
            else:
                word = format_word(self.c_word_as_string(relation), verbose_form)
            relation_list.append(word)
            fg_free_relation(relation)
        return relation_list

    def meridian(self, int which_cusp=0, as_int_list = False):
        """
        Returns a word representing a conjugate of the current
        meridian for the given cusp.  Guaranteed to commute with the
        longitude for the same cusp.

        >>> G = Manifold('m125').fundamental_group()
        >>> G.meridian(0)
        'aaba'
        >>> G.meridian(-1)  # The last cusp
        'baaba'
        """
        try:
            which_cusp = range(self.num_cusps)[which_cusp]
        except IndexError:
            raise IndexError('The specified cusp (%s) does not '
                             'exist.'%which_cusp)
        if as_int_list:
            return self.c_word_as_int_list(
               fg_get_meridian(self.c_group_presentation, which_cusp))
        else:
            return self.c_word_as_string(
               fg_get_meridian(self.c_group_presentation, which_cusp))

    def longitude(self, int which_cusp=0, as_int_list = False):
        """
        Returns a word representing a conjugate of the current
        longitude for the given cusp.  Guaranteed to commute with the
        meridian for the same cusp.  Note: for Klein bottle cusps,
        the longitude must be defined carefully.

        >>> G = Manifold('m004').fundamental_group()
        >>> G.longitude(0)
        'aBAbABab'
        >>> G.longitude()   # shortcut for the above.  
        'aBAbABab'
        """
        try:
            which_cusp = range(self.num_cusps)[which_cusp]
        except IndexError:
            raise IndexError('The specified cusp (%s) does not '
                             'exist.'%which_cusp)

        if as_int_list:
            return self.c_word_as_int_list(
               fg_get_longitude(self.c_group_presentation, which_cusp))
        else:
            return self.c_word_as_string(
               fg_get_longitude(self.c_group_presentation, which_cusp))

    def peripheral_curves(self, as_int_list = False):
        """
        Returns a list of meridian-longitude pairs for all cusps.

        >>> G = Manifold('m125').fundamental_group()
        >>> G.peripheral_curves()
        [('aaba', 'abb'), ('baaba', 'Ba')]
        """
        return [ (self.meridian(n, as_int_list),
                  self.longitude(n, as_int_list))
                 for n in range(self.num_cusps) ]

    def magma_string(self):
        """
        Returns a string which will define this group within MAGMA.
        """
        return ('Group<' + ','.join(self.generators()) + '|' +
                ', '.join(self.relators(verbose_form = True)) + '>')

    def gap_string(self):
        """
        Returns a string which will define this group within GAP.
        """
        gens = ', '.join(self.generators())
        gen_names = ', '.join(['"' + x + '"' for x in self.generators()])
        relators = ', '.join(self.relators(verbose_form = True))
        assignments = ''.join(
            ['%s := F.%d; ' % (x, i+1)
             for (i, x) in enumerate(self.generators())]
            )
        return ('CallFuncList(function() local F, %s; '
                'F := FreeGroup(%s); %s  return F/[%s]; end,[])'
                % (gens, gen_names, assignments, relators)
                )

    def _gap_init_(self):
        return self.gap_string()

    def _magma_init_(self, magma):
        return self.magma_string()

    def sage(self):
        """
        Returns the corresponding Sage FinitelyPresentedGroup
        """
        if not _within_sage:
            raise RuntimeError("Not within Sage")
        F = FreeGroup(self.generators())
        rels = [F(R) for R in self.relators(as_int_list=True)]
        return F/rels

class FundamentalGroup(CFundamentalGroup):
    """
    A FundamentalGroup represents a presentation of the fundamental
    group of a SnapPea Triangulation.  Group elements are described as
    words in the generators a,b,..., where the inverse of a is denoted
    A.  Words are represented by Python strings (and the concatenation
    operator is named "+", according to Python conventions).

    Instantiate as T.fundamental_group(), where T is a Triangulation.

    Methods:
        num_generators() --> number of generators
        num_relators()   --> number of relators
        generators()     --> list of generators
        relators()       --> list of relators
        meridian(n)      --> word representing the meridian on cusp #n
        longitude(n)     --> word representing the longitude on cusp #n
    """

if _within_sage:
    FundamentalGroup.__bases__ += (sage.structure.sage_object.SageObject,)

# Holonomy Groups
cdef class CHolonomyGroup(CFundamentalGroup):
    def _matrices(self, word):
        """
        Returns (M,O,L) where M = SL2C(word), O = O31(word), and L is
        the complex length.
        """
        cdef MoebiusTransformation M 
        cdef O31Matrix O
        cdef int *c_word
        cdef c_FuncResult result
        cdef int i, j
        word_list = self._word_as_list(word)
        c_word = self.c_word_from_list(word_list)
        result = fg_word_to_matrix(self.c_group_presentation, c_word, O, &M)
        if result == 0:
            sl2 = matrix(
                [[self._number_(Complex2Number(M.matrix[i][j]))
                  for j in range(2)] for i in range(2)] )
            o31 = matrix(
                [[self._number_(Real2Number(<Real>O[i][j]))
                  for j in range(4)] for i in range(4)] )
            L = self._number_(Complex2Number(complex_length_mt(&M)))
            return sl2, o31, L
        else:
            return None

    def SL2C(self, word):
        """
        Return the image of the element represented by the input word
        under some SL(2,C) representation that lifts the holonomy
        representation.  Note: the choice of lift is not guaranteed to
        vary continuously when filling coefficients are changed.
        """
        return self._matrices(word)[0]

    def O31(self, word):
        """
        Return the image of the element represented by the input word
        under the holonomy representation, where Isom(H^3) is
        identified with SO(3,1).
        """
        return self._matrices(word)[1]

    def complex_length(self, word):
        """
        Return the complex length of the isometry represented by the
        input word.
        """
        return self._matrices(word)[2]

class HolonomyGroup(CHolonomyGroup):
    """
    A HolonomyGroup is a FundamentalGroup with added structure
    consisting of a holonomy representation into O(3,1), and an
    arbitrarily chosen lift of the holonomy representation to SL(2,C).
    The holonomy is determined by the shapes of the tetrahedra, so a
    HolonomyGroup is associated to a Manifold, while a Triangulation
    only has a FundamentalGroup.  Methods are provided to evaluate the
    representations on a group element.

    A FundamentalGroup represents a presentation of the fundamental
    group of a SnapPea Triangulation.  Group elements are described as
    words in the generators a,b,..., where the inverse of a is denoted
    A.  Words are represented by python strings (and the concatenation
    operator is named '+', according to Python conventions).

    Instantiate via M.fundamental_group(), where M is a Manifold.
    """
    @staticmethod
    def _number_(number):
        return number

    @classmethod
    def use_field_conversion(cls, func):
        cls._number_ = staticmethod(func)

if _within_sage:
    HolonomyGroup.__bases__ += (sage.structure.sage_object.SageObject,)
    HolonomyGroup.use_field_conversion(lambda n : n.sage())
