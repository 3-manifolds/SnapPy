# get_triangulation

split_filling_info = re.compile('(.*?)((?:\([0-9 .+-]+,[0-9 .+-]+\))*$)')
is_torus_bundle = re.compile('b([+-no])([+-])([lLrR]+)$')
is_link_complement1_pat = '(?P<crossings>[0-9]+)[\^](?P<components>[0-9]+)[_](?P<index>[0-9]+)$'
is_link_complement2_pat = '(?P<crossings>[0-9]+)[_](?P<index>[0-9]+)[\^](?P<components>[0-9]+)$'
is_link_complement3_pat = '[lL](?P<components>[0-9]{1})(?P<crossings>[0-9]{2})(?P<index>[0-9]+)$'
is_link_complement1 = re.compile(is_link_complement1_pat)
is_link_complement2 = re.compile(is_link_complement2_pat)
is_link_complement3 = re.compile(is_link_complement3_pat)
rolfsen_link_regexs = [is_link_complement1, is_link_complement2, is_link_complement3]
is_HT_knot = re.compile('(?P<crossings>[0-9]+)(?P<alternation>[an])(?P<index>[0-9]+)$')
is_braid_complement = re.compile('[Bb]raid[:]?(\[[0-9, \-]+\])$')
is_int_DT_exterior = re.compile('DT[:]? *(\[[0-9, \-\(\)]+\](?: *, *\[[01, ]+\])?)$')
is_alpha_DT_exterior = re.compile('DT[:\[] *([a-zA-Z]+(?:\.[01]+)?)[\]]?$')
twister_word = '\[*([abcABC_\d!*]*|[abcABC_\d!,]*)\]*'
is_twister_bundle = re.compile('Bundle\(S_\{(\d+),(\d+)\},'+twister_word+'\)')
is_twister_splitting = re.compile('Splitting\(S_\{(\d+),(\d+)\},'+twister_word+','+twister_word+'\)')
is_isosig = re.compile('([a-zA-Z0-9\+\-]+)$')
is_decorated_isosig = decorated_isosig.isosig_pattern

# Hooks so that global module can monkey patch in modified versions
# of the Triangulation and Manifold classes.

_triangulation_class = Triangulation
_manifold_class = Manifold

def bundle_from_string(desc):
    desc = desc.replace(' ', '')
    m = is_twister_bundle.match(desc)
    if m:
        g, n, monodromy = m.groups()
        g, n = int(g), int(n)
        monodromy = monodromy.replace(',', '*').replace('_', '')
        surface = twister.Surface( (g, n) )
        return surface.bundle(monodromy, return_type='string')
    
def splitting_from_string(desc):
    desc = desc.replace(' ', '')
    m = is_twister_splitting.match(desc)
    if m:
        g, n, gluing, handles = m.groups()
        g, n = int(g), int(n)
        gluing = gluing.replace(',', '*').replace('_', '')
        handles = handles.replace(',', '*').replace('_', '')
        surface = twister.Surface( (g, n) )        
        return surface.splitting(gluing, handles, return_type='string')

#Orientability.orientable = 0
rev_spec_dict = {(5, 0) : 'm',
                 (5, 1) : 'm',
                 (6, 0) : 's',
                 (7, 0) : 'v',
                 (6, 1) : 'x',
                 (7, 1) : 'y'}

triangulation_help =  """ 
    A %s is specified by a string, according to the
    conventions detailed in its docstring.  
    """

cdef c_Triangulation* triangulation_from_bytes(bytestring) except ? NULL:
    cdef c_Triangulation* c_triangulation = NULL
    cdef TerseTriangulation c_terse
    cdef int N=0, n=0, m=1
    byteseq = bytearray(bytestring)
    c_terse.num_tetrahedra = N = byteseq[0]
    c_terse.glues_to_old_tet = <Boolean*>malloc(2*N*sizeof(Boolean))
    c_terse.which_old_tet = <int*>malloc((N+1)*sizeof(int))
    c_terse.which_gluing = <Permutation*>malloc((N+1)*sizeof(Permutation))
    bit = 0
    for n from 0 <= n < 2*N:
        if byteseq[m] & (1 << bit):
            c_terse.glues_to_old_tet[n] = 1
        else:
            c_terse.glues_to_old_tet[n] = 0
        bit += 1
        if bit%8 == 0:
            m += 1
            bit = 0
    if bit:
        m += 1
    for n from 0 <= n < 1 + N:
        c_terse.which_old_tet[n] = int(byteseq[m])
        m += 1
    for n from 0 <= n < 1 + N:
        c_terse.which_gluing[n] = int(byteseq[m])
        m += 1
    c_triangulation = terse_to_tri(&c_terse)
    free(c_terse.glues_to_old_tet)
    free(c_terse.which_old_tet)
    free(c_terse.which_gluing)
    return c_triangulation
        
cdef int set_cusps(c_Triangulation* c_triangulation, fillings) except -1:
    if c_triangulation == NULL:
        return 0
    if len(fillings) > 0:
        num_cusps = get_num_cusps(c_triangulation) 
        if len(fillings) > num_cusps:
            raise ValueError('The number of fillings specified exceeds '
                             'the number of cusps.')
        for i in range(len(fillings)):
            meridian, longitude = fillings[i]
            is_complete = (meridian == 0 and longitude == 0)
            set_cusp_info(c_triangulation, <int>i,
                          is_complete,
                          Object2Real(meridian),
                          Object2Real(longitude))
    return 0

# Testing code for get_triangulation
def get_triangulation_tester():
    """
    >>> get_triangulation_tester()
    L13n9331(0,0)(0,0)(0,0) 16.64369585 Z + Z + Z
    m003(0,0) 2.02988321 Z/5 + Z
    m004(0,0) 2.02988321 Z
    v1205(2,3) 4.70744340 Z/40
    x012(0,0)(0,0) 3.54972978 Z/2 + Z
    y123(0,0) 5.02755480 Z
    L13n9331(3,4)(2,3)(2,1) 14.60215339 Z/53
    K7_1(0,0) 3.57388254 Z
    6_1(0,0) 3.16396323 Z
    5^2_1(3,4)(1,-2) 2.73300075 Z/3
    8^3_3(0,0)(0,0)(0,0) 8.96736085 Z + Z + Z
    4_1(0,0) 2.02988321 Z
    12n123(0,0) 18.15036328 Z
    16n1235(0,0) 21.29383093 Z
    b++RL(0,0) 2.02988321 Z
    b-+RRL(0,0) 2.40690959 Z/3 + Z
    b+-RL(0,0) 2.02988321 Z/5 + Z
    b--RRL(0,0) 2.40690959 Z/3 + Z
    Braid:[1, 2, -1, -2](0,0)(0,0) 4.05976643 Z + Z
    DT:[(8, 10, -14), (2, 6, 20), (-4, 22, 24, 12, 26, 18, 16)](0,0)(0,0)(0,0) 16.64369585 Z + Z + Z
    DT[4,6,2](0,0) 0.0 Z
    DT[mcccgdeGacjBklfmih](0,0)(0,0)(0,0) 16.64369585 Z + Z + Z
    DT:mcccgdeGacjBklfmih(0,0)(0,0)(0,0) 16.64369585 Z + Z + Z
    a0*B1(0,0) 2.02988321 Z
    b1*A0 a0*B1(1,0) 0.0 Z/2
    dLQbcccdxwb(0,0) 2.56897060 Z/3 + Z
    L13n9331(0,0)(0,0)(0,0) Z + Z + Z
    m003(0,0) Z/5 + Z
    m004(0,0) Z
    v1205(2,3) Z/40
    x012(0,0)(0,0) Z/2 + Z
    y123(0,0) Z
    L13n9331(3,4)(2,3)(2,1) Z/53
    K7_1(0,0) Z
    6_1(0,0) Z
    5^2_1(3,4)(1,-2) Z/3
    8^3_3(0,0)(0,0)(0,0) Z + Z + Z
    4_1(0,0) Z
    12n123(0,0) Z
    16n1235(0,0) Z
    b++RL(0,0) Z
    b-+RRL(0,0) Z/3 + Z
    b+-RL(0,0) Z/5 + Z
    b--RRL(0,0) Z/3 + Z
    Braid:[1, 2, -1, -2](0,0)(0,0) Z + Z
    DT:[(8, 10, -14), (2, 6, 20), (-4, 22, 24, 12, 26, 18, 16)](0,0)(0,0)(0,0) Z + Z + Z
    DT[4,6,2](0,0) Z
    DT[mcccgdeGacjBklfmih](0,0)(0,0)(0,0) Z + Z + Z
    DT:mcccgdeGacjBklfmih(0,0)(0,0)(0,0) Z + Z + Z
    a0*B1(0,0) Z
    b1*A0 a0*B1(1,0) Z/2
    dLQbcccdxwb(0,0) Z/3 + Z
    """

    M = database.HTLinkExteriors['L13n9331']
    specs = [M._to_string(), 'm003', 'm004', 'v1205(2,3)', 'x012', 'y123',
         'L13n9331(3,4)(2,3)(2,1)', 'K7_1', '6_1',
         '5^2_1(3,4)(1,-2)', '8_3^3', 'L104001', '12n123', '16n1235',
         'b++RL', 'b-+RRL', 'b+-RL', 'b--RRL',
         'Braid[1,2,-1,-2]', 'DT:'+repr(M.DT_code()), 'DT[4,6,2]',
         'DT['+M.DT_code(alpha=True) + ']',
         'DT:'+M.DT_code(alpha=True),
         'Bundle(S_{1,1}, [a0, B1])', 'Splitting(S_{1,0}, [b1, A0], [a0,B1])',
         'dLQbcccdxwb',
         ]

    for spec in specs:
        M = Manifold(spec)
        vol = M.volume()
        if abs(vol) < 0.1:
            vol = 0.0
        print M, vol, M.homology()

    for spec in specs:
        M = Triangulation(spec)
        print M, M.homology()

# Support for Hoste-Thistethwaite tables

# These dictionaries are used in accessing the tables.  The key is the
# number of crossings; the value is the number of knots with that many
# crossings.

Alternating_numbers = { 3:1, 4:1, 5:2, 6:3, 7:7, 8:18, 9:41, 10:123, 11:367,
                        12:1288, 13:4878, 14:19536, 15:85263, 16:379799 }

Nonalternating_numbers = { 8:3, 9:8, 10:42, 11:185, 12:888, 13:5110,
                           14:27436, 15:168030, 16:1008906 }

Alternating_offsets = {}
offset = 0
for i in range(3,17):
    Alternating_offsets[i] = offset
    offset +=  Alternating_numbers[i]
Num_Alternating = offset

Nonalternating_offsets = {}
offset = 0
for i in range(8,17):
    Nonalternating_offsets[i] = offset
    offset += Nonalternating_numbers[i]
Num_Nonalternating = offset

def extract_HT_knot(record, crossings, alternation):
    DT=[]
    size = (1+crossings)//2
    for byte in record[:size]:
        first_nybble = (byte & 0xf0) >> 4
        if first_nybble == 0: first_nybble = 16
        DT.append(2*first_nybble)
        second_nybble = byte & 0x0f
        if second_nybble == 0: second_nybble = 16
        DT.append(2*second_nybble)
    if alternation == 'n':
        signs = record[-2]<<8 | record[-1]
        mask = 0x8000
        for i in range(crossings):
            if (signs & (mask >> i)) == 0:
                DT[i] = -DT[i]
    return DT[:crossings]

def get_HT_knot_DT(crossings, alternation, index):
    size = (1 + crossings)//2
    index -= 1
    if ( alternation == 'a'
         and crossings in Alternating_numbers.keys()
         and 0 <= index < Alternating_numbers[crossings] ):
        offset = 8*(Alternating_offsets[crossings] +  index)
        Alternating_table.seek(offset)
        data = Alternating_table.read(size)
        record = struct.unpack('%dB'%size, data)
    elif ( alternation == 'n'
         and crossings in Nonalternating_numbers.keys()
         and 0 <= index < Nonalternating_numbers[crossings] ):
        offset = 10*(Nonalternating_offsets[crossings] +  index)
        Nonalternating_table.seek(offset)
        data = Nonalternating_table.read(size+2)
        record = struct.unpack('%dB'%(size+2), data)
    else:
        raise ValueError('You have specified a Hoste-Thistlethwaite '
                         'knot with an \n'
                         'inappropriate index or number of crossings.')

    DT = extract_HT_knot(record, crossings, alternation)
    return DT

def get_HT_knot_by_index(alternation, index, manifold_class):
    DT=[]
    crossings = 16
    if alternation == 'a':
        for i in range(3,17):
            if Alternating_offsets[i] > index:
                crossings = i-1
                break
        index_within_crossings = index - Alternating_offsets[crossings]
    elif alternation == 'n':
        for i in range(8,17):
            if Nonalternating_offsets[i] > index:
                crossings = i-1
                break
        index_within_crossings = index - Nonalternating_offsets[crossings]
    name = '%d'%crossings + alternation + '%d'%(index_within_crossings + 1)
    return manifold_class(name)

#   Iterators

class Census:
    """
    Base class for manifold Iterators/Sequences.
    """
    # subclasses redefine this
    length = 0

    def __init__(self, indices=(0,0,0)):
        myslice = slice(*indices)
        self.start, self.stop, self.step = myslice.indices(self.length)
        self.index = self.start

    def __iter__(self):
        return self

    def __contains__(self, item):
        raise NotImplementedError("This census does not support manifold lookup")

    def next(self):
        if self.index >= self.stop:
            raise StopIteration
        self.index = self.index + self.step
        return self[self.index-self.step]

    __next__ = next

    def __len__(self):
        return self.length
    
    # subclasses override this
    def __getitem__(self, n):
        pass

#  Cusped Census

Orientable_lengths = (301, 962, 3552, 12846, 301+962+3552+12846)
Nonorientable_lengths = (114, 259, 887, 114+259+887) 

five_tet_orientable = (
   3, 4, 6, 7, 9, 10, 11, 15, 16, 17, 19, 22, 23, 26, 27, 29, 30, 32, 33,
   34, 35, 36, 37, 38, 39, 40, 43, 44, 45, 46, 47, 49, 52, 53, 54, 55, 58,
   59, 60, 61, 62, 64, 66, 67, 69, 70, 71, 72, 73, 74, 76, 77, 78, 79, 80,
   81, 82, 83, 84, 85, 87, 89, 90, 93, 94, 95, 96, 98, 99, 100, 102, 103,
   104, 105, 106, 108, 110, 111, 112, 113, 115, 116, 117, 118, 119, 120,
   121, 122, 123, 125, 129, 130, 135, 136, 137, 139, 140, 141, 142, 143,
   144, 145, 146, 147, 148, 149, 150, 151, 154, 155, 157, 159, 160, 161,
   162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 175, 178,
   179, 180, 181, 182, 183, 184, 185, 186, 188, 189, 190, 191, 192, 193,
   194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 206, 207, 208, 209,
   210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223,
   224, 225, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 238, 239,
   240, 241, 242, 243, 244, 246, 247, 248, 249, 250, 251, 252, 253, 255,
   256, 257, 258, 259, 260, 261, 262, 263, 266, 267, 269, 270, 271, 272,
   273, 275, 276, 277, 278, 279, 280, 281, 282, 284, 285, 286, 287, 288,
   289, 290, 291, 292, 293, 294, 295, 296, 297, 299, 300, 302, 303, 304,
   305, 306, 307, 310, 312, 316, 318, 319, 320, 322, 326, 328, 329, 335,
   336, 337, 338, 339, 340, 342, 343, 345, 346, 349, 350, 351, 352, 356,
   357, 358, 359, 360, 361, 363, 364, 366, 367, 368, 369, 370, 371, 372,
   373, 374, 376, 378, 385, 388, 389, 390, 391, 392, 393, 395, 397, 400,
   401, 402, 403, 410, 412)

five_tet_nonorientable = (
   0, 1, 2, 5, 8, 12, 13, 14, 18, 20, 21, 24, 25, 28, 31, 41, 42, 48, 50,
   51, 56, 57, 63, 65, 68, 75, 86, 88, 91, 92, 97, 101, 107, 109, 114, 124,
   126, 127, 128, 131, 132, 133, 134, 138, 152, 153, 156, 158, 174, 176,
   177, 187, 204, 205, 226, 237, 245, 254, 264, 265, 268, 274, 283, 298,
   301, 308, 309, 311, 313, 314, 315, 317, 321, 323, 324, 325, 327, 330,
   331, 332, 333, 334, 341, 344, 347, 348, 353, 354, 355, 362, 365, 375,
   377, 379, 380, 381, 382, 383, 384, 386, 387, 394, 396, 398, 399, 404,
   405, 406, 407, 408, 409, 411, 413, 414)

class CuspedCensus(Census):
    """
    Base class for Iterators/Sequences for manifolds in the SnapPea
    Cusped Census.
    """
    five_length, six_length, seven_length, eight_length, length = Orientable_lengths
    orientability = Orientability.index('orientable')
    path = str(manifold_path)

    def __init__(self, indices=(0, length, 1)):
        Census.__init__(self, indices)
        self.Census_Morwen8 = tarfile.open(os.path.join(manifold_path, 'morwen8.tgz'), 'r:*')

    # Override
    def lookup(self, n):
        return five_tet_orientable[n]

    def __getitem__(self, n):
        cdef c_Triangulation* c_triangulation
        cdef Manifold result
        if isinstance(n, slice):
            return self.__class__(n.indices(self.length))
        if n < 0:
            n = self.length - n
        if n < self.five_length:
            num_tet = 5
            census_index = self.lookup(n)
        elif n - self.five_length < self.six_length:
            num_tet = 6
            census_index = n - self.five_length
        elif n - self.five_length - self.six_length < self.seven_length:
            num_tet = 7
            census_index = n - self.five_length - self.six_length
        elif (self.orientability == Orientability.index('orientable')
              and (n - self.five_length - self.six_length -
                   self.seven_length < self.eight_length)):
              census_index = (n - self.five_length - self.six_length
                              - self.seven_length)
              # Make this work without passing the spec to Manifold()
              num = repr(census_index)
              spec =  "t" + "0"*(5 - len(num)) + num
              tarpath = "morwen8/" + spec
              try:
                  filedata = self.Census_Morwen8.extractfile(tarpath).read()
                  c_triangulation = read_triangulation_from_string(filedata)
              except: 
                  raise IOError('The Morwen 8 tetrahedra manifold %s '
                                'was not found.'% spec)
              result = Manifold(spec='empty')
              result.set_c_triangulation(c_triangulation)
              return result              
              ###return Manifold('t%d' % census_index)
        else:
            raise IndexError('Index is out of range.')
        c_triangulation = GetCuspedCensusManifold(
            self.path, num_tet, self.orientability, census_index)
        if c_triangulation == NULL:
            print(num_tet, census_index)
            raise RuntimeError('SnapPea failed to read the census manifold.')
        result = Manifold(spec='empty')
        result.set_c_triangulation(c_triangulation)
        return result

class ObsOrientableCuspedCensus(CuspedCensus):
    """
    Obsolete.
    """

class ObsNonorientableCuspedCensus(CuspedCensus):
    """
    Obsolete.
    """
    five_length, six_length, seven_length, length = Nonorientable_lengths
    orientability = Orientability.index('nonorientable')

    def __init__(self, indices=(0, length, 1)):
        Census.__init__(self, indices)

    def lookup(self, n):
        return five_tet_nonorientable[n]

# Closed Census (Obsolete)

class ObsOrientableClosedCensus(Census):
    """
    Obsolete.
    """
    data = None
    orientability = Orientability.index('orientable')
    path = str(manifold_path)

    def __init__(self, indices=(0,11031,1)):
        if ObsOrientableClosedCensus.data is None:
            datafile = os.path.join(closed_census_directory,
                                    'ClosedOrientableDistinct.txt')
            closed_orientable = open(datafile)
            ObsOrientableClosedCensus.data = closed_orientable.readlines()
            closed_orientable.close()
        self.length = len(ObsOrientableClosedCensus.data)
        Census.__init__(self, indices)

    def __getitem__(self,n):
        cdef c_Triangulation* c_triangulation
        cdef Manifold result
        if isinstance(n, slice):
            return self.__class__(n.indices(self.length))
        volume, num_tet, index, m, l = ObsOrientableClosedCensus.data[n].split()
        c_triangulation = GetCuspedCensusManifold(
            self.path, int(num_tet), self.orientability, int(index))
        if c_triangulation == NULL:
            print(num_tet, index)
            raise RuntimeError('SnapPea failed to read the census manifold.')
        result = Triangulation(spec='empty')
        result.set_c_triangulation(c_triangulation)
        result.dehn_fill(( int(m),int(l)) )
        return result.with_hyperbolic_structure()

class ObsNonorientableClosedCensus(Census):
    """
    Obsolete.
    """
    data = None
    orientability = Orientability.index('nonorientable')
    path = str(manifold_path)

    def __init__(self, indices=(0,17,1)):
        if ObsNonorientableClosedCensus.data is None:
            datafile = os.path.join(closed_census_directory,
                                    'ClosedNonorientableDistinct.txt')
            closed_nonorientable = open(datafile)
            ObsNonorientableClosedCensus.data = closed_nonorientable.readlines()
            closed_nonorientable.close()
        self.length = len(ObsNonorientableClosedCensus.data)
        Census.__init__(self, indices)

    def __getitem__(self,n):
        cdef c_Triangulation* c_triangulation
        cdef Manifold result
        if isinstance(n, slice):
            return self.__class__(n.indices(self.length))
        volume, num_tet, index, m, l = ObsNonorientableClosedCensus.data[n].split()
        c_triangulation = GetCuspedCensusManifold(
            self.path, int(num_tet), self.orientability, int(index))
        if c_triangulation == NULL:
            print(num_tet, index)
            raise RuntimeError('SnapPea failed to read the census manifold.')
        result = Triangulation(spec='empty')
        result.set_c_triangulation(c_triangulation)
        result.dehn_fill( (int(m),int(l)) )
        return result.with_hyperbolic_structure()

# Knot tables

class KnotExteriors(Census):
    """
    Base class for Iterators/Sequences for knots from the
    Hoste-Thistlethwaite tables.
    """
    length = sum(Alternating_numbers.values())
    alternation = 'a'

    def __init__(self, indices=(0, sum(Alternating_numbers.values()), 1)):
        Census.__init__(self, indices)

    def __getitem__(self, n):
        if isinstance(n, slice):
            return self.__class__(n.indices(self.length))
        else:
            return get_HT_knot_by_index(self.alternation, n, _manifold_class)

class AlternatingKnotExteriors(KnotExteriors):
    """
    Iterator/Sequence for Alternating knot exteriors from the
    Hoste-Thistlethwaite tables.   Goes through 16 crossings. 
    """

class NonalternatingKnotExteriors(KnotExteriors):
    """
    Iterator/Sequence for nonAlternating knot exteriors from the
    Hoste-Thistlethwaite tables.  Goes through 16 crossings. 
    """
    length = sum(Nonalternating_numbers.values())
    alternation = 'n'

    def __init__(self, indices=(0, sum(Nonalternating_numbers.values()), 1)):
        Census.__init__(self, indices)

census_knot_numbers = [0, 0, 1, 2, 4, 22, 43, 129]

class ObsCensusKnots(Census):
    """
    Obsolete
    """
    length = sum(census_knot_numbers)

    def __init__(self, indices=(0, sum(census_knot_numbers), 1)):
        Census.__init__(self, indices)
        self.Census_Knots = tarfile.open(census_knot_archive, 'r:*')

    def __repr__(self):
        return 'Knots in S^3 which appear in the SnapPea Census'
    
    def __getitem__(self, n):
        if isinstance(n, slice):
            return self.__class__(n.indices(self.length))
        else:
            total = 0
            for m in range(2,8):
                if total + census_knot_numbers[m] <= n :
                    total += census_knot_numbers[m]
                    continue
                else:
                    name = 'K%s_%s'%(m, n - total + 1)
                    break
            if name:
                tarpath =  'CensusKnots/%s'%name
                try:
                    filedata = self.Census_Knots.extractfile(tarpath).read()
                    c_triangulation = read_triangulation_from_string(filedata)
                except: 
                    raise IOError, "The census knot %s was not found."%name
                result =  Triangulation('empty')
                result.set_c_triangulation(c_triangulation)
                result.set_name(name)
                return result.with_hyperbolic_structure()
            else:
                raise IndexError('There are only 201 census knots.')

class ObsLinkExteriors(Census):
    """
    Obsolete.
    """
    # num_links[component][crossings] is the number of links with
    # specified number of components and crossings.
    num_links= [ [], [0, 0, 0, 1, 1, 2, 3, 7, 21, 49, 166, 552],
                 [0, 0, 0, 0, 1, 1, 3, 8, 16, 61, 183],
                 [0, 0, 0, 0, 0, 0, 3, 1, 10, 21, 74],
                 [0, 0, 0, 0, 0, 0, 0, 0, 3, 1, 24],
                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3]]

    max_crossings = 11

    def __init__(self, components, indices=(0,10000,1)):
         self.Christy_links = tarfile.open(link_archive, 'r:*')

         if not (1 <= components < len(self.num_links) ):
            raise IndexError('SnapPy has no data on links with '
                             '%s components.' % components)

         self.components = components

         self.length = sum(self.num_links[components])

         Census.__init__(self, (indices[0],
                                min(self.length, indices[1]),
                                indices[2]))

    def __repr__(self):
        return ('Christy census of %s-component link complements '
                'in S^3' % self.components)

    def __getitem__(self,j):
        if isinstance(j, slice):
            return self.__class__(j.indices(self.length))
        so_far = 0
        for k in range(self.max_crossings + 1):
            n =  self.num_links[self.components][k]
            so_far = so_far + n
            if so_far > j:
                l = j - so_far + n + 1
                filename = 'L%d%.2d%.3d' % (self.components, k, l)
                if self.components > 1:
                    name = "%d^%d_%d" % (k, self.components, l)
                else:        
                    name = "%d_%d" % (k,  l)
                tarpath =  'ChristyLinks/%s'%filename
                try:
                    filedata = self.Christy_links.extractfile(tarpath).read()
                    c_triangulation = read_triangulation_from_string(filedata)
                except: 
                    raise IOError('The link complement %s was not '
                                  'found.'%filename)
                result =  Triangulation('empty')
                result.set_c_triangulation(c_triangulation)
                result.set_name(name)
                return result.with_hyperbolic_structure()              

#----------------------------------------------------------------
#
#  Morwen's link table (data taken from Snap).  Now obsolete.
#
#----------------------------------------------------------------

class MorwenLinks:
    def __init__(self, num_components, num_crossings = None):
        raise ValueError(
            'The MorwenLinks class has been replaced by HTLinkExteriors. '
            'Please refer to the documentation for HTLinkExteriors.'
            )

def left_pad_string(str,  length, pad=' '):
    return pad*(length - len(str)) + str

class ObsMorwenLinks(Census):
    """
    Obsolete.
    """
    def __init__(self, num_components, num_crossings = None):
        if num_components < 0:
            return
        
        crossings = range(4, 15) if num_crossings == None else [num_crossings]
        files = []
        for c in crossings:
            n = left_pad_string('%d' % c, 2, '0')
            files.append('hyperbolic_data_%sa.gz' % n)
            if c > 5:
                files.append('hyperbolic_data_%sn.gz' % n)

        self.files = files

        self.DT_codes = []
        second_char = string.ascii_lowercase[num_components-1]
        for file in files:
            data = gzip.open(morwen_link_directory + os.sep + file).readlines()
            for line in data:
                strline = line.decode()
                if strline[1] == second_char:
                    self.DT_codes.append(strline.strip())

        self.length =  len(self.DT_codes)
        Census.__init__(self, indices=(0, len(self.DT_codes), 1))

    def __getitem__(self, n):
        if isinstance(n, slice):
            SC = self.__class__(-1)
            SC.DT_codes = self.DT_codes[n]
            SC.length = len(SC.DT_codes)
            Census.__init__(SC, indices=(0, len(SC.DT_codes), 1))
            return SC
        else:
            name = str('DT:%s'%self.DT_codes[n])
            result = Manifold(name)
            result.set_name(name)
            return result

# Creating fibered manifolds from braids

cdef c_Triangulation* get_fibered_manifold_associated_to_braid(num_strands,
                                                               braid_word):
    if num_strands < 2:
        raise ValueError('Braids must have at least 2 strands.')
    allowed_letters = ([x for x in range(1,num_strands)] +
                       [x for x in range(-num_strands+1, 0)])
    if False in [b in allowed_letters for b in braid_word]:
        raise ValueError('The braid word is invalid.')

    cdef int* word
    cdef c_Triangulation* c_triangulation

    n = len(braid_word)
    word = <int*>malloc(n*sizeof(int))
    for  i, w in enumerate(braid_word):
        word[i] = w
    c_triangulation = fibered_manifold_associated_to_braid(num_strands, n, word)
    free(word)
    name = to_byte_str('Braid:' + repr(braid_word))
    set_triangulation_name(c_triangulation,name)
    return c_triangulation

# Link Editor support

strandtype = {'X': KLPStrandX,     'Y': KLPStrandY}
signtype =   {'R': KLPHalfTwistCL, 'L': KLPHalfTwistCCL}

cdef c_Triangulation* get_triangulation_from_PythonKLP(pythonklp, remove_finite_vertices=True) except *:
    cdef KLPProjection P
    cdef c_Triangulation  *c_triangulation

    P.num_crossings, P.num_free_loops, P.num_components, Pycrossings = pythonklp
    P.crossings = <KLPCrossing *>malloc(P.num_crossings * sizeof(KLPCrossing))
    if P.crossings == NULL:
        raise RuntimeError('Could not allocate a crossing table.')

    for i from 0 <= i < P.num_crossings:
        cr_dict = Pycrossings[i]
        neighbor = cr_dict['Xbackward_neighbor']
        P.crossings[i].neighbor[<int>KLPStrandX][<int>KLPBackward] = &P.crossings[neighbor] 
        neighbor =  cr_dict['Xforward_neighbor']
        P.crossings[i].neighbor[<int>KLPStrandX][<int>KLPForward] = &P.crossings[neighbor]
        neighbor = cr_dict['Ybackward_neighbor']
        P.crossings[i].neighbor[<int>KLPStrandY][<int>KLPBackward] = &P.crossings[neighbor]
        neighbor = cr_dict['Yforward_neighbor']
        P.crossings[i].neighbor[<int>KLPStrandY][<int>KLPForward] = &P.crossings[neighbor]

        strand = cr_dict['Xbackward_strand']
        P.crossings[i].strand[<int>KLPStrandX][<int>KLPBackward] = strandtype[strand]
        strand = cr_dict['Xforward_strand']
        P.crossings[i].strand[<int>KLPStrandX][<int>KLPForward] = strandtype[strand]
        strand = cr_dict['Ybackward_strand']
        P.crossings[i].strand[<int>KLPStrandY][<int>KLPBackward] = strandtype[strand]
        strand = cr_dict['Yforward_strand']
        P.crossings[i].strand[<int>KLPStrandY][<int>KLPForward] = strandtype[strand]

        sign = cr_dict['sign']
        P.crossings[i].handedness = signtype[sign[0]]
        P.crossings[i].component[<int>KLPStrandX] = cr_dict['Xcomponent']
        P.crossings[i].component[<int>KLPStrandY] = cr_dict['Ycomponent']

    c_triangulation = triangulate_link_complement(&P, remove_finite_vertices);
    free(P.crossings)

    if c_triangulation == NULL:
        raise RuntimeError('Could not create the triangulation.')

    # The triangulation must have a name or SnapPea will segfault when
    #   trying to copy the triangulation.
    tri_name = to_byte_str('unnamed link')
    set_triangulation_name(c_triangulation, tri_name)
    return c_triangulation

