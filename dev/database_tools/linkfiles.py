# These classes are extracted from PLink, and simplified by removing all
# references to Tkinter 
from math import sqrt
DT_alphabet = '_abcdefghijklmnopqrstuvwxyzZYXWVUTSRQPONMLKJIHGFEDCBA'

class Vertex:
    """
    A vertex in a PL link diagram.
    """

    def __init__(self, x, y):
        self.x, self.y = int(x), int(y)
        self.in_arrow = None
        self.out_arrow = None

    def __repr__(self):
        return '(%d,%d)'%(self.x, self.y)

    def __eq__(self, other):
        """
        Vertices are equivalent if they are sufficiently close.
        Use the "is" operator to test if they are identical.
        """
        return abs(self.x - other.x) + abs(self.y - other.y) < Vertex.epsilon

    def point(self):
        return self.x, self.y

    def is_endpoint(self):

        return self.in_arrow == None or self.out_arrow == None
    
    def is_isolated(self):
        return self.in_arrow == None and self.out_arrow == None

    def reverse(self):
        self.in_arrow, self.out_arrow = self.out_arrow, self.in_arrow

    def update_arrows(self):
        if self.in_arrow: self.in_arrow.vectorize()
        if self.out_arrow: self.out_arrow.vectorize()

    def erase(self):
        """
        Prepare the vertex for the garbage collector.
        """
        self.in_arrow = None
        self.out_arrow = None

class Arrow:
    """
    An arrow in a PL link diagram.
    """
    
    def __init__(self, start, end):
        self.start, self.end = start, end
        self.lines = []
        self.cross_params = []
        if self.start != self.end:
            self.start.out_arrow = self
            self.end.in_arrow = self
            self.vectorize()
        
    def __repr__(self):
        return '%s-->%s'%(self.start, self.end)

    def __xor__(self, other):
        """
        Returns the barycentric coordinate at which self crosses other.
        """
        D = float(other.dx*self.dy - self.dx*other.dy)
        if D == 0:
            return None
        xx = other.start.x - self.start.x
        yy = other.start.y - self.start.y
        s = (yy*self.dx - xx*self.dy)/D
        t = (yy*other.dx - xx*other.dy)/D
        if 0 < s < 1 and 0 < t < 1:
            return t
        else:
            return None

    def vectorize(self):
        self.dx = float(self.end.x - self.start.x)
        self.dy = float(self.end.y - self.start.y)
        self.length = sqrt(self.dx*self.dx + self.dy*self.dy) 

    def reverse(self, crossings=[]):
        self.end, self.start = self.start, self.end
        self.vectorize()

    def set_start(self, vertex, crossings=[]):
        self.start = vertex
        if self.end:
            self.vectorize()

    def set_end(self, vertex, crossings=[]):
        self.end = vertex
        if self.start:
            self.vectorize()
            
    def erase(self):
        """
        Prepare the arrow for the garbage collector.
        """
        self.start = None
        self.end = None

class Crossing:
    """
    A pair of crossing arrows in a PL link diagram.
    """
    def __init__(self, over, under):
        self.over = over
        self.under = under
        self.KLP = {}    # See the SnapPea file link_projection.h
        self.hit1 = None # Markings for computing DT codes
        self.hit2 = None
        self.comp1 = None
        self.comp2 = None
        self.locate()

    def __repr__(self):
        return '%s over %s at (%d,%d)'%(self.over, self.under, self.x, self.y)

    def __hash__(self):
        # Since we redefine __eq__ we need to define __hash__
        return id(self)

    def __eq__(self, other):
        """
        Crossings are equivalent if they involve the same arrows.
        """
        if self.over in other and self.under in other:
            return True
        else:
            return False
        
    def __contains__(self, arrow):
        if arrow == None or arrow == self.over or arrow == self.under:
            return True
        else:
            return False
        
    def locate(self):
        t = self.over ^ self.under
        if t:
            self.x = int(self.over.start.x + t*self.over.dx)
            self.y = int(self.over.start.y + t*self.over.dy)
        else:
            self.x, self.y = None, None

    def sign(self):
        try:
            D = self.under.dx*self.over.dy - self.under.dy*self.over.dx
            if D > 0: return 'RH'
            if D < 0: return 'LH'
        except:
            return 0

    def strand(self, arrow):
        sign = self.sign()
        if arrow not in self:
            return None
        elif ( (arrow == self.over and sign == 'RH') or
               (arrow == self.under and sign =='LH') ):
            return 'X'
        else:
            return 'Y'

    def reverse(self):
        self.over, self.under = self.under, self.over

    def height(self, arrow):
        if arrow == self.under:
            return self.under ^ self.over
        elif arrow == self.over:
            return self.over ^ self.under
        else:
            return None

    def hit(self, count):
        if self.hit1 is None:
            self.hit1 = count
        elif self.hit2 is None:
            self.hit2 = count
        else:
            raise ValueError('Too many hits!')

    def clear_hits(self):
        self.hit1, self.hit2 = None, None

    def comp(self, component):
        if self.comp1 is None:
            self.comp1 = component
        elif self.comp2 is None:
            self.comp2 = component
        else:
            raise ValueError('Too many component hits!')

    def clear_comps(self):
        self.comp1, self.comp2 = None, None


class ECrossing:
    """
    A pair: (Crossing, Arrow), where the Arrow is involved in the Crossing.
    The ECrossings correspond 1-1 with edges of the link diagram.
    """ 
    def __init__(self, crossing, arrow):
        if arrow not in crossing:
            raise ValueError
        self.crossing = crossing
        self.arrow = arrow
        self.strand = self.crossing.strand(self.arrow)

    def pair(self):
        return (self.crossing, self.arrow)

    def goes_over(self):
        if self.arrow == self.crossing.over:
            return True
        return False


class LinkProjection:
    """
    A SnapPea link projection, instantiated from a link file.
    Can construct either a pythonified SnapPea KLPProjection
    or a DT code.
    """
    def __init__(self, filename):
        self.Arrows = []
        self.Vertices = []
        self.Crossings = []
        self.CrossPoints = []

        loadfile = open(filename, 'r')
        lines = [line for line in loadfile.readlines() if len(line) > 1]
        num_lines = len(lines)
        if not lines.pop(0).startswith('% Link Projection'):
            raise ValueError('This is not a SnapPea link projection file')
        try:
            num_components = int(lines.pop(0))
            for n in range(num_components):
                lines.pop(0) # We don't need this
            num_vertices = int(lines.pop(0))
            for n in range(num_vertices):
                x, y = lines.pop(0).split()
                X, Y = int(x), int(y)
                self.Vertices.append(Vertex(X, Y))
            num_arrows = int(lines.pop(0))
            for n in range(num_arrows):
                s, e = lines.pop(0).split()
                S, E = self.Vertices[int(s)], self.Vertices[int(e)]
                self.Arrows.append(Arrow(S, E))
            num_crossings = int(lines.pop(0))
            for n in range(num_crossings):
                u, o = lines.pop(0).split()
                U, O = self.Arrows[int(u)], self.Arrows[int(o)]
                self.Crossings.append(Crossing(O, U))
        except:
            raise ValueError('Failed while parsing line %d'%(
                num_lines - len(lines)))
        loadfile.close()

    def crossing_components(self):
        """
        Returns a list of lists of ECrossings, one per component,
        where the corresponding crossings are ordered consecutively
        through the component.  Requires that all components be closed.
        """
        for vertex in self.Vertices:
            if vertex.is_endpoint():
                raise ValueError('Found an endpoint!')
        result = []
        arrow_components = self.arrow_components()
        for component in arrow_components:
            crosses=[]
            for arrow in component:
                arrow_crosses = [(c.height(arrow), c, arrow) 
                                for c in self.Crossings if arrow in c]
                arrow_crosses.sort()
                crosses += arrow_crosses
            result.append([ECrossing(c[1],c[2]) for c in crosses]) 
        return result

    def arrow_components(self, include_isolated_vertices=False):
        """
        Returns a list of lists of arrows, one per component of the diagram.
        """
        pool = [v.out_arrow  for v in self.Vertices if not v.is_endpoint()]
        pool += [v.out_arrow for v in self.Vertices if v.in_arrow == None]
        result = []
        while len(pool):
            first_arrow = pool.pop()
            if first_arrow == None:
                continue
            component = [first_arrow]
            while component[-1].end != component[0].start:
                next = component[-1].end.out_arrow
                if next == None:
                    break
                pool.remove(next)
                component.append(next)
            result.append(component)
        if include_isolated_vertices:
            for vertex in [v for v in self.Vertices if v.is_isolated()]:
                result.append([Arrow(vertex, vertex, color=vertex.color)])
        return result

    def SnapPea_KLPProjection(self):
        """
        Constructs a python simulation of a SnapPea KLPProjection
        (Kernel Link Projection) structure.  See Jeff Weeks' SnapPea
        file link_projection.h for definitions.  Here the KLPCrossings
        are modeled by dictionaries.  This method requires that all
        components be closed.  A side effect is that the KLP attributes
        of all crossings are updated.

        The following excerpt from link_projection.h describes the
        main convention:

        If you view a crossing (from above) so that the strands go in the
        direction of the postive x- and y-axes, then the strand going in
        the x-direction is the KLPStrandX, and the strand going in the
        y-direction is the KLPStrandY.  Note that this definition does not
        depend on which is the overstrand and which is the understrand::
        
                                   KLPStrandY
                                       ^
                                       |
                                   ----+---> KLPStrandX
                                       |
                                       |
        """
        try:
            crossing_components = self.crossing_components()
        except ValueError:
            return None
        num_crossings = len(self.Crossings)
        num_free_loops = 0
        num_components = len(crossing_components)
        id = lambda x: self.Crossings.index(x.crossing)
        for component in crossing_components:
            this_component = crossing_components.index(component)
            N = len(component)
            for n in range(N):
                this = component[n]
                previous = component[n-1]
                next = component[(n+1)%N]
                this.crossing.KLP['sign'] = sign = this.crossing.sign()
                if this.strand == 'X':
                    this.crossing.KLP['Xbackward_neighbor'] = id(previous)
                    this.crossing.KLP['Xbackward_strand'] = previous.strand
                    this.crossing.KLP['Xforward_neighbor']  = id(next)
                    this.crossing.KLP['Xforward_strand'] = next.strand
                    this.crossing.KLP['Xcomponent'] = this_component
                else:
                    this.crossing.KLP['Ybackward_neighbor'] = id(previous)
                    this.crossing.KLP['Ybackward_strand'] = previous.strand
                    this.crossing.KLP['Yforward_neighbor']  = id(next)
                    this.crossing.KLP['Yforward_strand'] = next.strand
                    this.crossing.KLP['Ycomponent'] = this_component
            if N == 0:
                num_free_loops += 1
        KLP_crossings = [crossing.KLP for crossing in self.Crossings]
        return num_crossings, num_free_loops, num_components, KLP_crossings

    def DT_code(self, alpha=True):
        """
        Return the Dowker-Thistlethwaite code as a list of tuples
        of even integers.

        If alpha is set to True, return the alphabetical
        Dowker-Thistlethwaite code as used by Oliver Goodman's Snap.
        """
        components = self.crossing_components()
        lookup = {}
        # sort the components by increasing length
        components.sort(key=lambda x: len(x))
        for crossing in self.Crossings:
            crossing.clear_hits()
            crossing.clear_comps()
        # mark which component each crossing belongs to
        for component in components:
            for ecrossing in component:
                ecrossing.crossing.comp(component)
        count = 1
        chunks = []
        prefix_ints = [len(self.Crossings), len(components)]
        while len(components) > 0:
            this_component = components.pop(0)
            # Choose the first crossing, by Morwen's rules:
            # If any crossing on this component has been hit
            # find the first one with an odd label and then
            # start at its predecessor.
            odd_hits = [ec for ec in this_component if
                        ec.crossing.hit1 is not None and
                        ec.crossing.hit1%2 == 1]
            if len(odd_hits) > 0:
                odd_hits.sort(key=lambda x : x.crossing.hit1)
                n = this_component.index(odd_hits[0])
                this_component = this_component[n-1:] + this_component[:n-1]
            # Label the crossings on this component
            odd_count = 0
            these_crossings = set()
            for ecrossing in this_component:
                crossing = ecrossing.crossing
                these_crossings.add(crossing)
                if count%2 == 0 and ecrossing.goes_over():
                    crossing.hit(-count)
                else:
                    crossing.hit(count)
                if count%2 == 1:
                    odd_count += 1
                count += 1
            chunks.append(odd_count)
            # Check for crossings shared with unfinished components.
            odd_shared = [c for c in self.Crossings if
                          c.hit1 != None and
                          c.hit1%2 == 1 and 
                          c.comp1 != c.comp2 and 
                          c.comp2 in components]
            if len(odd_shared) > 0:
                # Choose the next component, by Morwen's rules:
                # Use the component containing the partner of the
                # first odd-numbered crossing that is shared with
                # another commponent (if there are any).
                odd_shared.sort(key=lambda x : x.hit1)
                next_component = odd_shared[0].comp2
                components.remove(next_component)
                components.insert(0, next_component)
        # build the Dowker-Thistlethwaite code
        even_codes = [None]*len(self.Crossings)
        for crossing in self.Crossings:
            if crossing.hit1%2 != 0:
                even_codes[(crossing.hit1 - 1)//2] = crossing.hit2
            else:
                even_codes[(crossing.hit2 - 1)//2] = crossing.hit1
        if not alpha:
            result = []
            for chunk in chunks:
                result.append(tuple(even_codes[:chunk]))
                even_codes = even_codes[chunk:]
            return result
        else:
            alphacode = ''.join(tuple([DT_alphabet[x>>1] for x in even_codes]))
            prefix_ints += chunks
            if prefix_ints[0] > 26:
                raise ValueError('Too many crossings!')
            prefix = ''.join(tuple([DT_alphabet[n] for n in prefix_ints]))
            return prefix + alphacode
