### Required modules:
# Some standard modules.
from __future__ import print_function
import os
from random import choice
from itertools import combinations
# Some custom modules.
import snappy
from plink import LinkManager
from .twister_core import build_bundle, build_splitting, twister_version

# Python 3 compatibility
try:
	basestring
except NameError:  # Python 3
	basestring = unicode = str

surface_database_path = os.path.join(os.path.dirname(__file__), 'surfaces')
surface_database = set(os.listdir(surface_database_path))
version = twister_version()

def _get_surface(surface):
	if isinstance(surface, tuple) and len(surface) == 2 and isinstance(surface[0], int) and isinstance(surface[1], int):
		return LP_surface(surface[0], surface[1])
	
	if isinstance(surface, basestring):
		# If surface is actually the contents of a surface file.
		if surface.startswith('# A Twister surface file'):
			return surface
		
		# If surface is in the database.
		if surface in surface_database:
			surface = os.path.join(surface_database_path, surface)
		
		try:
			lines = open(surface, 'r').readlines() 
		except IOError:
			raise IOError('Unable to open %s' % surface)
		
		contents = ''.join(lines)
		if contents.startswith('% Virtual Link Projection\n'):
			# We use PLink to handle virtual link projections.
			LM = LinkManager()
			LM._from_string(contents)
			return LM.Twister_surface_file()
		if contents.startswith('# A Twister surface file'):
			return contents
		
		raise TypeError('Not a Twister surface file.')
	
	raise TypeError('Surfaces can only be loaded from a string or pair of integers (genus,boundary).')

def _parse_line(line):
	data = line.split('#')
	if data[0] == '': return None
	return data[0].split(',')

def _parse_surface(surface_contents):
	# We now compute some properties of the surface. These are useful for error checking.
	# We could do more error checking to make sure that the surface is valid. But this
	# would just be repeating code in twister_core.
	curves = {'annulus':[], 'rectangle':[], 'macro':[]}
	num_squares = 0;
	for line in surface_contents.split('\n'):
		data = _parse_line(line)
		if data is not None:
			if len(data) > 2:
				curves[data[0]].append((data[1], data[2], data[3:]))
				if data[3:] != []:
					largest_num = max(map(lambda x: abs(int(x)) + 1, data[3:]))
					if largest_num > num_squares:
						num_squares = largest_num
	
	vertices = [set([str(i) + d]) for i in range(num_squares) for d in ['NE', 'NW', 'SE', 'SW']]
	
	def join(v, a, b):
		X = [s for s in v if a in s or b in s]
		if len(X) == 2:
			v.remove(X[0])
			v.remove(X[1])
			v.append(X[0] | X[1])
	
	for annulus in curves['annulus']:
		data = annulus[2]
		for i in range(len(data)):
			j = (i+1) % len(data)
			sign1 = data[i][0] == '-'
			sign2 = data[j][0] == '-'
			square1 = str(abs(int(data[i])))
			square2 = str(abs(int(data[j])))
			corner1a = 'NW' if data[i][0] == '-' else 'NE'
			corner2a = 'SW' if data[j][0] == '-' else 'NW'
			corner1b = 'NE' if data[i][0] == '-' else 'SE'
			corner2b = 'SE' if data[j][0] == '-' else 'SW'
			
			join(vertices, square1 + corner1a, square2 + corner2a)
			join(vertices, square1 + corner1b, square2 + corner2b)
	
	for rectangle in curves['rectangle']:
		data = rectangle[2]
		for i in range(len(data)-1):
			j = (i+1) % len(data)
			square1 = str(abs(int(data[i])))
			square2 = str(abs(int(data[j])))
			corner1a = 'NW' if data[i][0] == '-' else 'NE'
			corner2a = 'SW' if data[j][0] == '-' else 'NW'
			corner1b = 'NE' if data[i][0] == '-' else 'SE'
			corner2b = 'SE' if data[j][0] == '-' else 'SW'
			
			join(vertices, square1 + corner1a, square2 + corner2a)
			join(vertices, square1 + corner1b, square2 + corner2b)
	
	end_edges = []
	for rectangle in curves['rectangle']:
		data = rectangle[2]
		
		square1 = str(abs(int(data[0])))
		corner1a = 'SW' if data[0][0] == '-' else 'SW'
		corner1b = 'SE' if data[0][0] == '-' else 'NW'
		end_edges.append((square1 + corner1a, square1 + corner1b))
		
		square2 = str(abs(int(data[-1])))
		corner2a = 'NW' if data[0][0] == '-' else 'NE'
		corner2b = 'NE' if data[0][0] == '-' else 'SE'
		end_edges.append((square2 + corner2a, square2 + corner2b))
	
	boundary_components = [set([e]) for e in end_edges]
	for e1, e2 in combinations(end_edges, r=2):
		if any(v1 in v and v2 in v for v in vertices for v1 in e1 for v2 in e2):
			join(boundary_components, e1, e2)
	
	return num_squares, curves, vertices, boundary_components

def _compute_intersections(curves):
	curve_names = [name for (name, inverse_name, data) in curves['annulus']] + [name for (name, inverse_name, data) in curves['rectangle']]
	data = dict((name, set([abs(int(x)) for x in data])) for (name, inverse_name, data) in curves['annulus'] + curves['rectangle'])
	
	rows = [[''] + curve_names]
	for name in curve_names:
		row = [name]
		for name2 in curve_names:
			row.append(len(data[name].intersection(data[name2])) if name != name2 else 0)
		rows.append(row)
	
	return rows

class Surface:
	''' Represents a squared surface with a collection of curves on it. '''
	def __init__(self, surface):
		''' Creates a new surface from either:
		- the name of a surface in the surface database,
		- a path to a surface file,
		- a path to a plink virtual knot file,
		- the contents of a surface file, or
		- the pair (genus, boundary). For information about how this surface
		\t is created see the doc string of twister.LP_surface(). '''
		self.surface_contents = _get_surface(surface)
		self.num_squares, self.curves, vertices, self.boundary_components = _parse_surface(self.surface_contents)
		
		self.num_boundary = len(self.boundary_components)
		self.num_edges = sum(len(x[2]) for x in self.curves['annulus']) + sum(len(x[2]) + 1 for x in self.curves['rectangle'])
		self.num_vertices = len(vertices)
		self.Euler_characteristic = self.num_vertices - self.num_edges + self.num_squares
		self.genus = (2 - self.num_boundary - self.Euler_characteristic) // 2
		self.intersection_matrix = _compute_intersections(self.curves)
	
	def info(self, verbose=False):
		annuli = ', '.join(a[0] for a in self.curves['annulus'])
		inv_annuli = ', '.join(a[1] for a in self.curves['annulus'])
		rectangles = ', '.join(r[0] for r in self.curves['rectangle'])
		inv_rectangles = ', '.join(r[1] for r in self.curves['rectangle'])
		macros = ', '.join('%s:%s' % (r[0], r[1]) for r in self.curves['macro'])
		matrix = '\n'.join(str(row) for row in self.intersection_matrix)
		print('A Twister surface of genus %d with %d boundary component(s)' % (self.genus, self.num_boundary))
		print('Loops: %s' % annuli)
		if verbose: print('Inverse names: %s' % inv_annuli)
		print('Arcs: %s' % rectangles)
		if verbose: print('Inverse names: %s' % inv_rectangles)
		if verbose: print('Macros: %s' % macros)
		if verbose: print('Intersection matrix:\n%s' % matrix)
	
	def random_word(self, n, twists=True, half_twists=True, macros=True, inverses=True):
		''' Returns a random word of length n in the generators of Mod(S). Setting twists, half_twists or macros to False
		prevents them from appearing in the word. If all generators are disallowed then the empty word is returned. If 
		inverses is set to False then no inverse of a generator will be used. '''
		
		generators = []
		if twists: generators += [curve[0] for curve in self.curves['annulus']]
		if twists and inverses: generators += [curve[1] for curve in self.curves['annulus']]
		if half_twists: generators += [curve[0] for curve in self.curves['rectangle']]
		if half_twists and inverses: generators += [curve[1] for curve in self.curves['rectangle']]
		if macros: generators += [curve[0] for curve in self.curves['macro']]
		
		if generators == []: return ''
		
		return ''.join([choice(generators) for i in range(n)])
	
	def bundle(self, monodromy, name=None, optimize=True, warnings=True, debugging_level=0, return_type='manifold'):
		''' Generate a surface bundle over a circle with fibre this surface from mapping class group data using the Twister
		program of Bell, Hall and Schleimer.
		
		Arguments:
		Required:
		 monodromy - build a surface bundle with specified monodromy
		Optional:
		 name - name of the resulting manifold
		 optimize - try to reduce the number of tetrahedra (default True)
		 warnings - print Twister warnings (default True)
		 debugging_level - specifies the amount of debugging information to be shown (default 0)
		 return_type - specifies how to return the manifold, either as a 'manifold' (default), 'triangulation' or 'string'
		
		Monodromy is a word of annulus and rectangle names (or
		their inverses).  These are read from left to right and determine a
		sequence of (half) Dehn twists.  When prefixed with an "!" the name
		specifies a drilling.  For example, "a*B*a*B*A*A*!a*!b" will perform 6
		twists and then drill twice.
		
		Examples:
		
		The figure eight knot complement:
		>>> M = twister.Surface((1,1)).bundle(monodromy='a_0*B1')
		
		The minimally twisted six chain link:
		>>> M = twister.Surface('S_1_1').bundles(monodromy='!a*!b*!a*!b*!a*!b')
		>>> M.set_peripheral_curves('shortest_meridians', 0)
		>>> M.dehn_fill((1,0),0) '''
		
		if name is None: name = monodromy
		tri, messages = build_bundle(name, self.surface_contents, monodromy, optimize, True, warnings, debugging_level)
		
		# You might want to change what is done with any warning / error messages.
		# Perhaps they should be returned in the next block?
		if messages != '': print(messages)
		
		if tri is None: 
			return None
		
		return_type = return_type.lower()
		if return_type == 'manifold':
			return snappy.Manifold(tri)
		if return_type == 'triangulation':
			return snappy.Triangulation(tri)
		if return_type == 'string':
			return tri
		
		raise TypeError('Return type must be \'manifold\', \'triangulation\' or \'string\'.')
	
	def splitting(self, gluing, handles, name=None, optimize=True, warnings=True, debugging_level=0, return_type='manifold'):
		''' Generate a manifold with Heegaard splitting this surface from mapping class group data using the Twister
		program of Bell, Hall and Schleimer.
		
		Arguments:
		Required:
		 gluing - the gluing used to join the upper and lower compression bodies
		 handles - where to attach 2-handles
		Optional:
		 name - name of the resulting manifold
		 optimize - try to reduce the number of tetrahedra (default True)
		 warnings - print Twister warnings (default True)
		 debugging_level - specifies the amount of debugging information to be shown (default 0)
		 return_type - specifies how to return the manifold, either as a 'manifold' (default), 'triangulation' or 'string'
		
		Gluing is a word of annulus and rectangle names (or
		their inverses).  These are read from left to right and determine a
		sequence of (half) Dehn twists.  When prefixed with an "!" the name
		specifies a drilling.  For example, "a*B*a*B*A*A*!a*!b" will perform 6
		twists and then drill twice.
		
		Handles is again a word of annulus names (or inverses).  For example,
		'a*c*A' means attach three 2-handles, two above and one below.
		
		Examples:
		
		The genus two splitting of the solid torus:
		>>> M = twister.Surface('S_2').splitting(gluing='', handles='a*B*c') '''
		
		if name is None: name = gluing + ' ' + handles
		tri, messages = build_splitting(name, self.surface_contents, gluing, handles, optimize, True, warnings, debugging_level)
		
		# You might want to change what is done with any warning / error messages.
		# Perhaps they should be returned in the next block?
		if messages != '': print(messages)
		
		if tri is None: 
			return None
		
		return_type = return_type.lower()
		if return_type == 'manifold':
			return snappy.Manifold(tri)
		if return_type == 'triangulation':
			return snappy.Triangulation(tri)
		if return_type == 'string':
			return tri
		
		raise TypeError('Return type must be \'manifold\', \'triangulation\' or \'string\'.')

### Some standard types of surface that can be automatically generated.

def LP_surface(genus, boundary, make_prefix_unique=True):
	''' Returns the contents of a surface file for the surface S_{genus, boundary}.
	We generally follow the naming convention given in Figure 13 
	of the Labruere and Paris paper "Presentations for the punctured 
	mapping class groups in terms of Artin groups".
	
	When genus == 1, the loop a_n is dropped as it is isotopic to the
	loop a_0.
	
	The case of genus == 0, given here, is a small modification of
	Figure 11 in [LP].
	
	You should not use this function directly, rather just call
	Surface((genus, boundary)). '''
	
	contents = ['# A Twister surface file produced by LP_surface.']
	contents.append('# with generating set for MCG(S_{%d,%d}) following Figure 13 of Labruere and' % (genus, boundary))
	contents.append('# Paris paper "Presentations for the punctured mapping class groups')
	contents.append('# in terms of Artin groups".')
	contents.append('#')
	
	padded_length = len(str(max(2*genus, boundary))) if make_prefix_unique else 0
	
	square_count = 0
	
	if genus == 0:
		if boundary == 0:
			square_count += 1
			contents.append('annulus,a1,A1,-0,+1#')
			contents.append('annulus,a2,A2,-1,+0#')
		else:
			for i in range(boundary):
				square_count += 1
				start = 'rectangle,t' + str(i + 1).zfill(padded_length) + ',T' + str(i + 1).zfill(padded_length)
				connections = ',-' + str(square_count - 1) + ',+' + str(square_count)
				contents.append(start + connections + '#')
				
				if i == boundary - 1:
					start = 'annulus,a' + str(i + 1).zfill(padded_length) + ',A' + str(i + 1).zfill(padded_length)
					connections = ',-' + str(square_count) + ',+' + '0'
					contents.append(start + connections + '#')
				else:
					square_count += 1
					start = 'annulus,a' + str(i + 1).zfill(padded_length) + ',A' + str(i + 1).zfill(padded_length)
					connections = ',-' + str(square_count - 1) + ',+' + str(square_count)
					contents.append(start + connections + '#')
	else:
		c_loop = -1  # Just a marker.
		
		# Start with the a0 loop.
		next_line = 'annulus,a' + '0'.zfill(padded_length) + ',A' + '0'.zfill(padded_length) + ',+0'
		contents.append(next_line + '#')
		
		# And the b1 loop.
		start = 'annulus,b' + '1'.zfill(padded_length) + ',B' + '1'.zfill(padded_length)
		connections = ',-0'
		
		if boundary > 1:
			for i in range(boundary - 1):
				square_count += 1
				connections = connections + ',-' + str(square_count)
				
				square_count += 1
				
				start2 = 'annulus,a' + str(i + 1).zfill(padded_length) + ',A' + str(i + 1).zfill(padded_length)
				connections2 = ',+' + str(square_count - 1) + ',+' + str(square_count)
				contents.append(start2 + connections2 + '#')
				
				start3 = 'rectangle,t' + str(i + 1).zfill(padded_length) + ',T' + str(i + 1).zfill(padded_length)
				connections3 = ',-' + str(square_count)
				contents.append(start3 + connections3 + '#')
		
		elif boundary == 1:
			square_count += 1
			connections = connections + ',-' + str(square_count)
			
			start2 = 'rectangle,t' + '1'.zfill(padded_length) + ',T' + '1'.zfill(padded_length)
			connections2 = ',+' + str(square_count)
			contents.append(start2 + connections2 + '#')
		
		elif boundary == 0 and genus > 1:
			square_count += 1
			connections = connections + ',+' + str(square_count)
		
		# Add in an extra a loop (if needed) to isolate the half-twists from the rest of the surface.
		if genus > 1 and boundary > 0:
			square_count += 1
			connections = connections + ',-' + str(square_count)
			
			start2 = 'annulus,a' + str(boundary).zfill(padded_length) + ',A' + str(boundary).zfill(padded_length)
			connections2 = ',+' + str(square_count)
			contents.append(start2 + connections2 + '#')
			
			# Add in the start point for the b2 ... chain.
			square_count += 1
			connections = connections + ',+' + str(square_count)
		
		contents.append(start + connections + '#')
		
		# Now construct the rest of the b loops.
		for i in range(2, 2 * genus - 1):
			square_count += 1
			start = 'annulus,b' + str(i).zfill(padded_length) + ',B' + str(i).zfill(padded_length)
			connections = ',-' + str(square_count - 1) + ',+' + str(square_count)
			contents.append(start + connections + '#')
			if i == 2: c_loop = len(contents)
		
		# Don't forget, the last one is different.
		if genus > 1:
			square_count += 1
			start = 'annulus,b' + str(2 * genus - 1).zfill(padded_length) + ',B' + str(2 * genus - 1).zfill(padded_length)
			connections = ',-' + str(square_count - 1)
			contents.append(start + connections + '#')
		
		# Modify the 3rd b loop to add in the c loop (if needed).
		if c_loop > -1:
			next_line = 'annulus,c,C,-' + str(square_count)
			contents.append(next_line + '#')
			contents[c_loop] = contents[c_loop][:-1] + ',+' + str(square_count) + '#'
	
	return '\n'.join(contents)

def code_to_sign_sequence(code):
	''' Produces a sign sequence for a given Dowker code. '''
	
	def first_non_zero(L): return min(i for i in range(len(L)) if L[i])
	
	N = len(code)
	signs = [abs(n) - 1 for n in code]
	
	pairs = list(zip(range(0, 2*N, 2), signs))
	pairs_dict = dict([(x, y) for x, y in pairs] + [(y, x) for x, y in pairs])
	full_code = [pairs_dict[i] for i in range(2*N)]
	
	seq = full_code * 2  # seq is two copies of full DT involution on crossings numbered 0 to 2N-1.
	emb, A = [0] * 2*N, [0] * 2*N  # zero emb and A. A will only ever contain zeroes and ones.
	
	# Set initial conditions.
	A[0], A[seq[0]] = 1, 1
	emb[0], emb[seq[0]] = 1, -1
	
	# Determine the possible phi's
	all_phi = [[0] * 2*N for i in range(2*N)]
	for i in range(2*N):
		all_phi[i][i] = 1
		for j in range(i, i+2*N):
			all_phi[i][j % (2*N)] = 1 if i == j else -all_phi[i][(j-1) % (2*N)] if i <= seq[j] <= seq[i] else all_phi[i][(j-1) % (2*N)]
	
	while any(A):
		i = first_non_zero(A)  # let i be the index of the first non-zero member of A
		
		D = [1] * 2*N
		D[i:seq[i]+1] = [0] * (seq[i] - i + 1)
		while any(D):
			x = first_non_zero(D)  # let x be the index of the first non-zero member of D
			D[x] = 0
			
			if ((i <= seq[x] <= seq[i] and emb[x] != 0 and all_phi[i][x] * all_phi[i][seq[x]] * emb[i] != emb[x]) or ((seq[x] < i or seq[i] < seq[x]) and all_phi[i][x] * all_phi[i][seq[x]] != 1)) and x < i:  # This extra AND conditions shouldn't be needed.
				raise ValueError('Not a realisable DT-code.')
			
			if seq[i] < seq[x] or seq[x] < i:
				D[seq[x]] = 0
			elif emb[x] == 0: # emb[x] is already defined
				assert D[seq[x]] == 0
				emb[x] = all_phi[i][x] * all_phi[i][seq[x]] * emb[i]
				emb[seq[x]] = -emb[x]
				if abs(seq[x]-seq[x-1]) % (2*N) != 1:
					A[x] = 1
					A[seq[x]] = 1
		
		A[i], A[seq[i]] = 0, 0
	
	return [emb[i] for i in range(0, 2*N, 2)]  # Note [emb[pairs_dict[i]] for i in range(0, 2*N, 2)] is also a valid code.

def DT_drilling_surface(code, make_prefix_unique=True):
	''' Returns a Surface which can be used to construct a triangulation of the
	complement of the knot with Dowker--Thistlethwaite code 'code'. 
	
	Example: DT_drilling_surface([6,10,8,12,4,2]) returns a Surface that can be
	used to construct the 6_1 knot complement. '''
	
	# Note: There are smaller surface files that can be constructed to do this,
	# but they are significantly more complicated and have less symmetry.
	code = list(code)
	signs = code_to_sign_sequence(code)
	
	contents = ['# A Twister surface file']
	contents.append('#')
	contents.append('# Surface file for knot with Dowker code:')
	contents.append('#    ' + ','.join(map(str, code)))
	contents.append('#')
	contents.append('# To build this knot complement, make a Heegaard splitting with a 2-handle attached above')
	contents.append('# and below every annulus and drill every rectangle exactly once, making sure to drill all')
	contents.append('# of the "x" rectangles before drilling ANY of the "y" rectangles.')
	contents.append('#')
	
	num_crossings = len(code)
	padded_length = len(str(2*num_crossings)) if make_prefix_unique else 0
	num_squares = num_crossings * 5
	
	pairs = list(zip(range(1, 2*num_crossings+1, 2), code))
	pairs_dict = dict([(x, abs(y)) for (x, y) in pairs] + [(abs(y), x) for (x, y) in pairs])
	signs_dict = dict([(x, y > 0) for (x, y) in pairs] + [(abs(y), y > 0) for (x, y) in pairs])
	
	drill_names = []
	handle_names = []
	
	# First build all the rectangles.
	for i in range(1, 2*num_crossings+1):
		if i % 2 == 1:
			k = i // 2
			cells = ['-' + str(k * 5 + 1), '+' + str(k * 5 + 2), '+' + str(k * 5 + 3)]
		else:
			k = pairs_dict[i] // 2
			if signs[k] == 1:
				cells = ['-' + str(k * 5 + 4), '-' + str(k * 5 + 2), '+' + str(k * 5 + 0)]
			else:
				cells = ['-' + str(k * 5 + 0), '-' + str(k * 5 + 2), '+' + str(k * 5 + 4)]
		
		over = (i % 2 == 0) ^ signs_dict[i]
		name = ('y' if (i % 2 == 0) ^ signs_dict[i] else 'x') + str(i).zfill(padded_length)
		inverse_name = name.swapcase()
		contents.append(','.join(['rectangle', name, inverse_name] + cells) + '#')
		drill_names.append('!' + name)
	
	# Then the required annuli.
	for i in range(1, 2*num_crossings+1):
		j = i % (2*num_crossings) + 1
		if i % 2 == 1:
			k = i // 2
			l = pairs_dict[j] // 2
			if signs[l] == 1:
				cells = ['-' + str(k * 5 + 3), '+' + str(l * 5 + 4)]
			else:
				cells = ['-' + str(k * 5 + 3), '-' + str(l * 5 + 4)]
		else:
			k = pairs_dict[i] // 2
			l = j // 2
			if signs[k] == 1:
				cells = ['-' + str(k * 5 + 0), '+' + str(l * 5 + 1)]
			else:
				cells = ['+' + str(k * 5 + 0), '+' + str(l * 5 + 1)]
		
		name = 'a' + '_' + str(i).zfill(padded_length)
		inverse_name = name.swapcase()
		contents.append(','.join(['annulus', name, inverse_name] + cells) + '#')
		handle_names.append(name)
		handle_names.append(inverse_name)
	
	contents.append('# We also give a macro to drill all of the rectangles in an acceptable order.')
	contents.append('# And one for attaching all handles too.')
	drill_names.sort()  # Make sure to sort the rectangle drill names so that all of the "x"'s appear before ANY of the "y"'s.
	contents.append('macro,s,' + '*'.join(drill_names) + '#')
	contents.append('macro,h,' + '*'.join(handle_names) + '#')
	
	return Surface('\n'.join(contents))

def DT_handles_surface(code, make_prefix_unique=True):
	''' Returns a Surface which can be used to construct a triangulation of the
	complement of the knot with Dowker code 'code'. 
	
	Example: DT_handles_surface([6,10,8,12,4,2]) returns a Surface that can be
	used to construct the 6_1 knot complement. '''
	
	# Note: There are smaller surface files that can be constructed to do this,
	# but they are significantly more complicated and have less symmetry.
	code = list(code)
	signs = code_to_sign_sequence(code)  # Get the orientation of each intersection
	
	contents = ['# A Twister surface file']
	contents.append('#')
	contents.append('# Surface file for knot with Dowker code:')
	contents.append('#    ' + ','.join(map(str, code)))
	contents.append('#')
	contents.append('# To build this knot complement, make a Heegaard splitting with a 2-handle attached above')
	contents.append('# and below the sequence of annuli specified in the macro "h".')
	contents.append('#')
	
	num_crossings = len(code)
	
	pairs = list(zip(range(1, 2*num_crossings+1, 2), code))
	pairs_dict = dict([(x, abs(y)) for (x, y) in pairs] + [(abs(y), x) for (x, y) in pairs])
	signs_dict = dict([(x, y > 0) for (x, y) in pairs] + [(abs(y), y > 0) for (x, y) in pairs])
	
	overcrossing = [signs_dict[i+1] ^ (i % 2 == 1) for i in range(2*num_crossings)]
	where_switch = list(filter(lambda i: not signs_dict[i+1] ^ signs_dict[((i+1) % (2*num_crossings)) + 1], range(2*num_crossings)))
	last_true = max(where_switch) + 1
	
	num_annuli = len(where_switch)
	squares_in_crossings = 4 * num_crossings
	squares_in_links = 2 * num_annuli
	squares_in_rectangles = 2 * num_annuli
	
	num_squares = squares_in_crossings + squares_in_links + squares_in_rectangles
	padded_length = len(str(num_annuli)) if make_prefix_unique else 0
	
	annuli_count = 0
	squares = []
	back_squares = []
	handle_names = []
	for j in range(2*num_crossings):
		i = ((j + last_true) % (2*num_crossings))
		if i % 2 == 0:
			k = i // 2
			squares.append('-' + str(4*k+0))
			squares.append('+' + str(4*k+1))
			back_squares.append('+' + str(4*k+3))
			back_squares.append('-' + str(4*k+2))
		else:
			k = pairs_dict[i+1] // 2
			if signs[k] == +1:
				squares.append('-' + str(4*k+3))
				squares.append('+' + str(4*k+0))
				back_squares.append('+' + str(4*k+2))
				back_squares.append('-' + str(4*k+1))
			else:
				squares.append('-' + str(4*k+1))
				squares.append('+' + str(4*k+2))
				back_squares.append('+' + str(4*k+0))
				back_squares.append('-' + str(4*k+3))
		
		if i in where_switch:
			row = 'annulus,a' + str(annuli_count).zfill(padded_length) + ',A' + str(annuli_count).zfill(padded_length) + ','
			row += '+' + str(squares_in_crossings + 2*annuli_count) + ','
			row += ','.join(squares) + ','
			row += '-' + str(squares_in_crossings + 2*((annuli_count+1) % num_annuli)) + ','
			row += '-' + str(squares_in_crossings + squares_in_links + 2*annuli_count) + ','
			row += '+' + str(squares_in_crossings + squares_in_links + 2*annuli_count + 1) + ','
			row += '+' + str(squares_in_crossings + 2*((annuli_count+1) % num_annuli) + 1) + ','
			row += ','.join(back_squares[::-1]) + ','
			row += '-' + str(squares_in_crossings + 2*annuli_count + 1)
			row += '#'
			
			contents.append(row)
			
			handle_names.append(('a' if overcrossing[i] else 'A') + str(annuli_count).zfill(padded_length))
			
			annuli_count += 1
			squares = []
			back_squares = []
	
	for i in range(num_annuli):
		row = 'rectangle,t' + str(i).zfill(padded_length) + ',T' + str(i).zfill(padded_length) + ','
		row += '+' + str(squares_in_crossings + squares_in_links + 2*i) + ','
		row += '-' + str(squares_in_crossings + squares_in_links + 2*i + 1)
		row += '#'
		contents.append(row)
		
	
	num_squares = num_crossings * 5
	
	contents.append('# We also give a macro for attaching the required handles.')
	contents.append('macro,h,' + '*'.join(handle_names) + '#')
	
	return Surface('\n'.join(contents))
