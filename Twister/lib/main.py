### Required modules:
# Some standard modules.
from __future__ import print_function
import os
# Some custom modules.
import snappy
from .twister_core import build_bundle, build_splitting
# Python 3 compatibility
try:
	basestring
except NameError: # Python 3
	basestring = unicode = str

surface_database_path = os.path.join(os.path.dirname(__file__), 'surfaces')
surface_filter = lambda path: os.path.splitext(path)[1] == '.sur'
surface_database = set(filter(surface_filter, os.listdir(surface_database_path)))

def get_surface(file):
	''' To load the contents of a surface file either from the surface database
	or from a path to a file. '''
	
	if file in surface_database:
		file = os.path.join(surface_database_path, file)
	
	try:
		return ''.join(open(file, 'r'))
	except IOError:
		raise IOError('Unable to open %s' % file)

def twister(surface=None, surface_contents=None, monodromy=None, gluing=None, handles=None, name=None, 
	optimize=True, peripheral_curves=True, warnings=True, debugging_level=0, with_hyperbolic_structure=True):
	
	'''
	Generate a manifold from mapping class group data using the Twister
	program of Bell, Hall and Schleimer.
	
	Arguments:
	 surface - the contents of a Twister surface file, or the pair (genus, punctures)
	 monodromy - build a surface bundle with specified monodromy
	 gluing, handles - build a Heegaard splitting with specified gluing and handles
	 name - name of the resulting manifold
	 with_hyperbolic_structure - return a Manifold (if True) or a Triangulation
	 peripheral_curves - install canonical peripheral curves (default True)
	 optimize - try to reduce the number of tetrahedra (default True)
	 warnings - print Twister warnings (default True)
	 debugging_level - sets the amount of debugging information to be shown (default 0)
	
	The surface argument must be provided, and the final manifold is built
	starting with (surface x I).
	
	"monodromy" and "gluing" are words of annulus and rectangle names (or
	their inverses).  These are read from left to right and determine a
	sequence of (half) Dehn twists.  When prefixed with an "!" the name
	specifies a drilling.  For example, "a*B*a*B*A*A*!a*!b" will perform 6
	twists and then drill twice.
	
	"handles" is again a word of annulus names (or inverses).  For example,
	'a*c*A' means attach three 2-handles, two above and one below.
	
	Examples:
	
	The figure eight knot complement:
	>>> M = twister(surface=(1,1), monodromy='a_0*B_1')
	
	The genus two splitting of the solid torus:
	>>> M = twister(surface='S_2.sur', handles='a*B*c')
	
	The minimally twisted six chain link:
	>>> M = twister(surface='S_1_1.sur','!a*!b*!a*!b*!a*!b')
	>>> M.set_peripheral_curves('shortest_meridians', 0)
	>>> M.dehn_fill((1,0),0)
	'''
	
	if surface is None and surface_contents is None:
		raise ValueError('No surface file or surface file contents specified.')
	
	if surface is not None and surface_contents is None:
		if isinstance(surface, basestring):
			surface_contents = get_surface(surface)
		else:
			if len(surface) != 2: 
				raise ValueError('Surface must either be a surface file or a pair (genus, punctures).')
			
			genus, punctures = surface
			surface_contents = LP_surface(genus, punctures)
	
	if monodromy is not None and (gluing is not None or handles is not None):
		raise ValueError('Please specify *either* a bundle *or* a Heegaard splitting.')
	if monodromy is None and gluing is None and handles is None:
		raise ValueError('Please specify at least one of {monodromy, gluing, handles}.')
	
	if monodromy is not None:
		if name is None: name = monodromy
		tri, messages = build_bundle(name, surface_contents, monodromy, optimize, peripheral_curves, warnings, debugging_level)
	else:
		if gluing is None: gluing = ''
		if handles is None: handles = ''
		if name is None: name = gluing + ' ' + handles
		tri, messages = build_splitting(name, surface_contents, gluing, handles, optimize, peripheral_curves, warnings, debugging_level)
	
	# You might want to change what is done with any warning / error messages.
	# Perhapse they should be returned in the next block?
	if messages != '': print(messages)
	
	if tri is None: 
		return None
	else:
		if with_hyperbolic_structure:
			return snappy.Manifold(tri)  
		else:
			return snappy.Triangulation(tri)

### Some standard types of surface that can be automatically generated.

def LP_surface(genus, punctures, make_prefix_unique=True):
	'''
	Returns a list, each entry of which is a line in a surface file
	for the surface S_{genus, punctures}. We generally follow the
	naming convention given in Figure 13 of the Labruere and Paris
	paper "Presentations for the punctured mapping class groups in
	terms of Artin groups".
	
	When genus == 1, the loop a_n is dropped as it is isotopic to the
	loop a_0.
	
	The case of genus == 0, given here, is a small modification of
	Figure 11 in [LP].
	
	Example: construct_file_contents(3,2) returns the standard surface
	file for S_{3,2} broken by line into a list.
	'''
	
	contents = ['#']
	contents.append('# Generating set for MCG(S_{' + str(genus) +',' + str(punctures) + '}) following Figure 13 of Labruere and')
	contents.append('# Paris paper "Presentations for the punctured mapping class groups')
	contents.append('# in terms of Artin groups".')
	contents.append('#')
	contents.append('0#')
	
	padded_length = len(str(max(2*genus, punctures))) if make_prefix_unique else 0
	
	square_count = 0
	
	if genus == 0:
		if punctures == 0:
			square_count += 1
			contents.append('annulus,a_1,A_1,-0,+1#')
			contents.append('annulus,a_2,A_2,-1,+0#')
		else:
			for i in range(punctures):
				square_count += 1
				start = 'rectangle,t_' + str(i + 1).zfill(padded_length) + ',T_' + str(i + 1).zfill(padded_length)
				connections = ',-' + str(square_count - 1) + ',+' + str(square_count)
				contents.append(start + connections + '#')
				
				if i == punctures - 1:
					start = 'annulus,a_' + str(i + 1).zfill(padded_length) + ',A_' + str(i + 1).zfill(padded_length)
					connections = ',-' + str(square_count) + ',+' + '0'
					contents.append(start + connections + '#')
				else:
					square_count += 1
					start = 'annulus,a_' + str(i + 1).zfill(padded_length) + ',A_' + str(i + 1).zfill(padded_length)
					connections = ',-' + str(square_count - 1) + ',+' + str(square_count)
					contents.append(start + connections + '#')
	else:
		c_loop = -1  # Just a marker.
		
		# Start with the a_0 loop.
		next_line = 'annulus,a_' + '0'.zfill(padded_length) + ',A_' + '0'.zfill(padded_length) + ',+0'
		contents.append(next_line + '#')
		
		# And the b_1 loop.
		start = 'annulus,b_' + '1'.zfill(padded_length) + ',B_' + '1'.zfill(padded_length)
		connections = ',-0'
		
		if punctures > 1:
			for i in range(punctures - 1):
				square_count += 1
				connections = connections + ',-' + str(square_count)
				
				square_count += 1
				
				start2 = 'annulus,a_' + str(i + 1).zfill(padded_length) + ',A_' + str(i + 1).zfill(padded_length)
				connections2 = ',+' + str(square_count - 1) + ',+' + str(square_count)
				contents.append(start2 + connections2 + '#')
				
				start3 = 'rectangle,t_' + str(i + 1).zfill(padded_length) + ',T_' + str(i + 1).zfill(padded_length)
				connections3 = ',-' + str(square_count)
				contents.append(start3 + connections3 + '#')
		
		elif punctures == 1:
			square_count += 1
			connections = connections + ',-' + str(square_count)
			
			start2 = 'rectangle,t_' + '1'.zfill(padded_length) + ',T_' + '1'.zfill(padded_length)
			connections2 = ',+' + str(square_count)
			contents.append(start2 + connections2 + '#')
		
		elif punctures == 0 and genus > 1:
			square_count += 1
			connections = connections + ',+' + str(square_count)
		
		# Add in an extra a loop (if needed) to isolate the half-twists from the rest of the surface.
		if genus > 1 and punctures > 0:
			square_count += 1
			connections = connections + ',-' + str(square_count)
			
			start2 = 'annulus,a_' + str(punctures).zfill(padded_length) + ',A_' + str(punctures).zfill(padded_length)
			connections2 = ',+' + str(square_count)
			contents.append(start2 + connections2 + '#')
			
			# Add in the start point for the b_2 ... chain.
			square_count += 1
			connections = connections + ',+' + str(square_count)
		
		contents.append(start + connections + '#')
		
		# Now construct the rest of the b loops.
		for i in range(2, 2 * genus - 1):
			square_count += 1
			start = 'annulus,b_' + str(i).zfill(padded_length) + ',B_' + str(i).zfill(padded_length)
			connections = ',-' + str(square_count - 1) + ',+' + str(square_count)
			contents.append(start + connections + '#')
			if i == 2: c_loop = len(contents)
		
		# Don't forget, the last one is different.
		if genus > 1:
			square_count += 1
			start = 'annulus,b_' + str(2 * genus - 1).zfill(padded_length) + ',B_' + str(2 * genus - 1).zfill(padded_length)
			connections = ',-' + str(square_count - 1)
			contents.append(start + connections + '#')
		
		# Modify the 3rd b loop to add in the c loop (if needed).
		if c_loop > -1:
			next_line = 'annulus,c,C,-' + str(square_count)
			contents.append(next_line + '#')
			contents[c_loop] = contents[c_loop][:-1] + ',+' + str(square_count) + '#'
	
	# Write in the number of squares.
	contents[5] = str(square_count + 1) + '#'
	
	return '\n'.join(contents)

def code_to_sign_sequence(code):
	''' Produces a sign sequence for a given Dowker code.
	Note that this sequence may include 0's which may be filled with 
	either a +1 or -1 and so are filled with +1's by default.'''
	
	def first_non_zero(L): return min(i for i in range(len(L)) if L[i])
	
	N = len(code)
	signs = map(abs, code)
	signs = map(lambda n: n-1, signs)
	
	pairs = zip(range(0, 2*N, 2), signs)
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
				return None  # Something bad has happened, sequence is not realizable.
			
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
	
	return map(lambda n: n if n != 0 else 1, [emb[2*i] for i in range(N)])  # Note [emb[pairs_dict[2*i]] for i in range(N)] is also a valid code.

def DT_drilling_surface(code, make_prefix_unique=True):
	''' Returns a list, each entry of which is a line in a surface file
	for a surface which can be used to construct a triangulation of the
	complement of the knot with Dowker code 'code' and intersection signs 'signs'. 
	
	Example: construct_file_contents([6,10,8,12,4,2]) returns 
	a surface file that can be used to construct the 6_1 knot complement, broken 
	by line into a list.
	'''
	
	# Note: There are smaller surface files that can be constructed to do this,
	# but they are significantly more complicated and have less symmetry.
	
	signs = code_to_sign_sequence(code)
	contents = []
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
	contents.append(str(num_squares) + '#')
	
	pairs = zip(range(1, 2*num_crossings+1, 2), code)
	pairs_dict = dict([(x, abs(y)) for (x, y) in pairs] + [(abs(y), x) for (x, y) in pairs])
	signs_dict = dict([(x, y > 0) for (x, y) in pairs] + [(abs(y), y > 0) for (x, y) in pairs])
	
	drill_names = []
	handle_names = []
	
	# First build all the rectangles.
	for i in range(1, 2*num_crossings+1):
		if i % 2 == 1:
			k = i / 2
			cells = ['-' + str(k * 5 + 1), '+' + str(k * 5 + 2), '+' + str(k * 5 + 3)]
		else:
			k = pairs_dict[i] / 2
			if signs[k] == 1:
				cells = ['-' + str(k * 5 + 4), '-' + str(k * 5 + 2), '+' + str(k * 5 + 0)]
			else:
				cells = ['-' + str(k * 5 + 0), '-' + str(k * 5 + 2), '+' + str(k * 5 + 4)]
		
		over = (i % 2 == 0) ^ signs_dict[i]
		name = ('y' if (i % 2 == 0) ^ signs_dict[i] else 'x') + '_' + str(i).zfill(padded_length)
		inverse_name = name.swapcase()
		contents.append(','.join(['rectangle', name, inverse_name] + cells) + '#')
		drill_names.append('!' + name)
	
	# Then the required annuli.
	for i in range(1, 2*num_crossings+1):
		j = i % (2*num_crossings) + 1
		if i % 2 == 1:
			k = i / 2
			l = pairs_dict[j] / 2
			if signs[l] == 1:
				cells = ['-' + str(k * 5 + 3), '+' + str(l * 5 + 4)]
			else:
				cells = ['-' + str(k * 5 + 3), '-' + str(l * 5 + 4)]
		else:
			k = pairs_dict[i] / 2
			l = j / 2
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
	
	return '\n'.join(contents)

def DT_handles_surface(code, make_prefix_unique=True):
	''' Returns a list, each entry of which is a line in a surface file
	for a surface which can be used to construct a triangulation of the
	complement of the knot with Dowker code 'code'. 
	
	Example: construct_file_contents([6,10,8,12,4,2]) returns 
	a surface file that can be used to construct the 6_1 knot complement, broken 
	by line into a list.
	'''
	
	# Note: There are smaller surface files that can be constructed to do this,
	# but they are significantly more complicated and have less symmetry.
	
	signs = code_to_sign_sequence(code)  # Get the orientation of each intersection
	
	contents = []
	contents.append('#')
	contents.append('# Surface file for knot with Dowker code:')
	contents.append('#    ' + ','.join(map(str, code)))
	contents.append('#')
	contents.append('# To build this knot complement, make a Heegaard splitting with a 2-handle attached above')
	contents.append('# and below the sequence of annuli specified in the macro "h".')
	contents.append('#')
	
	num_crossings = len(code)
	
	pairs = zip(range(1, 2*num_crossings+1, 2), code)
	pairs_dict = dict([(x, abs(y)) for (x, y) in pairs] + [(abs(y), x) for (x, y) in pairs])
	signs_dict = dict([(x, y > 0) for (x, y) in pairs] + [(abs(y), y > 0) for (x, y) in pairs])
	
	overcrossing = [signs_dict[i+1] ^ (i % 2 == 1) for i in range(2*num_crossings)]
	where_switch = filter(lambda i: not signs_dict[i+1] ^ signs_dict[((i+1) % (2*num_crossings)) + 1], range(2*num_crossings))
	last_true = max(where_switch) + 1
	
	num_annuli = len(where_switch)
	squares_in_crossings = 4 * num_crossings
	squares_in_links = 2 * num_annuli
	squares_in_rectangles = 2 * num_annuli
	
	num_squares = squares_in_crossings + squares_in_links + squares_in_rectangles
	padded_length = len(str(num_annuli)) if make_prefix_unique else 0
	
	contents.append(str(num_squares) + '#')
	
	annuli_count = 0
	squares = []
	back_squares = []
	handle_names = []
	for j in range(2*num_crossings):
		i = ((j + last_true) % (2*num_crossings))
		if i % 2 == 0:
			k = i / 2  # = (i+1) / 2
			squares.append('-' + str(4*k+0))
			squares.append('+' + str(4*k+1))
			back_squares.append('+' + str(4*k+3))
			back_squares.append('-' + str(4*k+2))
		else:
			k = pairs_dict[i+1] / 2
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
			row = 'annulus,a_' + str(annuli_count).zfill(padded_length) + ',A_' + str(annuli_count).zfill(padded_length) + ','
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
			
			handle_names.append(('a' if overcrossing[i] else 'A') + '_' + str(annuli_count).zfill(padded_length))
			
			annuli_count += 1
			squares = []
			back_squares = []
	
	for i in range(num_annuli):
		row = 'rectangle,t_' + str(i).zfill(padded_length) + ',T_' + str(i).zfill(padded_length) + ','
		row += '+' + str(squares_in_crossings + squares_in_links + 2*i) + ','
		row += '-' + str(squares_in_crossings + squares_in_links + 2*i + 1)
		row += '#'
		contents.append(row)
		
	
	num_squares = num_crossings * 5
	
	contents.append('# We also give a macro for attaching the required handles.')
	contents.append('macro,h,' + '*'.join(handle_names) + '#')
	
	return '\n'.join(contents)
