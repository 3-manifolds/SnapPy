import re
valid_arc_name_characters = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_'  # Straight from Twister's global.h.

def output_error(s):
	print "Error: %s" % s
	exit(1)

def output_warning(s):
	print "Warning: % s" % s

def load_file(path):
	''' Loads a Twister surface file. '''
	try:
		f = open(path, 'r')
	except IOError:
		output_error('Unable to open "' + path + '".')
	
	contents = []
	
	parsed_line = ''
	for line in f:
		end_marker = line.find('#')
		parsed_line += line[:len(line) if end_marker == -1 else end_marker]
		if end_marker != -1 and parsed_line != '':
			contents.append(parsed_line.replace(' ', ''))
			parsed_line = ''
	
	f.close()
	return contents

def get_info(s):
	try:
		if s[0] == '+':
			return abs(int(s)), +1
		elif s[0] == '-':
			return abs(int(s)), -1
		else:
			return abs(int(s)), +1
	except ValueError:
		output_error('Invalid square information "' + s + '".')

def validate_file(path):
	''' Checks that a given file is in the correct Twister surface file format. '''
	# Could check for redundant bigons.
	contents = load_file(path)
	
	# First check the first line is the number of squares.
	try:
		num_squares = int(contents[0])
	except ValueError:
		output_error('First line must be the number of squares, not "' + contents[0] + '".')
	
	names = []
	macro_names = []
	macros = []
	squares_used = []
	
	for index, line in enumerate(contents):
		if index == 0: continue  # Skip line 0.
		
		# Check that the type is good.
		data = line.split(',')
		if len(data) < 4: output_error('Missing information on parsed line ' + str(index+1))
		if data[0] not in ['macro', 'annulus', 'rectangle']: output_error('Invalid type "' + data[0] + '" on parsed line ' + str(index+1) + '.')
		
		if data[0] == 'macro':
			macro_names.append(data[1])
			macro_names.append(data[2])
			macros.append(set([x for x in data[3].replace('!','').split('*') if x != '']))
			macros.append(set([x for x in data[3].replace('!','').split('*') if x != '']))
		else:
			names.append(data[1])
			names.append(data[2])
			for i in range(3, len(data)):
				num, sign = get_info(data[i])
				if num >= num_squares: output_error('Invalid square referenced "' + data[i] + '" on parsed line ' + str(index+1) + '.')
				squares_used.append(data[i])
	
	if len(set(squares_used)) != len(squares_used): output_error('Square used twice.')
	if len(squares_used) != 2*num_squares: output_error('Too many squares requested, only ' + str(len(squares_used) / 2) + ' out of ' + str(num_squares) + ' used.')
	
	# Check the curve (inverse) names are unique.
	if len(set(names+macro_names)) != len(names+macro_names): output_error('Duplicated name.')
	
	# And that they are all valid names (i.e. made only from valid letters).
	for name in names+macro_names: 
		if name.translate(None, valid_arc_name_characters) != '': output_error('Invalid name "' + name + '".')
	
	# Check that the macros are resolvable.
	reachable_names = set(names)
	old_reachable_names = set()
	range_len_macros = range(len(macros))
	while reachable_names != old_reachable_names:
		old_reachable_names = set(reachable_names)
		for i in range_len_macros:
			macros[i] -= reachable_names
			if macros[i] == set(): reachable_names.add(macro_names[i])
	
	for macro, macro_name in zip(macros, macro_names):
		if macro != set(): output_error('Unable to resolve macro "' + macro_name + '".')
	
	return True

def determine_info(path):
	''' Determines the genus, number of boundary components and Euler characteristic of the surface described by the specified surface file. '''
	validate_file(path)
	
	contents = load_file(path)
	num_squares = int(contents[0])
	
	# Squares look like this:
	#      N
	#  N-------E
	#  |     / |
	# W|  +/+  |E
	#  | /     |
	#  W-------S
	#      S
	
	vertices = [set([a+str(i)]) for i in range(num_squares) for a in ['N', 'E', 'S', 'W']]
	edges = [set([a+str(i)]) for i in range(num_squares) for a in ['N', 'E', 'S', 'W']]
	
	# How to glue together 2 classes.
	def glue(X, a, b):
		A = filter(lambda n: a in n, X)[0]
		B = filter(lambda n: b in n, X)[0]
		if A != B:
			C = A.union(B)
			X.remove(A)
			X.remove(B)
			X.append(C)
	
	for index, line in enumerate(contents):
		if index == 0: continue  # Skip line 0.
		data = line.split(',')
		if data[0] == 'macro': continue
		
		for i in range(len(data) - 3 if data[0] == 'annulus' else len(data) - 4):
			num, sign = get_info(data[i + 3])
			next_num, next_sign = get_info(data[((i+1) % (len(data) - 3)) + 3])
			
			glue(edges, ('E' if sign == 1 else 'N')+str(num), ('W' if next_sign == 1 else 'S')+str(next_num))
			glue(vertices, ('E' if sign == 1 else 'N')+str(num), ('N' if next_sign == 1 else 'W')+str(next_num))
			glue(vertices, ('S' if sign == 1 else 'E')+str(num), ('W' if next_sign == 1 else 'S')+str(next_num))
	
	Euler_characteristic = len(vertices) - len(edges) + num_squares  # V - E + F.
	
	edge_ends = {'N':['N', 'E'], 'S':['W', 'S'], 'E':['E', 'S'], 'W':['N', 'W']}
	
	unglued = [X for Y in filter(lambda n: len(n) == 1, edges) for X in Y]  # Get all the boundary edges by searching for ones that are unglued.
	unglued_ends = [[edge_ends[X[0]][0] + X[1:], edge_ends[X[0]][1] + X[1:]] for X in unglued]
	to_connect = [[vertices.index(filter(lambda n: X in n, vertices)[0]), vertices.index(filter(lambda n: Y in n, vertices)[0])] for X, Y in unglued_ends]
	
	boundary = 0
	while to_connect != []:
		boundary += 1
		X = to_connect.pop()
		while X[0] != X[-1]:
			for Y in to_connect:
				if X[-1] not in Y: continue
				to_connect.remove(Y)
				X.append(Y[1] if Y[0] == X[-1] else Y[0])
				break
		# Do something with X? Nah.
	
	genus = (2 - (Euler_characteristic + boundary)) / 2
	
	return genus, boundary, Euler_characteristic

def validate_command(command, path):
	
	# Check that the braces match.
	c = 0
	for letter in command:
		if letter == '(': c += 1
		elif letter == ')': c -= 1
		if c < 0: output_error('Unbalanced braces.')
	if c != 0: output_error('Unbalanced braces.')
	
	# Check every '^' is followed by a power.
	data = command.split('^')
	for d in data[1:]:
		power_ends = min(x for x in [d.find('*'), d.find('('), d.find(')'), len(d)] if x != -1)
		try:
			power = int(d[:power_ends])
		except ValueError:
			output_error('A "^" must be immediately followed by a power.')
	
	return True

def parse_command(path, command):
	''' Expands the macros and powers in a command. '''
	
	validate_command(command, path)
	validate_file(path)
	
	contents = load_file(path)
	
	names = []
	macro_names = []
	macros = {}
	inverse = {}
	
	for line in contents[1:]:
		data = line.split(',')
		if data[0] == 'macro':
			macro_names.append(data[1])
			macro_names.append(data[2])
			inverse[data[1]] = data[2]
			inverse[data[2]] = data[1]
			macros[data[1]] = data[3]
			macros[data[2]] = '---'
		else:
			inverse[data[1]] = data[2]
			inverse[data[2]] = data[1]
			names.append(data[1])
			names.append(data[2])
	
	def get_inverse(n):
		if n[0] == '!':
			return '!' + inverse[n[1:]]
		else:
			return inverse[n]
	
	# Fill in missing inverses.
	for i in macros:
		if macros[i] == '---': macros[i] = '*'.join(map(get_inverse, macros[inverse[i]].split('*')[::-1]))
	
	# A '!' must be followed by a valid_arc_name_character.
	if re.search('![^'+valid_arc_name_characters+']', command) != None: output_error('"!" must be followed by a curve name to drill.')
	
	# Make the formatting nicer.
	# Remove all spaces.
	command = command.replace(' ', '')
	# Add in extra *'s.
	command = command.replace('(', '*(').replace(')', ')*')
	# Remove any unneeded *'s.
	command = re.sub('\*+', '*', command)
	command = command.replace('(*', '(').replace('*)', ')')
	command = command.replace('*^', '^')
	command = command.strip('*')
	
	# Add in any missing parentheses.
	i = command.find('^')
	while i != -1:
		if command[i-1] != ')':
			c = 0
			for j in range(-1, i)[::-1]:
				if command[j] == ')': c += 1
				if command[j] == '(': c -= 1
				if command[j] == '*' and c == 0: break
			command = command[:j+1] + '(' + command[j+1:i] + ')' + command[i:]
		i = command.find('^', i+2)
	
	# Add in any missing powers and add in '^1'.
	i = command.find(')')
	while i != -1:
		if i+1 == len(command) or command[i+1] != '^':
			command = command[:i+1] + '^1' + command[i+1:]
		i = command.find(')', i+2)
	
	# Now expand any parentheses.
	while '(' in command:
		end = command.find(')')
		start = command.rfind('(', 0, end)
		centre = command[start+1:end]
		power_ends = min(x for x in [command.find('*', end+2), command.find('(', end+2), command.find(')', end+2), len(command)] if x != -1)
		power = int(command[end+2:power_ends])
		
		if centre == '':
			command = command[:start] + command[power_ends:]
		elif power >= 0:
			command = command[:start] + '*'.join([centre] * power) + command[power_ends:]
		else:
			command = command[:start] + '(' + '*'.join(get_inverse(x) for x in centre.split('*')[::-1] if x != '') + ')^' + str(abs(power)) + command[power_ends:]
	
	# Expand the macros in the command using repeated substitution.
	old_command = ''
	while command != old_command:
		old_command = command
		for i in macros:
			command = command.replace(i, macros[i])
	
	# Remove any unneeded *'s.
	command = re.sub('\*+', '*', command)
	command = command.strip('*')
	
	return command

if __name__ == "__main__":
	import sys
	args = sys.argv
	print parse_command(args[1], args[2])
	genus, boundary, Euler_characteristic = determine_info(args[1])
	
	print "This file defines a surface with:"
	print "\tGenus: ", genus
	print "\tBoundary components: ", boundary
	print "\tEuler characteristic: ", Euler_characteristic
	

	