
def construct_file_contents(code, signs):
	''' Returns a list, each entry of which is a line in a surface file
	for a surface which can be used to construct a triangulation of the
	complement of the knot with Dowker code 'code' and intersection signs 'signs'. 
	
	Example: construct_file_contents([6,10,8,12,4,2], [-1,+1,+1,-1,+1,+1]) returns 
	a surface file that can be used to construct the 6_1 knot complement, broken 
	by line into a list.
	'''
	
	# Note: There are more efficient surface files that can be constructed to do this,
	# but they are significantly more complicated and have less symmetry.
	contents = []
	contents.append('#')
	contents.append('# Surface file for knot with:')
	contents.append('# Dowker code: '+ ','.join(map(str, code)))
	contents.append('# Sign code: '+ ','.join(map(str, signs)))
	contents.append('#')
	contents.append('# To build this knots complement, make a Heegaard splitting with a 2-handle attached above')
	contents.append('# and below every annulus and drill every rectangle exactly once, making sure to drill all')
	contents.append('# of the "x" rectangles before drilling ANY of the "y" rectangles.')
	contents.append('#')
	
	num_crossings = len(code)
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
		name = ('y' if (i % 2 == 0) ^ signs_dict[i] else 'x') + '_' + str(i)
		drill_names.append('!' + name)
		inverse_name = ('Y' if (i % 2 == 0) ^ signs_dict[i] else 'X') + '_' + str(i)
		row = ','.join(['rectangle', name, inverse_name] + cells) + '#'
		contents.append(row)
	
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
		
		name = 'a' + '_' + str(i)
		inverse_name = 'A' + '_' + str(i)
		handle_names.append(name)
		handle_names.append(inverse_name)
		row = ','.join(['annulus', name, inverse_name] + cells) + '#'
		contents.append(row)
	
	contents.append('# We also give a macro to drill all of the rectangles in an acceptable order.')
	contents.append('# And one for attaching all handles too.')
	drill_names.sort()  # Make sure to sort the rectangle drill names so that all of the "x"'s appear before ANY of the "y"'s.
	contents.append('macro,s,S,' + '*'.join(drill_names) + '#')
	contents.append('macro,t,T,' + '*'.join(handle_names) + '#')
	
	return contents

def write_surface_to_file(contents, file):
	''' Writes the provided contents to the specified file.'''
	f = open(file, 'w')
	for line in contents:
		f.write(line + '\n')
	
	f.close()

if __name__ == "__main__":
	import sys
	code, signs = map(int, sys.argv[1].split(',')), map(int, sys.argv[2].split(','))
	write_surface_to_file(construct_file_contents(code, signs), 'Knot_S_0_' + str(len(code)) + '.sur')
	
	# The WRONG way to find the signs: try all of them.
	# Although this works, it is exponential in the number of crossings!
	
	# ------------------------------------------------------------------
	# code = map(int, sys.argv[1].split(','))
	
	# from Twister import determine_info
	# for i in range(len(code)):
		# s = [1 if (i == x) else -1 for x in range(len(code))]
		# write_surface_to_file(construct_file_contents(code, s), 'Knot_S_0_' + str(len(code)) + '.sur')
		# genus, boundary, Euler_characteristic = determine_info('Knot_S_0_' + str(len(code)) + '.sur')
		# print genus, s
		
	
	# for i in range(2**(len(code)-1)):
		# s = [1 if (i & 2**x != 0) else -1 for x in range(len(code))]
		# write_surface_to_file(construct_file_contents(code, s), 'Knot_S_0_' + str(len(code)) + '.sur')
		# genus, boundary, Euler_characteristic = determine_info('Knot_S_0_' + str(len(code)) + '.sur')
		# print genus, s
		# if genus == 0:
			# print '!!!', s  # Write out the correct sequence.
			# break
	# ------------------------------------------------------------------
	
