
def construct_file_contents(genus, punctures):
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
	
	square_count = 0
	
	if genus == 0:
		if punctures == 0:
			square_count += 1
			contents.append('annulus,a_1,A_1,-0,+1#')
			contents.append('annulus,a_2,A_2,-1,+0#')
		else:
			for i in range(punctures):
				square_count += 1
				start = 'rectangle,t_' + str(i + 1) + ',T_' + str(i + 1)
				connections = ',-' + str(square_count - 1) + ',+' + str(square_count)
				contents.append(start + connections + '#')
				
				if i == punctures - 1:
					start = 'annulus,a_' + str(i + 1) + ',A_' + str(i + 1)
					connections = ',-' + str(square_count) + ',+' + '0'
					contents.append(start + connections + '#')
				else:
					square_count += 1
					start = 'annulus,a_' + str(i + 1) + ',A_' + str(i + 1)
					connections = ',-' + str(square_count - 1) + ',+' + str(square_count)
					contents.append(start + connections + '#')
	else:
		c_loop = -1  # Just a marker.
		
		# Start with the a_0 loop.
		next_line = 'annulus,a_0,A_0,+0'
		contents.append(next_line + '#')
		
		# And the b_1 loop.
		start = 'annulus,b_1,B_1'
		connections = ',-0'
		
		if punctures > 1:
			for i in range(punctures - 1):
				square_count += 1
				connections = connections + ',-' + str(square_count)
				
				square_count += 1
				
				start2 = 'annulus,a_' + str(i + 1) + ',A_' + str(i + 1)
				connections2 = ',+' + str(square_count - 1) + ',+' + str(square_count)
				contents.append(start2 + connections2 + '#')
				
				start3 = 'rectangle,t_' + str(i + 1) + ',T_' + str(i + 1)
				connections3 = ',-' + str(square_count)
				contents.append(start3 + connections3 + '#')
		
		elif punctures == 1:
			square_count += 1
			connections = connections + ',-' + str(square_count)
			
			start2 = 'rectangle,t_1,T_1'
			connections2 = ',+' + str(square_count)
			contents.append(start2 + connections2 + '#')
		
		elif punctures == 0 and genus > 1:
			square_count += 1
			connections = connections + ',+' + str(square_count)
		
		
		# Add in an extra a loop (if needed) to isolate the half-twists from the rest of the surface.
		if genus > 1 and punctures > 0:
			square_count += 1
			connections = connections + ',-' + str(square_count)
			
			start2 = 'annulus,a_' + str(punctures) + ',A_' + str(punctures)
			connections2 = ',+' + str(square_count)
			contents.append(start2 + connections2 + '#')
			
			# Add in the start point for the b_2 ... chain.
			square_count += 1
			connections = connections + ',+' + str(square_count)
		
		contents.append(start + connections + '#')
		
		# Now construct the rest of the b loops.
		for i in range(2, 2 * genus - 1):
			# print contents
			square_count += 1
			start = 'annulus,b_' + str(i) + ',B_' + str(i)
			connections = ',-' + str(square_count - 1) + ',+' + str(square_count)
			contents.append(start + connections + '#')
			if i == 2: c_loop = len(contents)
		
		# Don't forget, the last one is different.
		if genus > 1:
			square_count += 1
			start = 'annulus,b_' + str(2 * genus - 1) + ',B_' + str(2 * genus - 1)
			connections = ',-' + str(square_count - 1)
			contents.append(start + connections + '#')
		
		# Modify the 3rd b loop to add in the c loop (if needed).
		if c_loop > -1:
			next_line = 'annulus,c,C,-' + str(square_count)
			contents.append(next_line + '#')
			contents[c_loop] = contents[c_loop][:-1] + ',+' + str(square_count) + '#'
	
	# Write in the number of squares.
	contents[5] = str(square_count + 1) + '#'
	
	return contents

def write_surface_to_file(contents, file):
	''' Writes the provided contents to the specified file.'''
	f = open(file, 'w')
	for line in contents:
		f.write(line + '\n')
	
	f.close()

if __name__ == "__main__":
	import sys
	genus, punctures = int(sys.argv[1]), int(sys.argv[2])
	
	write_surface_to_file(construct_file_contents(genus, punctures), 'LP_S_' + str(genus) + '_' + str(punctures) + '.sur')