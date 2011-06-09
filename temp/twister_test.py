import os, sys, re, tempfile
import snappy.CyTwister


def construct_file_contents(genus, punctures):
	''' Returns a list, each entry of which is a line in the standard surface file for the surface S_{genus, punctures}
		according to the naming convention described in figure 13 of Labruere and Paris' paper "Presentations for the 
		punctured mapping class groups in terms of Artin groups".
		
		In the case where genus == 1, the loop a_n is dropped as it is isotopic to the loop a_0.
		
		The case where genus == 0 is handled seperately.
		
		Example: construct_file_contents(3,2) returns the standard surface file for S_{3,2} broken by line into a list.
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
			square_count = square_count + 1
			contents.append('annulus,a_1,A_1,-0,+1#')
			contents.append('annulus,a_2,A_2,-1,+0#')
		else:
			for i in range(punctures):
				square_count = square_count + 1
				start = 'rectangle,t_' + str(i + 1) + ',T_' + str(i + 1)
				connections = ',-' + str(square_count - 1) + ',+' + str(square_count)
				contents.append(start + connections + '#')
				
				if i == punctures - 1:
					start = 'annulus,a_' + str(i + 1) + ',A_' + str(i + 1)
					connections = ',-' + str(square_count) + ',+' + '0'
					contents.append(start + connections + '#')
				else:
					square_count = square_count + 1
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
		
		for i in range(punctures - 1):
			square_count = square_count + 1
			connections = connections + ',-' + str(square_count)
			
			square_count = square_count + 1
			
			start2 = 'annulus,a_' + str(i + 1) + ',A_' + str(i + 1)
			connections2 = ',+' + str(square_count - 1) + ',+' + str(square_count)
			contents.append(start2 + connections2 + '#')
			
			start3 = 'rectangle,t_' + str(i + 1) + ',T_' + str(i + 1)
			connections3 = ',-' + str(square_count)
			contents.append(start3 + connections3 + '#')
		
		if punctures == 1:
			square_count = square_count + 1
			connections = connections + ',-' + str(square_count)
			
			start2 = 'rectangle,t_1,T_1'
			connections2 = ',+' + str(square_count)
			contents.append(start2 + connections2 + '#')
		
		if genus > 1:
			# Add in an extra a loop (if needed) to isolate the half-twists from the rest of the surface.
			square_count = square_count + 1
			connections = connections + ',-' + str(square_count)
			
			start2 = 'annulus,a_' + str(punctures) + ',A_' + str(punctures)
			connections2 = ',+' + str(square_count)
			contents.append(start2 + connections2 + '#')
			
			# Add in the start point for the b_2 ... chain.
			square_count = square_count + 1
			connections = connections + ',+' + str(square_count)
		
		contents.append(start + connections + '#')
		
		# Now construct the rest of the b loops.
		for i in range(2, 2 * genus - 1):
			square_count = square_count + 1
			start = 'annulus,b_' + str(i) + ',B_' + str(i)
			connections = ',-' + str(square_count - 1) + ',+' + str(square_count)
			contents.append(start + connections + '#')
			if i == 2: c_loop = len(contents)
		
		# Don't forget the last one is different.
		if genus > 1:
			square_count = square_count + 1
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

	
def write_surface_to_file(genus, punctures, file):
    file_contents = construct_file_contents(genus, punctures)
    open(file, "w").write("\n".join(file_contents) + "\n")
        
def triangulation_from_twister(genus, punctures, monodromy=None, gluing=None, handles=None):
    if genus < 0 or punctures < 0:
        raise ValueError, "Invaild topological description of a surface"
    if monodromy != None and (gluing != None or handles != None):
        raise ValueError, "Specify *either* a bundle *or* a heegaard splitting"

    surface_file_name = tempfile.mktemp() + ".sur"
    tri_file_name = tempfile.mktemp() + ".tri"

    write_surface_to_file(genus, punctures, surface_file_name)

    if monodromy:
        manifold_type = "--bundle"
        description = monodromy
        tail = []
    else:
        manifold_type = "--splitting"
        description = gluing if gluing else ""
        if not handles:
            handles = ""
        tail = ["--handles", handles]

    snappy.CyTwister.call_twister(
        ["Twister.out", "--surface", surface_file_name,
         manifold_type, description, "--name", tri_file_name] + tail)

    M = snappy.Manifold(tri_file_name)
    os.remove(surface_file_name)
    os.remove(tri_file_name)
    return M

while 1:
	M = triangulation_from_twister(1, 1, monodromy="a_0*B_1")
	print M.volume(), M.fundamental_group().relators()
	M = triangulation_from_twister(1, 1, handles="a_0*B_1")
	print M.fundamental_group().generators(), M.fundamental_group().relators()
	monodromy= "[" + "*".join(10*["a_0*b_1"]) + "]"
	M = triangulation_from_twister(1,1,monodromy=monodromy)
	print M.num_cusps(), M.volume()



    
