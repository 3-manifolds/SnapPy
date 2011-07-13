import os, sys, re, tempfile, time
import snappy.CyTwister

def construct_surface_file_contents(genus, punctures):
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
  
def write_surface_to_file(genus, punctures, file):
    file_contents = construct_surface_file_contents(genus, punctures)
    open(file, "w").write("\n".join(file_contents) + "\n")

# ----------------------------------

class TwisterError(Exception):
    pass

def twister_create_file( surface=None, monodromy=None, gluing=None, handles=None,
             name=None, peripheral_curves=True, optimize=True, warnings=True):
    """
    Creates a Manifold or Triangulation using Twister
    """

    tempfiles = []
    if isinstance(surface, basestring):
        surface_file_name = os.path.abspath(surface)
    else:
        genus, punctures = surface 
        surface_file_name = tempfile.mktemp() + ".sur"
        tempfiles.append(surface_file_name)
        write_surface_to_file(genus, punctures, surface_file_name)

    if not os.path.isfile(surface_file_name):
        raise ValueError, "Surface file %s can't be found" % surface_file_name

    if monodromy != None and (gluing != None or handles != None):
        raise ValueError, "Specify *either* a bundle *or* a heegaard splitting"

    tri_file_name, err_file_name = tempfile.mktemp() + ".tri", tempfile.mktemp() + ".err"
    tempfiles += [err_file_name]
    call_args = ["Twister.out", "--surface", surface_file_name, "--name", tri_file_name,
                 "--output", tri_file_name, "-ow", err_file_name]
    
    if monodromy:
        call_args += ["--bundle", monodromy]
    else:
        description = gluing if gluing else ""
        handles = handles if handles else ""
        call_args += ["--splitting", description, "--handles", handles]

    if not peripheral_curves:
        call_args.append( "-ml" ) 

    if not optimize:
        call_args.append( "--optimisations" )

    if not warnings:
        call_args.append( "--warnings" )

    snappy.CyTwister.call_twister(call_args)

    try:
        errs = open(err_file_name).read()
    except IOError:
        errs = ""

    if errs:
        err = errs.split("\n")[0]
        raise TwisterError, err 

        
    for file in tempfiles:
        if os.path.isfile(file):
            os.remove(file)

    return tri_file_name

def twister( surface=None, monodromy=None, gluing=None, handles=None,
             name=None, with_hyperbolic_structure=True,
             peripheral_curves=True, optimize=True, warnings=True):

    tri_file_name = twister_create_file(surface, monodromy, gluing, handles,
                                       name, peripheral_curves, optimize, warnings)
    
    ans = snappy.Manifold(tri_file_name) if with_hyperbolic_structure else snappy.Triangulation(tri_file_name)

    if name:
        ans.set_name(name)
    else:
        ans.set_name("Twister Manifold")

    os.remove(tri_file_name)
    return ans

def test_twister():
    M = twister(surface = (1, 1), monodromy="a_0*B_1")
    print M.volume(), M.fundamental_group().relators()

    M = twister(surface="../Twister/Surfaces/S_1_1.sur", monodromy="a*B")
    print M.volume(), M.fundamental_group().relators()

    M = twister(surface = (1, 1), handles="a_0*B_1")
    print M.fundamental_group().generators(), M.fundamental_group().relators()

    monodromy=  "*".join(10*["!a_0*!b_1"]) 
    M = twister(surface = (1,1),monodromy=monodromy)
    print M.num_cusps(), M.volume()

    M = twister(surface = (1, 1), monodromy="a_0*B_1*C", peripheral_curves=False, optimize=False, warnings=False)
    print M.volume(), M.fundamental_group().relators()

    M = twister(surface = (1, 1), monodromy="a_0*B_1", with_hyperbolic_structure=False)
    print type(M)

    M = twister("4braid.sur", gluing="", handles="e*E")
    print M.num_cusps(), M.fundamental_group()

    M = twister("4braid.sur", gluing="b*c*a*a*B*a*B*B", handles="e*E")
    print M.num_cusps(), M.is_two_bridge()

    M = twister("4braid.sur", gluing="B*a*a*B*B*a*B*B", handles="e*E")
    print M.num_cusps(), M.is_two_bridge()
                   
bundle_strings = [
    "Bundle( S(2,1) , [a_0, B_1, a_1,!b_2])",
    "Bundle( S(1,2),a_0*B_1*a_1*!b_1)",
    "Bundle( S(1,12),[a_0*B_1*a_2])",
    ]

bundle_pat = re.compile("Bundle\(S\((\d+),(\d+)\),\[*([abcABC_\d!,*]*)\]*\)")

def bundle_from_string(desc):
    desc = desc.replace(' ', '')
    m = bundle_pat.match(desc)
    if m:
        g, n, monodromy = m.groups()
        g, n = int(g), int(n)
        monodromy = monodromy.replace(",", "*")
        tri_file = twister_create_file(surface=(g,n), monodromy=monodromy)


splitting_pat = re.compile("Splitting\(S\((\d+),(\d+)\),\[*([abcABC_\d!,*]*)\]*,*\[*([abcABC_\d!,*]*)\]*\)")

splitting_strings = [
    "Splitting(S(2,0),[b_1*B_2*c*b_1*!c,b_2,A_0,C,B_2,b_3,b_2,c], a_0*c*B_3])",
    "Splitting(S(2,0),[b_1*b_2*b_3*c*a_0*b_2*b_2*b_3*b_1*a_0*B_1*A_0*c*c*b_1*b_3*B_2*b_1*b_2*b_3*c*a_0*b_2*b_2*b_3*b_1*a_0*B_1*A_0*c*c*b_1*b_3*B_2*b_1*b_2*b_3*c*a_0*b_2*b_2*b_3*b_1*a_0*B_1*A_0*c*c*b_1*b_3*B_2*b_1*b_2*b_3*c*a_0*b_2*b_2*b_3*b_1*a_0*B_1*A_0*c*c*b_1*b_3*B_2*b_1*b_2*b_3*c*a_0*b_2*b_2*b_3*b_1*a_0*B_1*A_0*c*c*b_1*b_3*B_2], [a_0*c*C])"
    ]

test = "Bundle(S(2,0),[b_1*b_2*b_3*c*a_0*b_2*b_2*b_3*b_1*a_0*B_1*A_0*c*c*b_1*b_3*B_2])"
    
def splitting_from_string(desc):
    desc = desc.replace(' ', '')
    m = splitting_pat.match(desc)
    if m:
        g, n, gluing, handles = m.groups()
        g, n = int(g), int(n)
        gluing, handles = gluing.replace(",", "*"), handles.replace(",", "*")
        print "glue", gluing
        print "handles", handles
        tri_file = twister_create_file(surface=(g,n),gluing=gluing, handles=handles)

if __name__ == "__main__":
    test_twister()
