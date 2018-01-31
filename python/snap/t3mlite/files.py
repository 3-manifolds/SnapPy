#$Id: files.py,v 1.6 2010/01/18 20:33:11 t3m Exp $
#   t3m - software for studying triangulated 3-manifolds
#   Copyright (C) 2002 Marc Culler, Nathan Dunfield and others
#
#   This program is distributed under the terms of the 
#   GNU General Public License, version 2 or later, as published by
#   the Free Software Foundation.  See the file GPL.txt for details.

from .arrow import eArrow
from .simplex import *
from .tetrahedron import Tetrahedron
import os, sys, re

# Nathan's code for importing and exporting snappea files.
# Converts a SnapPea file to MComplex.  Doesn't really use all the
# structure of the SnapPea file as it relies only on the fact that the
# gluing data for the ith pair of tetrahedra is given by the ith pair
# of lines like:
#
#      2    5    1   34 
#   3120 0321 0132 0132

def read_SnapPea_file(file_name=None, data = None):
    if data is None: 
        data = open(file_name).read().decode('ascii')
    count = 0

    neighbors_match = "^\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s*$"
    perm_match = "\s*([0123]{4,4})\s+([0123]{4,4})\s+([0123]{4,4})\s+([0123]{4,4})\s*$"
    snappea_re = re.compile(neighbors_match + perm_match, re.MULTILINE)
    
    fake_tets =[]
    
    curr_poss = 0
    while 1:
        m = snappea_re.search(data, curr_poss)
        if not m:
            break
        else:
            neighbors = [int(g) for g in m.group(1,2,3,4)]
            perms = []
            for perm in m.group(5,6,7,8):
                perm = [int(p) for p in [perm[0], perm[1], perm[2], perm[3]]]
                perms.append(perm)
            fake_tets.append( (neighbors, perms) )
            curr_poss = m.end(8)
    return fake_tets

#------------End function SnapPea to Mcomplex--------------------


# Exports an Mcomplex in SnapPea 2.0 format.
# ASSUMES THAT THE MANIFOLD IS ORIENTABLE AND THAT THE LINK OF
# ANY VERTEX HAS GENUS AT MOST ONE.

def write_SnapPea_file(mcomplex, fileobject):
    out = fileobject.write
    if hasattr(fileobject, 'name'):
        name = fileobject.name
    else:
        name = 'untitled'

    out("% Triangulation\n\n" + name + "\nnot_attempted 0.0\nunknown_orientability\nCS_unknown\n\n")

    torus_cusps = []
    for vertex in mcomplex.Vertices:
        g = vertex.link_genus()
        if g > 1:
            raise ValueError("Link of vertex has genus more than 1.")
        if g == 1:
            torus_cusps.append(vertex)

    # All torus cusps are unfilled
    
    out("%d 0\n" % len(torus_cusps))
    for i in torus_cusps:
        out( "   torus   0.000000000000   0.000000000000\n" )

    out("\n")

    # The num of tetrahedra

    out("%d\n" % len(mcomplex))
    
    # Output the tetraheda themselves.
    
    for tet in mcomplex.Tetrahedra:
        for face in TwoSubsimplices:
            out("    %d" % mcomplex.Tetrahedra.index( tet.Neighbor[face]))
        out("\n")
        for face in TwoSubsimplices:
            out(" %d%d%d%d" % tet.Gluing[face].tuple())

        out("\n")
        for vert in ZeroSubsimplices:
            vertex = tet.Class[vert]
            if vertex.link_genus() == 1:
                out("%d " % torus_cusps.index(vertex))
            else:
                out("-1 ")
        out("\n")
        for i in range(4):
            out("0 0 0 0  0 0 0 0   0 0 0 0   0 0 0 0\n")
        out("0.0 0.0\n\n")


# Nathan's code for importing and exporting geo files.
#
# converts Casson's u, v, w, x into 0, 1, 2, 3

conv = {"u" : V0, "v" : V1, "w" : V2, "x" : V3}
conv_back = {V0: "u", V1 : "v", V2 : "w", V3: "x"}

# Geo saves edges in the form "5ux".  This function returns
# (tet_num, starting_vertex, ending_vertex).  Because Casson
# starts his indexing of tets at 1 and Jeff starts at 0,
# we subtract 1.

def read_edge(edge):
    m = re.match("([0-9]+)([uvwx])([uvwx])", edge)
    return (int(m.group(1)) - 1, conv[m.group(2)], conv[m.group(3)])

#  Geo stores manifolds by storing the link around each edge.  This
#  function takes two successive edges in the link and glues the
#  corresponding tetrahedra together.

def read_geo_file(file_name, num_tet=None):
    data = open(file_name).readlines()
    if num_tet == None:
        num_tet = len(data) - 2
    tets = []
    for i in range(num_tet):
        tets.append(Tetrahedron())

    for line in data[1  : ]:
        line = line.decode('ascii')
        cycle = re.split("\s+", line[ : -1])[1 : ]
        for i in range(len(cycle)):
            t1, v1, v2 = read_edge(cycle[i])
            t2, w1, w2 = read_edge(cycle[(i+1)%len(cycle)]) #Yes, that's w2, w1
            a = eArrow(tets[t1], v1, v2)
            b = eArrow(tets[t2], w1, w2)
            a.glue(b)

    return Mcomplex(tets)

#---------Code to go from Mcomplex to Geo---------------------

def write_geo_file(mcomplex, fileobject):
    out = fileobject.write
    out("k\n")
    i = 1
    for edge in mcomplex.Edges:
        tet = edge.Corners[0].Tetrahedron
        edge_name = edge.Corners[0].Subsimplex
        init = Head[edge_name]
        fin = Tail[edge_name]
        a = eArrow(tet, init, fin).opposite()
        b = a.copy()
        out("%d\t%d%s%s " % (i, mcomplex.Tetrahedra.index(b.Tetrahedron) + 1,
                             conv_back[b.tail()], conv_back[b.head()]))
        b.next()
        while b != a:
            out("%d%s%s " % (mcomplex.Tetrahedra.index(b.Tetrahedron) + 1,
                             conv_back[b.tail()], conv_back[b.head()]))
            b.next()

        i = i + 1
        out("\n")

# writing a file for Matveev's program Spine

def write_spine_file(mcomplex, fileobject):
    out = fileobject.write
    for edge in mcomplex.Edges:
        n = edge.valence()
        A = edge.get_arrow()
        tets, global_faces, local_faces, back_local_faces = [], [], [], []
        for i in range(n):
            tets.append(A.Tetrahedron.Index + 1)
            global_faces.append(A.Tetrahedron.Class[A.Face].Index + 1)
            local_faces.append(A.Face)
            back_local_faces.append(comp(A.head()))
            A.next()

            
        signs = [1 if (tets[i], local_faces[i]) < (tets[(i + 1) % n], back_local_faces[(i + 1)%n]) else -1 for i in range(n)]
        ans= repr([signs[i]*global_faces[i] for i in range(n)])[1:-1].replace(",", "")
        out(ans + "\n")

    
        
            
        



__all__  = ('read_SnapPea_file',
            'write_SnapPea_file',
            'read_geo_file',
            'write_geo_file',
            'write_spine_file',
            )
