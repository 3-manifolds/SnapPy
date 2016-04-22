import math
import tkFileDialog
import snappy

def export_stl(self, model='klein', cutout=False,
               subdivisions = 5, shrink_factor = .9, cutoff_radius=.9):
	"""
	export_stl provides the same function as the export option in the 
	file menu of the Browser or Dirichlet menus. Takes arguments to 
	specify the model that the stl file will produce. Shrink Factor
	effects the hollowness of the cutout version of the models.
	"""
	facedict = self.face_list()
	Model = model.lower()
	shrink = shrink_factor
	if 0 < shrink < 1:
		pass
	else:
		raise Exception('Shrink factor must be between 0 and 1')
	if Model=='klein' or Model=='poincare':
		pass
	else:
		raise Exception('Model must be assigned either Klein or Poincare')
	if cutout == False:
		export_norm_stl(facedict, Model, subdivisions)
	elif cutout == True:
		export_cutout_stl(facedict, Model, shrink, subdivisions)
	else:
		raise Exception('cutout must be assigned either True or False')	

def export_norm_stl(facedict, model, subdivisions):		
	if model=='klein':
		klein_to_stl(facedict)
	elif model=='poincare':
		poincare_to_stl(facedict, subdivisions)
		
def export_cutout_stl(facedict, model, shrink, subdivisions):			
	if model=='klein':
		klein_cutout(facedict, shrink, subdivisions)
	elif model == 'poincare':
		poincare_cutout(facedict, shrink, subdivisions)	
    
def facet_stl(file,vertex1,vertex2,vertex3):
	a = (vertex3[0]-vertex1[0], vertex3[1]-vertex1[1],vertex3[2]-vertex1[2])
	b = (vertex2[0]-vertex1[0], vertex2[1]-vertex1[1],vertex2[2]-vertex1[2])
	normal = (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
	file.write('  facet normal %f %f %f\n' %normal)
	file.write('    outer loop\n')
	file.write('      vertex %f %f %f\n' %vertex1)
	file.write('      vertex %f %f %f\n' %vertex2)
	file.write('      vertex %f %f %f\n' %vertex3)
	file.write('    endloop\n')
	file.write('  endfacet\n')

def tri_div(triangles):
	new_triangles=[]
	for triangle in triangles:
		x=triangle[0]
		y=triangle[1]
		z=triangle[2]
		xy=midpoint(x,y)
		yz=midpoint(y,z)
		zx=midpoint(z,x)
		t1=[x,xy,zx]
		t2=[xy,yz,zx]
		t3=[zx,yz,z]
		t4=[xy,y,yz]
		new_triangles.append(t1)
		new_triangles.append(t2)
		new_triangles.append(t3)
		new_triangles.append(t4)
	triangles=new_triangles
	return triangles

def midpoint(vertex1,vertex2):
	x1=vertex1[0]
	x2=vertex2[0]
	y1=vertex1[1]
	y2=vertex2[1]
	z1=vertex1[2]
	z2=vertex2[2]
	midpoint= ((x1+x2)/2,(y1+y2)/2,(z1+z2)/2)
	return midpoint	
	
def transform(vertex,cutoff_radius=.9):
	x=vertex[0]
	y=vertex[1]
	z=vertex[2]
	D=x**2+y**2+z**2
	scale=1/(1+math.sqrt(max(0,1-D)))
	if scale >= cutoff_radius:
		scale= cutoff_radius
	xp = scale*x
	yp = scale*y
	zp = scale*z
	p_vertex=(xp,yp,zp)
	return p_vertex
	
def klein_to_stl(facedict):
	f = tkFileDialog.asksaveasfile(
    mode='w',
    title='Save Klein Model as STL file',
    defaultextension = '.stl',
    filetypes = [
        ("STL files", "*.stl *.STL", ""),
        ("All files", "")])
	f.write('solid\n')
	klein_faces = facedict
	for face in klein_faces:
		vertices = face['vertices']
		for i in range(len(vertices)-2):
			vertex1 = vertices[0]
			vertex2 = vertices[i+1]
			vertex3 = vertices[i+2]
			facet_stl(f,vertex1,vertex2,vertex3)
	f.write('endsolid')
	f.close()	
	
def poincare_to_stl(facedict, subdivisions):
	f = tkFileDialog.asksaveasfile(
    mode='w',
    title='Save Poincare Model as STL file',
    defaultextension = '.stl',
    filetypes = [
        ("STL files", "*.stl *.STL", ""),
        ("All files", "")])
	f.write('solid\n')
	klein_faces = facedict
	for face in klein_faces:
		vertices = face['vertices']
		for i in range(len(vertices)-2):
			v1 = vertices[0]
			v2 = vertices[i+1]
			v3 = vertices[i+2]
			triangle = [v1,v2,v3]
			triangles = []
			triangles.append(triangle)
			for i in range(subdivisions):
				triangles = tri_div(triangles)
			for triangle in triangles:
				Vertex1 = triangle[0]
				Vertex2 = triangle[1]
				Vertex3 = triangle[2]
				vertex1 = transform(Vertex1)
				vertex2 = transform(Vertex2)
				vertex3 = transform(Vertex3)
				facet_stl(f,vertex1,vertex2,vertex3)
	f.write('endsolid')
	f.close()
	
def klein_cutout(facedict,shrink):
	f = tkFileDialog.asksaveasfile(
    mode='w',
    title='Save Klein Model Cutout as STL file',
    defaultextension = '.stl',
    filetypes = [
        ("STL files", "*.stl *.STL", ""),
        ("All files", "")])
	f.write('solid\n')
	klein_faces = facedict
	point_list = []
	for face in klein_faces:
		vertices = face['vertices']
		center=[0,0,0]
		for vertex in vertices:
			x = vertex[0]
			y = vertex[1]
			z = vertex[2]
			center = [center[0]+x,center[1]+y,center[2]+z]	
		c1 = center[0]/len(vertices)
		c2 = center[1]/len(vertices)
		c3 = center[2]/len(vertices)
		center = [c1,c2,c3]
		new_vertices = []
		for vertex in vertices:
			x=vertex[0]
			y=vertex[1]
			z=vertex[2]
			dir_vec = [((c1-x)/3),((c2-y)/3),((c3-z)/3)]
			x0=x+dir_vec[0]
			y0=y+dir_vec[1]
			z0=z+dir_vec[2]
			new_vertex=(x0,y0,z0)
			new_vertices.append(new_vertex)
		new_inside_points = []
		for point in new_vertices:
			p1 = point[0]*shrink
			p2 = point[1]*shrink
			p3 = point[2]*shrink
			new_point=(p1,p2,p3)
			new_inside_points.append(new_point)
		for i in range(len(new_vertices)):
			vertex1 = new_vertices[i]
			if i!=len(new_vertices)-1:
				vertex2 = new_inside_points[i+1]
			else:
				vertex2 = new_inside_points[0]
			vertex3 = new_inside_points[i]
			facet_stl(f,vertex1,vertex2,vertex3)
		for i in range(len(new_vertices)):
			vertex1 = new_vertices[i]
			if i!=len(new_vertices)-1:
				vertex2 = new_vertices[i+1]
				vertex3 = new_inside_points[i+1]
			else:
				vertex2 = new_vertices[0]
				vertex3 = new_inside_points[0]
			facet_stl(f,vertex1,vertex2,vertex3)
		for i in range(len(vertices)):
			vertex1 = vertices[i]
			if i!=len(vertices)-1:
				vertex2 = new_vertices[i+1]
			else:
				vertex2 = new_vertices[0]
			vertex3 = new_vertices[i]
			facet_stl(f,vertex1,vertex2,vertex3)
			point_list.extend([vertex1,vertex2,vertex3])
		for i in range(len(vertices)):
			vertex1 = vertices[i]
			if i!=len(vertices)-1:
				vertex2 = vertices[i+1]
				vertex3 = new_vertices[i+1]
			else:
				vertex2 = vertices[0]
				vertex3 = new_vertices[0]
			facet_stl(f,vertex1,vertex2,vertex3)
			point_list.extend([vertex1,vertex2,vertex3])
	new_points=[]
	for point in point_list:
		p1 = point[0]*shrink
		p2 = point[1]*shrink
		p3 = point[2]*shrink
		new_point=(p1,p2,p3)
		new_points.append(new_point)
	for i in range(0,len(new_points)-1,3):
		vertex1=new_points[i]
		vertex2=new_points[i+1]
		vertex3=new_points[i+2]
		facet_stl(f,vertex1,vertex2,vertex3)
	f.write('endsolid')
	f.close()
		
def poincare_cutout(facedict, shrink, subdivisions):
	f = tkFileDialog.asksaveasfile(
    mode='w',
    title='Save Poincare Model Cutout as STL file',
    defaultextension = '.stl',
    filetypes = [
        ("STL files", "*.stl *.STL", ""),
        ("All files", "")])
	f.write('solid\n')
	klein_faces = facedict
	point_list = []
	for face in klein_faces:
		vertices = face['vertices']
		center=[0,0,0]
		for vertex in vertices:
			x = vertex[0]
			y = vertex[1]
			z = vertex[2]
			center = [center[0]+x,center[1]+y,center[2]+z]	
		c1 = center[0]/len(vertices)
		c2 = center[1]/len(vertices)
		c3 = center[2]/len(vertices)
		center = [c1,c2,c3]
		new_vertices = []
		for vertex in vertices:
			x=vertex[0]
			y=vertex[1]
			z=vertex[2]
			dir_vec = [((c1-x)/3),((c2-y)/3),((c3-z)/3)]
			x0=x+dir_vec[0]
			y0=y+dir_vec[1]
			z0=z+dir_vec[2]
			new_vertex=(x0,y0,z0)
			new_vertices.append(new_vertex)
		new_points = new_vertices
		for j in range(subdivisions):
			midpoints = [new_points[0]]
			for i in range(len(new_points)):
				if i!=len(new_points)-1:
					mid = midpoint(new_points[i], new_points[i+1])
					midpoints.extend([mid,new_points[i+1]])
				else:
					mid = midpoint(new_points[0], new_points[i])
					midpoints.extend([mid])
			new_points = midpoints
		for i in range(len(new_points)):
			new_points[i] = transform(new_points[i])
		new_inside_points = []
		for point in new_points:
			p1 = point[0]*shrink
			p2 = point[1]*shrink
			p3 = point[2]*shrink
			new_point=(p1,p2,p3)
			new_inside_points.append(new_point)
		for i in range(len(new_points)):
			vertex1 = new_points[i]
			if i!=len(new_points)-1:
				vertex2 = new_inside_points[i+1]
			else:
				vertex2 = new_inside_points[0]
			vertex3 = new_inside_points[i]
			facet_stl(f,vertex1,vertex2,vertex3)
		for i in range(len(new_points)):
			vertex1 = new_points[i]
			if i!=len(new_points)-1:
				vertex2 = new_points[i+1]
				vertex3 = new_inside_points[i+1]
			else:
				vertex2 = new_points[0]
				vertex3 = new_inside_points[0]
			facet_stl(f,vertex1,vertex2,vertex3)
		for i in range(len(vertices)):
			v1 = vertices[i]
			if i!=len(vertices)-1:
				v2 = new_vertices[i+1]
			else:
				v2 = new_vertices[0]
			v3 = new_vertices[i]
			triangle = [v1,v2,v3]
			triangles = []
			triangles.append(triangle)
			for i in range(subdivisions):
				triangles = tri_div(triangles)
			for triangle in triangles:
				Vertex1 = triangle[0]
				Vertex2 = triangle[1]
				Vertex3 = triangle[2]
				vertex1 = transform(Vertex1)
				vertex2 = transform(Vertex2)
				vertex3 = transform(Vertex3)
				facet_stl(f,vertex1,vertex2,vertex3)
				point_list.extend([vertex1,vertex2,vertex3])
		for i in range(len(vertices)):
			v1 = vertices[i]
			if i!=len(vertices)-1:
				v2 = vertices[i+1]
				v3 = new_vertices[i+1]
			else:
				v2 = vertices[0]
				v3 = new_vertices[0]
			triangle = [v1,v2,v3]
			triangles = []
			triangles.append(triangle)
			for i in range(subdivisions):
				triangles = tri_div(triangles)
			for triangle in triangles:
				Vertex1 = triangle[0]
				Vertex2 = triangle[1]
				Vertex3 = triangle[2]
				vertex1 = transform(Vertex1)
				vertex2 = transform(Vertex2)
				vertex3 = transform(Vertex3)
				facet_stl(f,vertex1,vertex2,vertex3)
				point_list.extend([vertex1,vertex2,vertex3])
	new_points=[]
	for point in point_list:
		p1 = point[0]*shrink
		p2 = point[1]*shrink
		p3 = point[2]*shrink
		new_point=(p1,p2,p3)
		new_points.append(new_point)
	for i in range(0,len(new_points)-1,3):
		vertex1=new_points[i]
		vertex2=new_points[i+1]
		vertex3=new_points[i+2]
		facet_stl(f,vertex1,vertex2,vertex3)
	f.write('endsolid')
	f.close()

snappy.DirichletDomain.export_stl = export_stl