
import snappy

def export_stl(self, filename, klein_vs_poincare = model, shrink_factor = .9,  cutout=False):
	"""
	Documentation goes in here.
	"""
	print(shrink_factor)
	print(cutout)
	face_dicts = self.face_list()
	#f = open(filename,'w')
	#f.write(str(face_dicts))
	#f.close()	4
	if model.lower()=='klein':
		pass
	elif model.lower()=='poincare':
		pass
	else:
		raise Exception('Thats really stupid')

snappy.DirichletDomain.export_stl = export_stl