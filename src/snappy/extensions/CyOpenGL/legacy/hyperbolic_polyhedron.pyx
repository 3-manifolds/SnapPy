class HyperbolicPolyhedron:
    """
    A hyperbolic polyhedron for display in OpenGL, either in the Klein
    model or the Poincare model.  Includes a representation of the
    sphere at infinity.  It is initialized with the SnapPea description
    of the faces of a Dirichlet domain, represented as a list of
    dictionaries.
    """

    def __init__(self, facedicts, model_var, sphere_var, togl_widget=None):
        self.facedicts = facedicts
        self.model = model_var
        self.sphere = sphere_var
        self.face_specular = [0.5, 0.5, 0.5, 1.0]
        self.front_shininess = 50.0
        self.back_shininess = 50.0
        self.S_infinity = WireframeSphere(color=[1.0, 1.0, 1.0, .2],
                                          front_specular=[0.5, 0.5, 0.5, 1.0],
                                          front_shininess=50.0,
                                          togl_widget=togl_widget)
        self.S_infinity.build_display_list(1.0, 30, 30)
        self.Klein_faces = []
        self.Poincare_faces = []
        for dict in facedicts:
            vertices = [vector3(vertex) for vertex in dict['vertices']]
            closest = vector3(dict['closest'])
            center = closest*(1/dict['distance']**2)
            color = hls_to_rgb(dict['hue'], 0.5, 1.0) + (1.0,)
            self.Klein_faces.append(
                KleinPolygon(vertices, closest,
                             color=color,
                             front_specular=self.face_specular,
                             back_specular=self.face_specular,
                             front_shininess=self.front_shininess,
                             back_shininess=self.back_shininess,
                             togl_widget=togl_widget))
            self.Poincare_faces.append(
                PoincarePolygon(vertices, center,
                                color=color,
                                front_specular=self.face_specular,
                                back_specular=self.face_specular,
                                front_shininess=self.front_shininess,
                                back_shininess=self.back_shininess,
                                togl_widget=togl_widget))
        for face in self.Klein_faces:
            face.build_display_list()
        for face in self.Poincare_faces:
            face.build_display_list()

    def draw(self, *args):
        model = self.model.get()
        if model == 'Klein':
            for face in self.Klein_faces:
                face.display()
        elif model == 'Poincare':
            for face in self.Poincare_faces:
                face.display()
        if self.sphere.get():
            self.S_infinity.display()
