from GLSLManifoldInsideView import *

from GLSLManifoldInsideView import _constant_uniform_bindings

class DirichletViewWidget(SimpleImageShaderWidget):
    def __init__(self, dirichlet_domain, master, *args, **kwargs):
        
        self.ui_uniform_dict = {
            'maxSteps' : ('int', 20),
            'maxDist' : ('float', 17),
            'subpixelCount': ('int', 1),
            'fov': ('float', 90),
            'edgeThickness': ('float', 0.005),
            'edgeThicknessCylinder' : ('float', 1.01)
            }

        self.ui_parameter_dict = {
            'insphere_scale' : ('float', 0.05),
            'cuspAreas' : ('float[]', manifold.num_cusps() * [ 1.0 ]),
            'translationVelocity' : ('float', 0.4),
            'rotationVelocity' : ('float', 0.4)
            }

        self.dirchlet_domain = dirchlet_domain
        
        self.current_key_pressed = None
        self.mouse = None
        self.boost_mouse = None

        self._initialize_raytracing_data()

        
