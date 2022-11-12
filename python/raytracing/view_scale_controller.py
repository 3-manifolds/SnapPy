import math

_material_fov_min = 20.0
_material_fov_max = 120.0
_material_fov_default = 90.0

class ViewScaleController:
    def __init__(self,
                 uniform_dict,
                 parameter_dict,
                 scale,
                 label0,
                 label1,
                 update_function):
        self.uniform_dict = uniform_dict
        self.parameter_dict = parameter_dict
        self.scale = scale
        self.label0 = label0
        self.label1 = label1
        self.update_function = update_function

        self.scale.set_callback(self.scale_command)

        default_value = _linear_remap(
            _material_fov_default,
            _material_fov_min, _material_fov_max,
            self.scale.left_end, self.scale.right_end)

        self.scale.set_value(default_value)
        self.scale_command(default_value)
                                                           
    def scale_command(self, value):
        fov = _linear_remap(
            value,
            self.scale.left_end, self.scale.right_end,
            _material_fov_min, _material_fov_max)

        self.label1.configure(text = '%.1f' % fov)
        
        self.parameter_dict['viewScaleParameter'][1] = float(fov)
        self.uniform_dict['viewScale'][1] = float(math.tan(fov / 360.0 * math.pi))
        self.update_function()
        

def _linear_remap(v, l0, r0, l1, r1):
    return l1 + (r1 - l1) / (r0 - l0) * (v - l0)
