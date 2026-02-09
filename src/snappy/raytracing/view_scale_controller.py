import math

_parameter_interval = [
    (20.0,120.0),    # material: field of view
    (-2.5,2.5),      # ideal: log of scale
    None ]           # hyperideal

_scale_function = [
    lambda parameter: math.tan(parameter / 360.0 * math.pi),
    lambda parameter: math.exp(parameter),
    None ]

_inverse_scale_function = [
    lambda scale: math.atan(scale) * 360.0 / math.pi,
    lambda scale: math.log(scale),
    None ]


class ViewScaleController:
    def __init__(self,
                 uniform_dict,
                 scale,
                 label0,
                 label1,
                 update_function):
        self.uniform_dict = uniform_dict
        self.scale = scale
        self.label0 = label0
        self.label1 = label1
        self.update_function = update_function

        self.update()

        self.scale.set_callback(self.scale_command)

    def scale_command(self, value):
        perspective_type = self.uniform_dict['perspectiveType'][1]
        parameter_interval = _parameter_interval[perspective_type]
        if not parameter_interval:
            return

        parameter = _linear_remap(
            value,
            (self.scale.left_end, self.scale.right_end),
            parameter_interval)

        self.uniform_dict['viewScale'][1] = float(
            _scale_function[perspective_type](parameter))
        self.update_label()
        self.update_function()

    def update_label(self):
        perspective_type = self.uniform_dict['perspectiveType'][1]

        if perspective_type == 0:
            view_scale = self.uniform_dict['viewScale'][1]
            fov = _inverse_scale_function[perspective_type](view_scale)
            self.label0.configure(text='FOV:')
            self.label1.configure(text='%.1f' % fov)
        elif perspective_type == 1:
            view_scale = self.uniform_dict['viewScale'][1]
            self.label0.configure(text='Euclidean length:')
            if view_scale > 1.0:
                self.label1.configure(text='%.1f' % view_scale)
            elif view_scale > 0.1:
                self.label1.configure(text='%.2f' % view_scale)
            else:
                self.label1.configure(text='%.3f' % view_scale)
        else:
            self.label0.configure(text='')
            self.label1.configure(text='')

    def update_scale(self):
        perspective_type = self.uniform_dict['perspectiveType'][1]
        parameter_interval = _parameter_interval[perspective_type]
        if not parameter_interval:
            self.scale.state(("disabled",))
            return

        self.scale.state(("!disabled",))
        view_scale = self.uniform_dict['viewScale'][1]
        parameter = _inverse_scale_function[perspective_type](view_scale)

        self.scale.set_value(
            _linear_remap(
                parameter,
                parameter_interval,
                (self.scale.left_end, self.scale.right_end)))

    def update(self):
        self.update_label()
        self.update_scale()


def _linear_remap(v, src_interval, dst_interval):
    l0, r0 = src_interval
    l1, r1 = dst_interval
    return l1 + (r1 - l1) / (r0 - l0) * (v - l0)
