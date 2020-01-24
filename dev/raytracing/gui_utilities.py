import tkinter as Tk_
import tkinter.ttk as ttk

import time

class UniformDictController:
    @staticmethod
    def create_horizontal_scale(container, uniform_dict, key,
                                row, from_, to, update_function,
                                column = 0,
                                title = None,
                                format_string = None,
                                index = None, component_index = None):
        if title:
            title_label = ttk.Label(container, text = title)
            title_label.grid(row = row, column = column, sticky = Tk_.NSEW)
            column += 1
        
        scale = ttk.Scale(container, from_ = from_, to = to,
                      orient = Tk_.HORIZONTAL)
        scale.grid(row = row, column = column, sticky = Tk_.NSEW)
        column += 1

        value_label = ttk.Label(container)
        value_label.grid(row = row, column = column, sticky = Tk_.NSEW)

        return UniformDictController(
            uniform_dict, key, scale, value_label, update_function,
            format_string = format_string,
            index = index, component_index = component_index)

    def __init__(self, uniform_dict, key,
                 scale, label,
                 update_function,
                 format_string = None,
                 index = None, component_index = None):

        self.uniform_dict = uniform_dict
        self.key = key
        self.scale = scale
        self.label = label
        self.update_function = update_function
        self.index = index
        self.component_index = component_index

        self.uniform_type, value = self.uniform_dict[self.key]
        
        if self.uniform_type in ['int', 'float']:
            self.scalar_type = self.uniform_type
            if not (index is None and component_index is None):
                raise Exception("int/float uniform does not support index")
        elif self.uniform_type == 'float[]':
            self.scalar_type = 'float'
            if index is None or not component_index is None:
                raise Exception("Need to specify index for float[] uniform")
        elif self.uniform_type == 'vec2[]':
            self.scalar_type = 'float'
            if index is None or component_index is None:
                raise Exception("Need to specify indices for vec2[] uniform")
        else:
            raise Exception("Unsupported uniform type %s" % self.uniform_type)
        
        if format_string:
            self.format_string = format_string
        else:
            if self.scalar_type == 'int':
                self.format_string = '%d'
            else:
                self.format_string = '%.2f'

        if self.scale:
            self.scale.configure(command = self.scale_command)
        self.update_scale_and_label()
        
    def get_value(self):
        if self.uniform_type == 'int':
            return int(self.uniform_dict[self.key][1])
        if self.uniform_type == 'float':
            return float(self.uniform_dict[self.key][1])
        if self.uniform_type == 'float[]':
            return float(self.uniform_dict[self.key][1][self.index])
        if self.uniform_type == 'vec2[]':
            return float(
                self.uniform_dict[self.key][1][self.index][self.component_index])
                 
    def set_value(self, value):
        if self.uniform_type == 'int':
            self.uniform_dict[self.key][1] = int(float(value))
        elif self.uniform_type == 'float':
            self.uniform_dict[self.key][1] = float(value)
        elif self.uniform_type == 'float[]':
            self.uniform_dict[self.key][1][self.index] = float(value)
        elif self.uniform_type == 'vec2[]':
            self.uniform_dict[self.key][1][self.index][self.component_index] = (
                float(value))

    def update_scale(self):
        if self.scale:
            self.scale.configure(value = self.get_value())

    def update_label(self):
        if self.label:
            self.label.configure(text = self.format_string % self.get_value())

    def update_scale_and_label(self):
        self.update_scale()
        self.update_label()

    def scale_command(self, value):
        self.set_value(value)
        self.update_label()
        if self.update_function:
            self.update_function()

class FpsLabelUpdater:
    def __init__(self, label):
        self.label = label
        self.num_iterations = 0
        self.total_time = 0.0
        self.last_time = time.time()

    def __call__(self, t):
        self.num_iterations += 1
        self.total_time += t
        if self.num_iterations > 50 or time.time() - self.last_time > 2.0:
            current_time = time.time()
            fps = self.num_iterations / (current_time - self.last_time)
            time_ms = 1000 * self.total_time / self.num_iterations

            self.label.configure(text = '%.1ffps (%dms)' % (fps, time_ms))
            self.last_time = current_time
            self.num_iterations = 0
            self.total_time = 0.0
    
