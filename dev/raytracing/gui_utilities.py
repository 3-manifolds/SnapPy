import tkinter as Tk_
import tkinter.ttk as ttk

import time

def attach_scale_and_label_to_uniform(uniform_dict,
                                      key, update_function, scale, label,
                                      format_string = None,
                                      index = None):
    uniform_type, value = uniform_dict[key]

    if index is None:
        if not uniform_type in ['int', 'float']:
            raise Exception("Unsupported uniform slider type")
    else:
        if uniform_type == 'int[]':
            uniform_type = 'int'
        elif uniform_type == 'float[]':
            uniform_type = 'float'
        elif uniform_type == 'vec2':
            uniform_type = 'float'
        else:
            raise Exception("Unsupported uniform slider type")

        value = value[index]

    if format_string is None:
        if uniform_type == 'int':
            format_string = '%d'
        else:
            format_string = '%.2f'

    scale.configure(value = value)
    label.configure(text = format_string % value)
    
    def scale_command(t,
                      uniform_dict = uniform_dict,
                      key = key,
                      format_string = format_string,
                      uniform_type = uniform_type,
                      update_function = update_function,
                      label = label,
                      index = index):

        value = float(t)
        if uniform_type == 'int':
            value = int(value)
        if index is None:
            uniform_dict[key] = (uniform_type, value)
        else:
            uniform_dict[key][1][index] = value
        label.configure(text = format_string % value)
        if update_function:
            update_function()

    scale.configure(command = scale_command)

    def update_scale(uniform_dict = uniform_dict,
                     key = key,
                     format_string = format_string,
                     uniform_type = uniform_type,
                     scale = scale,
                     label = label,
                     index = index):

        uniform_type, value = uniform_dict[key]
        if not index is None:
            value = value[index]
        if uniform_type == 'int':
            value = int(value)
        else:
            value = float(value)
        label.configure(text = format_string % value)
        scale.configure(value = value)

    return update_scale

def create_horizontal_scale_for_uniforms(
    window, uniform_dict, key, title, row, from_, to, update_function,
    index = None,
    format_string = None,
    column = 0):

    if title:
        title_label = ttk.Label(window, text = title)
        title_label.grid(row = row, column = column, sticky = Tk_.NSEW)
        column += 1

    scale = ttk.Scale(window, from_ = from_, to = to,
                      orient = Tk_.HORIZONTAL)
    scale.grid(row = row, column = column, sticky = Tk_.NSEW)
    column += 1

    value_label = ttk.Label(window)
    value_label.grid(row = row, column = column, sticky = Tk_.NSEW)
    
    return attach_scale_and_label_to_uniform(
        uniform_dict, key, update_function, scale, value_label,
        index = index,
        format_string = format_string)

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
    
