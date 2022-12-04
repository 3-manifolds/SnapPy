import tkinter
from tkinter import ttk

from .gui_utilities import UniformDictController, ScrollableFrame
from .geodesics import geodesic_index_to_color
from ..drilling.exceptions import WordAppearsToBeParabolic

class GeodesicsWindow(tkinter.Toplevel):
    def __init__(self, inside_viewer, *args, **kwards):
        self.inside_viewer = inside_viewer
        self.raytracing_view = inside_viewer.widget

        tkinter.Toplevel.__init__(self, class_='snappy')
        self.title('Geodesics')

        self.frame = ttk.Frame(self)
        self.frame.pack(expand = True, fill = tkinter.BOTH)
        self.frame.columnconfigure(0, weight = 1)

        top_frame = ttk.Frame(self.frame)
        top_frame.pack()

        left_top_frame = ttk.Frame(top_frame)
        left_top_frame.pack(side = tkinter.LEFT, padx = 20)

        self.length_button = ttk.Button(
            left_top_frame, text = "Add up to length", command=self.add_length_spectrum)
        self.length_button.grid(row = 0, column = 0)

        self.length_var = tkinter.StringVar(value=1.0)

        self.length_box = ttk.Spinbox(
            left_top_frame,
            from_=0.2, to=20.0, increment=0.2,
            textvariable = self.length_var,
            width = 4)
        self.length_box.grid(row = 0, column = 1)

        right_top_frame = ttk.Frame(top_frame)
        right_top_frame.pack(side = tkinter.LEFT, padx = 20)

        self.word_button = ttk.Button(
            right_top_frame, text = "Add word", command = self.add_word)
        self.word_button.grid(row = 0, column = 0)
        self.word_entry = ttk.Entry(
            right_top_frame)
        self.word_entry.grid(row = 0, column = 1)

        self.status_label = ttk.Label(self.frame)
        self.status_label.pack()

        self.scrollable_frame = ScrollableFrame(self.frame)
        self.scrollable_frame.pack(expand = True, fill = tkinter.BOTH)

        self.geodesics_frame = self.scrollable_frame.scrollable_frame
        self.geodesics_frame.columnconfigure(0, weight = 0)
        self.geodesics_frame.columnconfigure(1, weight = 0)
        self.geodesics_frame.columnconfigure(2, weight = 0)
        self.geodesics_frame.columnconfigure(3, weight = 0)
        self.geodesics_frame.columnconfigure(4, weight = 0)
        self.geodesics_frame.columnconfigure(5, weight = 1)
        self.geodesics_frame.columnconfigure(6, weight = 0)

        self.populate_geodesics_frame()

    def populate_geodesics_frame(self):
        for widget in self.geodesics_frame.grid_slaves():
            widget.destroy()

        row = 0

        checkbox_column = 0
        color_column = 1
        words_column = 2
        length_column = 3
        radius_column = 5

        l = ttk.Label(self.geodesics_frame, text = 'Show')
        l.grid(row = row, column = checkbox_column)

        l = ttk.Label(self.geodesics_frame, text = 'Color')
        l.grid(row = row, column = color_column)

        l = ttk.Label(self.geodesics_frame, text = 'Word(s)')
        l.grid(row = row, column = words_column)

        l = ttk.Label(self.geodesics_frame, text = 'Complex length')
        l.grid(row = row, column = length_column, columnspan = 2)

        l = ttk.Label(self.geodesics_frame, text = 'Radius')
        l.grid(row = row, column = radius_column)

        row += 1

        for geodesic in self.raytracing_view.geodesics.geodesics_sorted_by_length():
            UniformDictController.create_checkbox(
                self.geodesics_frame,
                self.raytracing_view.ui_parameter_dict,
                key = 'geodesicTubeEnables',
                index = geodesic.index,
                row = row,
                column = checkbox_column,
                update_function = self.geodesic_checkbox_clicked)

            l = ttk.Label(self.geodesics_frame, text = ', '.join(geodesic.words))
            l.grid(row = row, column = words_column)

            l = ttk.Label(self.geodesics_frame,
                          text = '%.8f' % geodesic.complex_length.real())
            l.grid(row = row, column = length_column)

            im_length = geodesic.complex_length.imag()
            abs_im_length = im_length.abs()

            if abs_im_length > 1e-10:
                s = '+' if im_length > 0 else '-'

                l = ttk.Label(self.geodesics_frame,
                              text = s + ' %.8f * I' % abs_im_length)
                l.grid(row = row, column = length_column + 1)


            color = geodesic_index_to_color(geodesic.index)

            l = tkinter.Label(self.geodesics_frame,
                              text = "Color",
                              fg = color_to_tkinter(color),
                              bg = color_to_tkinter(color))
            l.grid(row = row, column = color_column, padx = 5)

            scale = UniformDictController.create_horizontal_scale(
                self.geodesics_frame,
                self.raytracing_view.ui_parameter_dict,
                key = 'geodesicTubeRadii',
                index = geodesic.index,
                row = row,
                column = radius_column,
                left_end = 0.0,
                right_end = 1.0,
                update_function = self.raytracing_view.update_geodesic_data_and_redraw,
                format_string = '%.3f').scale

            # Need to color Scale - but the following code fails.
            # scale.configure(background = color_to_tkinter(color))

            row += 1

    def add_length_spectrum(self):
        self.status_label.configure(text = "")

        self.raytracing_view.geodesics.add_length_spectrum(
            float(self.length_box.get()))

        self.raytracing_view.resize_geodesic_params()

        self.populate_geodesics_frame()

    def add_word(self):
        word = self.word_entry.get()
        try:
            index = self.raytracing_view.geodesics.add_word(word)
        except WordAppearsToBeParabolic:
            self.status_label.configure(text = word + " is parabolic")
            return

        self.status_label.configure(text = "")
        
        self.raytracing_view.resize_geodesic_params()
        self.raytracing_view.enable_geodesic(index)
        if self.raytracing_view.disable_edges_for_geodesics():
            self.inside_viewer.update_edge_and_insphere_controllers()

        self.raytracing_view.update_geodesic_data_and_redraw()

        self.populate_geodesics_frame()

    def geodesic_checkbox_clicked(self):
        if self.raytracing_view.disable_edges_for_geodesics():
            self.inside_viewer.update_edge_and_insphere_controllers()
        self.raytracing_view.update_geodesic_data_and_redraw()

def color_to_tkinter(color):
    return "#%.3x%.3x%.3x" % tuple([min(max(int(x * 4095), 0), 4095)
                                    for x in color])
