import tkinter
from tkinter import ttk

from .gui_utilities import UniformDictController

class GeodesicsWindow(tkinter.Toplevel):
    def __init__(self, inside_viewer, *args, **kwards):
        self.inside_viewer = inside_viewer
        self.raytracing_view = inside_viewer.widget

        tkinter.Toplevel.__init__(
            self,
            master = inside_viewer,
            class_='snappy')

        self.title('Geodesics')

        self.frame = ttk.Frame(self)
        self.frame.pack(expand = True, fill = tkinter.BOTH)
        self.frame.columnconfigure(0, weight = 1)

        length_frame = ttk.Frame(self.frame)
        length_frame.grid(row = 0, column = 0)

        self.length_button = ttk.Button(
            length_frame, text = "Add up to length", command=self.add_length_spectrum)
        self.length_button.grid(row = 0, column = 0)

        self.length_var = tkinter.StringVar(value=1.0)

        self.length_box = ttk.Spinbox(
            length_frame,
            from_=0.2, to=20.0, increment=0.2,
            textvariable = self.length_var,
            width = 4)
        self.length_box.grid(row = 0, column = 1)

        self.word_button = ttk.Button(
            length_frame, text = "Add word", command = self.add_word)
        self.word_button.grid(row = 0, column = 2)
        self.word_entry = ttk.Entry(
            length_frame)
        self.word_entry.grid(row = 0, column = 3)

        self.geodesics_frame = None

        self.populate_geodesics_frame()

    def populate_geodesics_frame(self):
        if not self.geodesics_frame is None:
            self.geodesics_frame.destroy()

        self.geodesics_frame = ttk.Frame(self.frame)
        self.geodesics_frame.grid(row = 1, column = 0, sticky = tkinter.NSEW)

        self.geodesics_frame.columnconfigure(0, weight = 0)
        self.geodesics_frame.columnconfigure(1, weight = 0)
        self.geodesics_frame.columnconfigure(2, weight = 0)
        self.geodesics_frame.columnconfigure(3, weight = 0)
        self.geodesics_frame.columnconfigure(4, weight = 1)
        self.geodesics_frame.columnconfigure(5, weight = 0)

        row = 0

        checkbox_column = 0
        words_column = 1
        length_column = 2
        color_column = 3
        radius_column = 4

        l = ttk.Label(self.geodesics_frame, text = 'Show')
        l.grid(row = row, column = checkbox_column)

        l = ttk.Label(self.geodesics_frame, text = 'Word(s)')
        l.grid(row = row, column = words_column)

        l = ttk.Label(self.geodesics_frame, text = 'Complex length')
        l.grid(row = row, column = length_column)

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

            l = ttk.Label(self.geodesics_frame, text = str(geodesic.complex_length))
            l.grid(row = row, column = length_column)

            l = ttk.Label(self.geodesics_frame, text = "Color %d" % geodesic.index)
            l.grid(row = row, column = color_column)

            UniformDictController.create_horizontal_scale(
                self.geodesics_frame,
                self.raytracing_view.ui_parameter_dict,
                key = 'geodesicTubeRadii',
                index = geodesic.index,
                row = row,
                column = radius_column,
                left_end = 0.0,
                right_end = 1.0,
                update_function = self.raytracing_view.update_geodesic_data_and_redraw,
                format_string = '%.3f')

            row += 1

    def add_length_spectrum(self):

        self.raytracing_view.geodesics.add_length_spectrum(
            float(self.length_box.get()))

        self.raytracing_view.resize_geodesic_params()

        self.populate_geodesics_frame()

    def add_word(self):
        index = self.raytracing_view.geodesics.add_word(
            self.word_entry.get())

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
