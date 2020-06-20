#!/usr/bin/env python
from builtins import range
from .gui import *
from .CyOpenGL import (HoroballScene, OpenGLOrthoWidget,
                       GetColor, cyglSetStandardLighting)
from plink.ipython_tools import IPythonTkRoot
import os, sys

class HoroballViewer(ttk.Frame):
    def __init__(self, master, nbhd=None, which_cusp=0, cutoff=None,
                 title='Horoball Viewer',
                 prefs={'cusp_horoballs' : True,
                        'cusp_triangulation' : True,
                        'cusp_ford_domain' : True,
                        'cusp_labels' : True,
                        'cusp_parallelogram' : True,
                        'cusp_cutoff' : '0.1000'},
                 bgcolor=None,
                 main_window=None):
        self.prefs = prefs
        ttk.Frame.__init__(self, master)
        self.nbhd = nbhd
        self.empty = (self.nbhd is None)
        self.mouse_x = 0
        self.mouse_y = 0
        self.menubar = None
        self.main_window = main_window
        if cutoff is None:
            self.cutoff = float(prefs['cusp_cutoff'])
        else:
            self.cutoff = float(cutoff)
        self.which_cusp = which_cusp
        self.moving_cusp = 0
        self.cusp_moving = False
        self.last_slider_value = None
        self.busy_drawing = False
        self.title = title
        self.style = style = SnapPyStyle()
        self.bgcolor = bgcolor if bgcolor else style.groupBG
        self.pgram_var = pgram_var = Tk_.IntVar(self,
            value=prefs['cusp_parallelogram'])
        self.Ford_var = Ford_var = Tk_.IntVar(self,
            value=prefs['cusp_ford_domain'])
        self.tri_var = tri_var = Tk_.IntVar(self,
            value=prefs['cusp_triangulation'])
        self.horo_var = horo_var = Tk_.IntVar(self,
            value=prefs['cusp_horoballs'])
        self.label_var = label_var = Tk_.IntVar(self,
            value=prefs['cusp_labels'])
        self.flip_var = flip_var = Tk_.BooleanVar(self)
        self.columnconfigure(0, weight=1)
        self.rowconfigure(1, weight=1)
        self.top_frame = top_frame = ttk.Frame(self)
        top_frame.columnconfigure(1, weight=1)
        self.bottomframe = bottomframe = ttk.Frame(self)
        self.widget = widget = OpenGLOrthoWidget(master=bottomframe,
            width=600, height=500, fovy=3.0, depth=1, double=True, swapinterval=0,
            help = """
Use the mouse to drag the scene relative to the fundamental parallelogram.

Use the sliders to adjust the sizes of the horoballs. Color coding indicates who bumps whom.

To change the cutoff size, enter a number in the box and hit return.

Cusps which are "tied" change size in unison.

To view the scene from outside of the upper half-space, check the "Flip" checkbutton.

Use the View Options to select which components of the scene are drawn.
""")
        widget.bind('<ButtonPress-1>', self.click)
        widget.bind('<B1-Motion>', self.translate)
        if self.horo_var.get():
            widget.set_background(0.3, 0.3, 0.4)
        else:
            widget.set_background(1.0,1.0,1.0)
        widget.autospin_allowed = 0
        cyglSetStandardLighting()
        option_frame= ttk.Frame(top_frame)
        view_button = ttk.Menubutton(option_frame, text='View Options')
        self.view_menu = view_menu = Tk_.Menu(view_button, tearoff=0)
        view_menu.add_checkbutton(label='parallelogram', command=self.view_check,
                                variable=self.pgram_var)
        view_menu.add_checkbutton(label='Ford edges', command=self.view_check,
                                variable=self.Ford_var)
        view_menu.add_checkbutton(label='triangulation', command=self.view_check,
                                variable=self.tri_var)
        view_menu.add_checkbutton(label='horoballs', command=self.view_check,
                                variable=self.horo_var)
        view_menu.add_checkbutton(label='labels', command=self.view_check,
                                variable=self.label_var)
        view_button.config(menu=view_menu)
        view_button.grid(row=0, column=0, columnspan=2, sticky=Tk_.W, padx=0, pady=0)
        flip_button = ttk.Checkbutton(option_frame, text='Flip',
                                      variable = self.flip_var,
                                      takefocus=False,
                                      command=self.flip)
        flip_button.grid(row=1, column=0, sticky=Tk_.W, padx=0, pady=0)
        self.cutoff_label = ttk.Label(option_frame, text='Cutoff: ')
        self.cutoff_var = cutoff_var = Tk_.StringVar(self,
            value='%.4f'%self.cutoff)
        self.cutoff_entry = ttk.Entry(option_frame, width=6, takefocus=False,
                                      textvariable=cutoff_var)
        self.cutoff_entry.bind('<Return>', self.set_cutoff)
        self.cutoff_label.grid_forget()
        self.cutoff_entry.grid_forget()
        self.cutoff_label.grid(row=2, column=0, sticky=Tk_.EW)
        self.cutoff_entry.grid(row=2, column=1, sticky=Tk_.W, padx=(0,20), pady=2)
        self.slider_frame = slider_frame = ttk.Frame(top_frame)
        self.eye_label = ttk.Label(slider_frame, text='Eye')
        self.tie_label = ttk.Label(slider_frame, text='Tie')
        if self.nbhd and self.nbhd.num_cusps() > 1:
            self.eye_label.grid(row=0, column=0, sticky=Tk_.W, pady=0)
            self.tie_label.grid(row=0, column=1, sticky=Tk_.W, pady=0)
        ttk.Label(slider_frame, text='Cusp Position').grid(
            row=0, column=2, pady=0)
        ttk.Label(slider_frame, text='Volume').grid(
            row=0, column=3, pady=0, padx=0, sticky=Tk_.W)
        self.eye_var = Tk_.IntVar(self, value=self.which_cusp)
        self.cusp_sliders = []
        self.slider_frames = []
        self.tie_vars = []
        self.tie_buttons = []
        self.eye_buttons = []
        self.volume_labels = []
        slider_frame.grid_columnconfigure(0, weight=1)
        slider_frame.grid_columnconfigure(2, minsize=370, weight=0)
        self.build_sliders()
        option_frame.grid(row=0, column=0, padx=(10,5), pady=5)
        slider_frame.grid(row=0, column=1, padx=5, pady=(5,10))
        top_frame.grid(row=0, column=0, sticky=Tk_.NSEW, padx=0, pady=0)
        zoomframe = ttk.Frame(bottomframe)
        self.zoom = zoom = ttk.Scale(zoomframe, from_=0, to=100,
            orient=Tk_.VERTICAL, command=self.set_zoom)
        zoom.set(30)
        zoom.pack(side=Tk_.TOP, expand=Tk_.YES, fill=Tk_.Y)
        bottomframe.columnconfigure(0, weight=1)
        bottomframe.rowconfigure(0, weight=1)
        widget.grid(row=0, column=0, sticky=Tk_.NSEW)
        zoomframe.grid(row=0, column=1, sticky=Tk_.NS)
        bottomframe.grid(row=1, column=0, sticky=Tk_.NSEW)
        self.update_idletasks()
        self.build_menus()
        self.scene = HoroballScene(nbhd, pgram_var, Ford_var, tri_var,
            horo_var, label_var, flipped=self.flip_var.get(), cutoff=self.cutoff,
            which_cusp=self.which_cusp,togl_widget=self.widget)
        self.widget.redraw_impl = self.scene.draw
        if isinstance(master, Tk_.Toplevel):
            master.config(menu=self.menubar)
            # hacks needed on Sierra
            self.after(20, self.configure_sliders)
            self.after(50, self.rebuild)
        else:
            self.configure_sliders()

    def apply_prefs(self, prefs):
        for key in self.prefs:
            value = prefs.get(key, 'missing')
            if value != 'missing':
                self.prefs[key] = value
        self.pgram_var.set(prefs['cusp_parallelogram'])
        self.Ford_var.set(prefs['cusp_ford_domain'])
        self.tri_var.set(prefs['cusp_triangulation'])
        self.horo_var.set(prefs['cusp_horoballs'])
        self.label_var.set(prefs['cusp_labels'])
        self.cutoff = float(prefs['cusp_cutoff'])
        self.cutoff_var.set('%.4f'%self.cutoff)
        self.rebuild()

    def view_check(self):
        if self.horo_var.get():
            self.widget.set_background(0.3, 0.3, 0.4)
        else:
            self.widget.set_background(1.0, 1.0, 1.0)
        self.widget.redraw_if_initialized()

    def build_sliders(self):
        nbhd = self.nbhd
        if nbhd is None:
            return
        self.cusp_vars = []
        self.cusp_colors = []
        self.tie_vars = []
        num_cusps = nbhd.num_cusps()
        if num_cusps > 1:
            self.eye_label.grid(row=0, column=0, sticky=Tk_.E, pady=0)
            self.tie_label.grid(row=0, column=1, sticky=Tk_.E, pady=0)
        else:
            self.eye_label.grid_forget()
            self.tie_label.grid_forget()
        for n in range(num_cusps):
            disp = float(nbhd.stopping_displacement(which_cusp=n))
            nbhd.set_displacement(disp, which_cusp=n)
            if nbhd and nbhd.num_cusps() > 1:
                eye_button = ttk.Radiobutton(self.slider_frame, text='',
                    variable=self.eye_var, takefocus=False, value=n,
                    command=self.set_eye)
                self.eye_buttons.append(eye_button)
                eye_button.grid(row=n+1, column=0)
                tie_var = Tk_.IntVar(self)
                tie_var.set(nbhd.get_tie(n))
                self.tie_vars.append(tie_var)
                tie_button = ttk.Checkbutton(self.slider_frame, variable=tie_var,
                    takefocus=False, command=self.rebuild)
                tie_button.index = n
                tie_button.grid(row=n+1, column=1)
                self.tie_buttons.append(tie_button)
            R, G, B, A = GetColor(nbhd.original_index(n))
            self.cusp_colors.append('#%.3x%.3x%.3x'%(
                int(R*4095), int(G*4095), int(B*4095)))
            self.cusp_vars.append(Tk_.IntVar(self))
            self.slider_frames.append(Tk_.Frame(self.slider_frame, borderwidth=0))
            self.slider_frames[n].grid(row=n+1, column=2, sticky=Tk_.EW,
                                       padx=6, pady=1)
            slider = Tk_.Scale(self.slider_frames[n],
                               showvalue=0, from_=0, to=100,
                               width=11, length=200, orient=Tk_.HORIZONTAL,
                               background=self.cusp_colors[n],
                               troughcolor=self.bgcolor, borderwidth=1,
                               relief=Tk_.FLAT,
                               variable=Tk_.DoubleVar(self))
            slider.index = n
            slider.stamp = 0
            slider.bind('<ButtonPress-1>', self.start_radius)
            slider.bind('<ButtonRelease-1>', self.end_radius)
            slider.grid(padx=(0,20), pady=0, sticky=Tk_.W)
            self.cusp_sliders.append(slider)
            volume_label = ttk.Label(self.slider_frame, width=6)
            volume_label.grid(row=n+1, column=3, sticky=Tk_.W)
            self.volume_labels.append(volume_label)

    def new_scene (self, new_nbhd):
        self.nbhd = new_nbhd
        self.empty = (self.nbhd is None)
        self.set_ties()
        if new_nbhd and self.which_cusp >= new_nbhd.num_cusps():
            self.which_cusp = 0
        while self.volume_labels:
            label = self.volume_labels.pop()
            label.grid_forget()
            label.destroy()
        while self.cusp_sliders:
            slider = self.cusp_sliders.pop()
            slider.destroy()
        while self.slider_frames:
            frame = self.slider_frames.pop()
            frame.grid_forget()
            frame.destroy()
        while self.tie_buttons:
            button = self.tie_buttons.pop()
            button.grid_forget()
            button.destroy()
        while self.eye_buttons:
            button = self.eye_buttons.pop()
            button.grid_forget()
            button.destroy()
        self.eye_var.set(self.which_cusp)
        self.build_sliders()
        self.widget.tk.call(self.widget._w, 'makecurrent')
        self.scene = HoroballScene(new_nbhd, self.pgram_var,
            self.Ford_var, self.tri_var, self.horo_var, self.label_var,
            flipped=self.flip_var.get(), cutoff=self.cutoff,
            which_cusp=self.which_cusp, togl_widget=self.widget)
        assert(self.scene is not None)
        self.widget.redraw_impl = self.scene.draw
        self.configure_sliders()
        self.rebuild()

    def click(self, event):
        self.mouse_x = event.x
        self.mouse_y = event.y
        # Make sure that the scale is reasonable before dragging.
        self.set_zoom(self.zoom.get())

    def flip(self):
        flipped = self.flip_var.get()
        self.scene.flip(flipped)
        self.widget.flipped = flipped
        self.widget.redraw_if_initialized()

    def configure_sliders(self):
        nbhd = self.nbhd
        if self.nbhd is None:
            return
        slider_width = 30
        size = 330 - slider_width
        max_reach = nbhd.max_reach()
        for n in range(nbhd.num_cusps()):
            stopper_color = self.cusp_colors[nbhd.stopper(n)]
            stop = float(nbhd.stopping_displacement(n))
            length = int(stop*size/max_reach) + slider_width
            disp = float(nbhd.get_displacement(n))
            position = 100.0*disp/stop
            # print stop, length, disp position
            self.cusp_sliders[n].set(position)
            self.slider_frames[n].config(background=stopper_color)
            self.volume_labels[n].config(text='%.4f'%nbhd.volume(n))
            self.cusp_sliders[n].config(length=length,
                                        command=self.update_radius)
        self.update_idletasks()

    def translate(self, event):
        """
        Translate the HoroballScene.
        """
        X = self.scale*(event.x - self.mouse_x)
        Y = self.scale*(self.mouse_y - event.y)
        self.mouse_x, self.mouse_y = event.x, event.y
        self.scene.translate(X + Y*1j)
        self.widget.tkTranslate(event)

    # Subclasses may override this, e.g. if they use a help menu.
    def add_help(self):
        help = Button(self.top_frame, text = 'Help', width = 4,
                      borderwidth=0, highlightthickness=0,
                      background=self.bgcolor, command = self.widget.help)
        help.grid(row=0, column=5, sticky=E, pady=3)
        self.top_frame.columnconfigure(5, weight=1)

  # Subclasses may override this to provide menus.
    def build_menus(self):
        pass

  # Subclasses may override this to update menus, e.g. when embedded in a larger window.
    def update_menus(self, menubar):
        pass

    def close(self, event=None):
        self.destroy()

    def redraw(self):
        self.widget.redraw_if_initialized()

    def set_zoom(self, x):
        fovy = 1.0 + (100.0-float(x))/15.0
        self.widget.fovy = fovy
        height = self.widget.winfo_height()
        if height > 0:
            self.scale = fovy / height
        else:
            self.update_idletasks()
            self.after(50, self.set_zoom, x)
        self.widget.redraw_if_initialized()

    def rebuild(self, full_list=True):
        self.set_ties()
        self.configure_sliders()
        self.widget.make_current()
        self.scene.build_scene(which_cusp=self.which_cusp, full_list=full_list)
        self.widget.redraw_if_initialized()

    def start_radius(self, event):
        self.cusp_moving = True
        self.moving_cusp = index = event.widget.index
        self.last_slider_value = self.cusp_sliders[index].get()
        self.update_radius()

    def update_radius(self, event=None):
        index = self.moving_cusp
        value = self.cusp_sliders[index].get()
        if value == self.last_slider_value:
            return
        if self.busy_drawing:
            return
        self.last_slider_value = value
        stop = float(self.nbhd.stopping_displacement(index))
        disp = value*stop/100.0
        self.nbhd.set_displacement(disp, index)
        self.busy_drawing = True
        self.rebuild(full_list=False)
        self.busy_drawing = False

    def end_radius(self, event):
        self.cusp_moving = False
        self.rebuild()

    def set_eye(self):
        self.which_cusp = self.eye_var.get()
        self.rebuild()

    def set_ties(self):
        if self.nbhd == None:
            return
        if len(self.tie_vars) == self.nbhd.num_cusps():
            for n, var in enumerate(self.tie_vars):
                self.nbhd.set_tie(n, var.get())

    def set_cutoff(self, event=None):
        try:
            self.cutoff = float(self.cutoff_var.get())
            self.scene.set_cutoff(self.cutoff)
            self.rebuild()
        except:
            pass
        self.cutoff_var.set('%.4f'%self.cutoff)

    def delete_resource(self):
        try:
            self.scene.delete_resource()
        except AttributeError:
            pass

    def test(self):
        X = 100
        self.widget.event_generate('<Button-1>', x=X, y=300, warp=True)
        self.update_idletasks()
        for n in range(10):
            X += 30
            time.sleep(0.1)
            self.widget.event_generate('<B1-Motion>', x=X, y=300, warp=True)
        self.widget.event_generate('<ButtonRelease-1>', x=X+30, y=300, warp=True)
        self.update_idletasks()
        time.sleep(0.5)
        self.label_var.set(0)
        self.update_idletasks()
        time.sleep(0.5)
        self.cusp_sliders[0].set(50)
        self.update_idletasks()
        time.sleep(1.0)
        self.set_zoom(90)
        self.update_idletasks()
        time.sleep(0.5)


__doc__ = """
   The horoviewer module exports the HoroballViewer class, which is
   a Tkinter / OpenGL window for viewing cusp neighborhoods.
   """

__all__ = ['HoroballViewer']

if __name__ == '__main__':
    import snappy
    from snappy.gui import ViewerWindow
    if len(sys.argv) > 1:
        mfld = sys.argv[1]
    else:
        mfld = 'm125'
    M = snappy.Manifold(mfld)
    HV = ViewerWindow(HoroballViewer, M.cusp_neighborhood())
    HV.mainloop()
