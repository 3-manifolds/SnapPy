#!/usr/bin/env python

try:
    import Tkinter as Tk_
except ImportError: #Python 3
    import tkinter as Tk_
    
from snappy.CyOpenGL import *
from colorsys import hls_to_rgb

import os, sys

class HoroballViewer:

    def __init__(self, nbhd, which_cusp=0, cutoff=None,
                 root=None, title='Horoball Viewer',
                 prefs={'cusp_horoballs' : True,
                        'cusp_triangulation' : True,
                        'cusp_ford_domain' : True,
                        'cusp_labels' : True,
                        'cusp_parallelogram' : True,
                        'cusp_cutoff' : '0.1000'},
                 container=None):
        self.nbhd = nbhd
        if cutoff is None:
            self.cutoff = float(prefs['cusp_cutoff'])
        else:
            self.cutoff = float(cutoff)
        self.which_cusp = which_cusp
        self.moving_cusp = 0
        self.cusp_moving = False
        self.title = title
        if root is None:
            self.root = Tk_._default_root
        else:
            self.root = root
        if container:
            self.window = window = container
        else:
            self.window = window = Tk_.Toplevel(master=root, class_='snappy')
            window.withdraw()
            window.protocol("WM_DELETE_WINDOW", self.close)
            window.title(title)
        self.pgram_var = pgram_var = Tk_.IntVar(window,
                                                value=prefs['cusp_parallelogram'])
        self.Ford_var = Ford_var = Tk_.IntVar(window,
                                              value=prefs['cusp_ford_domain'])
        self.tri_var = tri_var = Tk_.IntVar(window,
                                            value=prefs['cusp_triangulation'])
        self.horo_var = horo_var = Tk_.IntVar(window,
                                              value=prefs['cusp_horoballs'])
        self.label_var = label_var = Tk_.IntVar(window,
                                                value=prefs['cusp_labels'])
        self.flip_var = flip_var = Tk_.BooleanVar(window)
        window.columnconfigure(0, weight=1)
        window.rowconfigure(1, weight=1)
        self.topframe = topframe = Tk_.Frame(window, borderwidth=0,
                                             relief=Tk_.FLAT)
        self.bottomframe = bottomframe = Tk_.Frame(window, borderwidth=0,
                                             relief=Tk_.FLAT)
        meridian, longitude = nbhd.translations(which_cusp)
        self.widget = widget = OpenGLOrthoWidget(master=bottomframe,
                                            width=600,
                                            height=600,
                                            fovy=3.0,
                                            depth=1,
                                            double=True,
                                            swapinterval=0,
                                            help = """
Use the mouse to drag the scene relative to the
fundamental parallelogram.  

Use the sliders to adjust the sizes of the
horoballs. Color coding indicates who bumps who.

To change the cutoff size, enter a number in
the box and hit return.

Cusps which are "tied" change size in unison.

To view the scene from outside of the upper
half-space, check the the "Flip" checkbutton.

The View menu controls which components of the
scene are visible.
""")
        self.scale = 3.0/600
        widget.bind('<ButtonPress-1>', self.click)
        widget.bind('<B1-Motion>', self.translate)
        widget.set_background(0.3, 0.3, 0.4)
        widget.autospin_allowed = 0
        self.GL = GL_context()
        self.GLU = GLU_context()
        flip_button = Tk_.Checkbutton(topframe, text='Flip',
                                      variable = self.flip_var,
                                      command = self.flip,
                                      highlightthickness=0)
        flip_button.grid(row=0, column=0, sticky=Tk_.E, padx=0, pady=0)
        Tk_.Label(topframe, text='Cutoff').grid(row=1, column=0, sticky=Tk_.E)
        self.cutoff_var = cutoff_var = Tk_.StringVar(window,
                                                     value='%.4f'%self.cutoff)
        cutoff_entry = Tk_.Entry(topframe,
                                 width=6,
                                 textvariable=cutoff_var)
        cutoff_entry.bind('<Return>', self.set_cutoff)
        cutoff_entry.grid(row=1, column=1, sticky=Tk_.W, padx=(0,20))
        Tk_.Label(topframe, text='Tie').grid(row=0, column=2,
                                             sticky=Tk_.W, pady=0)
        Tk_.Label(topframe, text='Cusp radius').grid(row=0, column=3, pady=0)
        self.cusp_sliders = []
        self.slider_frames = []
        self.tie_buttons = []
        self.build_sliders()
        topframe.grid_columnconfigure(3, weight=1)
        topframe.grid(row=0, column=0, sticky=Tk_.NSEW, padx=6, pady=3)
        zoomframe = Tk_.Frame(bottomframe, borderwidth=0, relief=Tk_.FLAT)
        self.zoom = zoom = Tk_.Scale(zoomframe, showvalue=0, from_=100, to=0,
                                     command = self.set_zoom, width=11,
                                     troughcolor='#f4f4f4', borderwidth=1,
                                     relief=Tk_.SUNKEN)
        zoom.set(30)
        spacer = Tk_.Frame(zoomframe, height=14, borderwidth=0, relief=Tk_.FLAT)
        zoom.pack(side=Tk_.TOP, expand=Tk_.YES, fill=Tk_.Y)
        spacer.pack()
        bottomframe.columnconfigure(0, weight=1)
        widget.grid(row=0, column=0, sticky=Tk_.EW)
        zoomframe.grid(row=0, column=1, sticky=Tk_.NS)
        bottomframe.grid(row=1, column=0, sticky=Tk_.NSEW)
#        self.configure_sliders(390)
        self.build_menus()
        self.mouse_x = 0
        self.mouse_y = 0
        self.movie_id=0
        if container is None:
            window.deiconify()
            window.update()  # Seems to avoid a race condition with togl
        window.bind('<Configure>', self.handle_resize)
        bottomframe.bind('<Configure>', self.togl_handle_resize)
        self.scene = HoroballScene(nbhd, pgram_var, Ford_var, tri_var,
                                   horo_var, label_var,
                                   flipped=self.flip_var.get(),
                                   cutoff=self.cutoff,
                                   which_cusp=self.which_cusp)
        self.widget.redraw = self.scene.draw
        window.update_idletasks()
        self.configure_sliders()

    def build_sliders(self):
        self.cusp_vars = []
        self.cusp_colors = []
        self.tie_vars = []
        self.tie_dict = {}

        for n in range(self.nbhd.num_cusps()):
            disp = self.nbhd.stopping_displacement(which_cusp=n)
            self.nbhd.set_displacement(disp, which_cusp=n)
            tie_var = Tk_.IntVar(self.window)
            self.tie_vars.append(tie_var)
            self.tie_dict[str(tie_var)] = n
            tie_var.trace('w', self.set_tie)
            tie_button = Tk_.Checkbutton(self.topframe, variable = tie_var)
            tie_button.index = n
            tie_button.grid(row=n+1, column=2, sticky=Tk_.W)
            self.tie_buttons.append(tie_button)
            R, G, B, A = GetColor(self.nbhd.original_index(n))
            self.cusp_colors.append('#%.3x%.3x%.3x'%(
                int(R*4095), int(G*4095), int(B*4095)))
            self.cusp_vars.append(Tk_.IntVar(self.window))
            self.slider_frames.append(
                Tk_.Frame(self.topframe, borderwidth=1, relief=Tk_.SUNKEN))
            self.slider_frames[n].grid(row=n+1, column=3,
                                       sticky=Tk_.EW, padx=6)
            slider = Tk_.Scale(self.slider_frames[n], 
                               showvalue=0, from_=-0, to=100,
                               width=11, length=200, orient=Tk_.HORIZONTAL,
                               background=self.cusp_colors[n],
                               borderwidth=0, relief=Tk_.FLAT,
                               variable=Tk_.DoubleVar(self.window))
            slider.index = n
            slider.stamp = 0
            slider.bind('<ButtonPress-1>', self.start_radius)
            slider.bind('<ButtonRelease-1>', self.end_radius)
            slider.pack(padx=0, pady=0, side=Tk_.LEFT)
            self.cusp_sliders.append(slider)

    def new_scene(self, new_nbhd):
        self.nbhd = new_nbhd
        if self.which_cusp > new_nbhd.num_cusps():
            self.which_cusp = 0
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
        self.build_sliders()
        self.scene = HoroballScene(new_nbhd, self.pgram_var,
                                   self.Ford_var, self.tri_var,
                                   self.horo_var, self.label_var,
                                   flipped=self.flip_var.get(),
                                   cutoff=self.cutoff,
                                   which_cusp=self.which_cusp)
        self.widget.redraw = self.scene.draw
        self.configure_sliders()
        self.rebuild(full_list=True)
         
    def click(self, event):
        self.mouse_x = event.x
        self.mouse_y = event.y

    def flip(self):
        flipped = self.flip_var.get()
        self.scene.flip(flipped)
        self.widget.flipped = flipped
        self.widget.tkRedraw()

    def handle_resize(self, event):
        self.window.update_idletasks()
        self.configure_sliders()
        
    def configure_sliders(self):#, size=0):
        # The frame width is not valid until the window has been rendered.
        # Supply the expected size if calling from __init__.
#        if size == 0:
#            size = float(self.slider_frames[0].winfo_width() - 10)
        size = float(self.slider_frames[0].winfo_width() - 10)
        max = self.nbhd.max_reach()
        for n in range(self.nbhd.num_cusps()):
            stopper_color = self.cusp_colors[self.nbhd.stopper(n)]
            stop = self.nbhd.stopping_displacement(which_cusp=n)
            disp = self.nbhd.get_displacement(which_cusp=n)
            length = int(stop*size/max)
            self.cusp_sliders[n].config(length=length)
            self.cusp_sliders[n].set(100.0*disp/stop)
            self.slider_frames[n].config(background=stopper_color)
            self.window.update_idletasks()
        self.widget.tkRedraw()

    def togl_handle_resize(self, event):
        self.widget.config(height=self.bottomframe.winfo_height())
        self.widget.tkRedraw()

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
        help = Button(self.topframe, text = 'Help', width = 4,
                      borderwidth=0, highlightthickness=0,
                      background="#f4f4f4", command = self.widget.help)
        help.grid(row=0, column=4, sticky=E, pady=3)
        self.topframe.columnconfigure(3, weight=1)

  # Subclasses may override this to provide menus.
    def build_menus(self):
        pass

    def close(self):
        self.widget.activate()
        self.scene.destroy()
        self.window.destroy()

    def set_zoom(self, x):
        fovy = 1.0 + float(x)/15.0
        self.widget.fovy = fovy
        self.scale = fovy/self.widget.winfo_height()
        self.widget.tkRedraw()

    def rebuild(self, full_list=True):
        self.configure_sliders()
        self.widget.activate()
        self.scene.build_scene(full_list)
        self.widget.tkRedraw()

    def start_radius(self, event):
        self.cusp_moving = True
        self.moving_cusp = event.widget.index
        self.update_radius()

    def update_radius(self):
        index = self.moving_cusp
        value = self.cusp_sliders[index].get()
        stop = self.nbhd.stopping_displacement(index)
        disp = value*stop/100.0
        self.nbhd.set_displacement(disp, index)
        self.rebuild(full_list=False)
        if self.cusp_moving:
            self.movie_id = self.window.after(100, self.update_radius)

    def end_radius(self, event):
        try:
            self.window.after_cancel(self.movie_id)
        except:
            pass
        self.cusp_moving = False
        self.rebuild(full_list=True)
        
    def set_tie(self, name, *args):
        index = self.tie_dict[name]
        value = self.tie_vars[index].get()
        self.nbhd.set_tie(index, value)
        self.rebuild()

    def set_cutoff(self, event):
        try:
            self.cutoff = float(self.cutoff_var.get())
            self.scene.set_cutoff(self.cutoff)
            self.rebuild()
        except:
            pass
        self.cutoff_var.set('%.4f'%self.cutoff)
        
        
__doc__ = """
   The horoviewer module exports the HoroballViewer class, which is
   a Tkinter / OpenGL window for viewing cusp neighborhoods.
   """

__all__ = ['HoroballViewer']

if __name__ == '__main__':
    import snappy
    if len(sys.argv) > 1:
        mfld = sys.argv[1]
    else:
        mfld = 'm125'
    M = snappy.Manifold(mfld)
    HV = HoroballViewer(M.cusp_neighborhood())
    HV.window.mainloop()


