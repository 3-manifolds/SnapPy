# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import sys, os
try:
    import Tkinter as Tk_
    import ttk
except ImportError:
    import tkinter as Tk_
    from tkinter import ttk
from snappy.polyviewer import PolyhedronViewer
from snappy.horoviewer import HoroballViewer

# The ttk.LabelFrame is designed to go in a standard window.
# If placed in a ttk.Notebook it will have the wrong background
# color, since the notebook has a darker background than a
# standard window.  This hack fixes that, by overlaying a label with
# the correct background.

GroupBG = WindowBG = '#000000'
ST_args = SM_args = {}

def init_style():
    """
    Initialize ttk style attributes.  Must be called after creating
    a Tk root window.
    """
    global ST_args, SM_args, GroupBG, WindowBG
    if sys.platform == 'darwin':
        WindowBG = 'SystemDialogBackgroundActive'
        GroupBG = ttk.Style().lookup('TLabelframe', 'background')
    else:
        WindowBG = GroupBG = ttk.Style().lookup('TLabelframe', 'background')
    ST_args = {
        'selectborderwidth' : 0,
        'highlightbackground' : WindowBG,
        'highlightcolor' : WindowBG,
        'readonlybackground' : WindowBG,
        'relief' : Tk_.FLAT,
        'state' : 'readonly'}
    SM_args = {
        'background' : WindowBG,
        'selectborderwidth' : 0,
        'highlightbackground' : WindowBG,
        'highlightcolor' : WindowBG,
        'relief' : Tk_.FLAT
        }


class NBLabelframeMac(ttk.Labelframe):
    def __init__(self, master, text=''):
        ttk.Labelframe.__init__(self, master, text=' ')
        self.overlay = Tk_.Label(self, text=text,
                bg=GroupBG,
                padx=12,
                anchor=Tk_.W,
                relief=Tk_.FLAT,
                borderwidth=1,
                highlightbackground=GroupBG,
                highlightcolor=GroupBG)
        self.overlay.place(relwidth=1, x=0, y=-2, bordermode="outside")

if sys.platform == 'darwin':
    NBLabelframe = NBLabelframeMac
else:
    NBLabelframe = ttk.Labelframe

class SelectableText(NBLabelframe):
    def __init__(self, master, labeltext=''):
        NBLabelframe.__init__(self, master, text=labeltext)
        self.var = Tk_.StringVar(master)
        self.value = Tk_.Entry(self, textvariable=self.var, **ST_args)
        self.value.pack(padx=2, pady=2)
        
    def set(self, value):
        self.var.set(value)
        self.value.selection_clear()

    def get(self):
        return self.var.get()

class SelectableMessage(NBLabelframe):
    def __init__(self, master, labeltext=''):
        NBLabelframe.__init__(self, master, text=labeltext)
        self.var = Tk_.StringVar(master)
        self.value = Tk_.Text(self, width=30, height=10,
                              **SM_args)
        self.value.bind('<KeyPress>', lambda event: 'break')
        self.value.bind('<<Paste>>', lambda event: 'break')
        self.value.bind('<<Copy>>', self.copy)
        self.value.pack(padx=2, pady=2)

    def set(self, value):
        self.value.delete('0.1', Tk_.END)
        self.value.insert(Tk_.INSERT, value)

    def get(self):
        return self.value.get('0.1', Tk_.END)

    def copy(self, event):
        self.value.selection_get(selection='CLIPBOARD')

class DirichletTab(PolyhedronViewer):
    def __init__(self, facedicts, root=None, title='Polyhedron Tab',
                 container=None):
        self.focus_var = Tk_.IntVar()
        self.window_master = None
        PolyhedronViewer.__init__(self, facedicts, root=root,
                                  title=title, container=container)
    def add_help(self):
        pass

    def build_menus(self):
        pass

    def close(self):
        pass

class CuspNeighborhoodTab(HoroballViewer):
    def __init__(self, nbhd, root=None, title='Polyhedron Tab',
                 container=None):
        self.focus_var = Tk_.IntVar()
        self.window_master = None
        HoroballViewer.__init__(self, nbhd, root=root,
                                  title=title, container=container)
    def add_help(self):
        pass

    def build_menus(self):
        pass

    def close(self):
        pass

class Browser:
    def __init__(self, manifold, master=None):
        self.manifold = manifold
        if master == None:
            master = Tk_._default_root
        self.window = window = Tk_.Toplevel(master)
        window.title(manifold.name())
        window.config(bg=GroupBG)
        window.protocol("WM_DELETE_WINDOW", self.close)
        if sys.platform == 'darwin':
            this_dir =  os.path.dirname(__file__)
            Tk_path = os.path.join(this_dir, 'togl',
                                   'darwin-tk' + str(Tk_.TkVersion))
#            master.tk.call('lappend', 'auto_path', Tk_path)
#            master.tk.call('package', 'require', 'mactoolbar')
#            window.tk.call('set', 'tk::mac::useCompatibilityMetrics', '0')
#            window.tk.call('mactoolbar::createbutton',
#                'garbage string', 'Hi', "It's SnapPy",
#                os.path.join(this_dir, 'info_icon.gif'),
#                lambda : None)
#            window.tk.call('mactoolbar::create', self.window._w)
            # This must come after creating the toolbar.
#            window.tk.call('tk::unsupported::MacWindowStyle',
#                'style', self.window._w, 'document',
#                ('standardDocument', 'unifiedTitleAndToolbar')
#                )
            #print self.tk.call( 'tk::unsupported::MacWindowStyle',
            #    'style', self.window._w)
        self.notebook = nb = ttk.Notebook(window)
        self.notebook.bind('<<NotebookTabChanged>>', self.update_current_tab)
        self.build_invariants()
        try:
            faces = manifold.dirichlet_domain().face_list()
        except RuntimeError:
            faces = []
        self.dirichlet_frame = Tk_.Frame(window)
        self.dirichlet_viewer = DirichletTab(
            facedicts=faces,
            root=window,
            container=self.dirichlet_frame)
        nb.add(self.dirichlet_frame, text='Dirichlet')
        try:
            nbhd = manifold.cusp_neighborhood()
        except RuntimeError:
            nbhd = None
        self.horoball_frame = Tk_.Frame(window)
        self.horoball_viewer = CuspNeighborhoodTab(
            nbhd=nbhd,
            root=window,
            container=self.horoball_frame)
        nb.add(self.horoball_frame, text='Cusp Nbhds')
        window.grid_columnconfigure(1, weight=1)
        window.grid_rowconfigure(0, weight=1)
        self.side_panel = self.build_side_panel()
        self.status = Tk_.StringVar(window)
        self.bottombar = Tk_.Frame(window, height=20, bg='white')
        bottomlabel = Tk_.Label(self.bottombar, textvar=self.status,
                                anchor=Tk_.W, relief=Tk_.FLAT, bg='white')
        bottomlabel.pack(fill=Tk_.BOTH, expand=True, padx=30)
        self.side_panel.grid(row=0, column=0, sticky=Tk_.NSEW, padx=0, pady=0)
        nb.grid(row=0, column=1, sticky=Tk_.NSEW, padx=0, pady=0)
        self.bottombar.grid(row=1, columnspan=2, sticky=Tk_.NSEW)
        self.update_invariants()
        # temporary ???
        window.geometry('700x600')

    def validate_coeff(self, P, W):
        tkname, cusp, curve = W.split(':')
        cusp, curve = int(cusp), int(curve)
        try:
            float(P)
        except ValueError:
            var = self.filling_vars[cusp][curve]
            var.set('%g'%self.manifold.cusp_info()[cusp].filling[curve])
            return False
        return True

    def validate_cutoff(self, P):
        try:
            self.length_cutoff = float(P)
            self.update_length_spectrum()
        except ValueError:
            self.window.after_idle( 
                self.cutoff_var.set, str(self.length_cutoff))
            return False
        return True
    
    def build_invariants(self):
        self.invariant_frame = frame = Tk_.Frame(self.window, bg=GroupBG)
        frame.columnconfigure(0, weight=1)
        frame.columnconfigure(1, weight=1)
        self.volume = SelectableText(frame, labeltext='Volume')
        self.volume.grid(row=0, column=0, padx=30, pady=5, sticky=Tk_.E)
        self.cs = SelectableText(frame, labeltext='Chern-Simons Invariant')
        self.cs.grid(row=1, column=0, padx=30, pady=5, sticky=Tk_.E)
        self.orblty = SelectableText(frame, labeltext='Orientability')
        self.orblty.grid(row=2, column=0, padx=30, pady=5, sticky=Tk_.E)
        self.homology = SelectableText(frame, labeltext='First Homology')
        self.homology.grid(row=3, column=0, padx=30, pady=5, sticky=Tk_.E)
        self.pi_one = SelectableMessage(frame, labeltext='Fundamental Group')
        self.pi_one.grid(row=0, column=1, rowspan=3,
                         padx=30, pady=5, sticky=Tk_.W+Tk_.N+Tk_.S)
        self.pi_one_options = Tk_.Frame(frame, bg=GroupBG)
        self.simplify_var = Tk_.BooleanVar(frame, value=True)
        self.simplify = Tk_.Checkbutton(
            self.pi_one_options,
            variable=self.simplify_var,
            text='simplified presentation',
            bg=GroupBG,
            command=self.compute_pi_one)
        self.simplify.pack(anchor=Tk_.W)
        self.minimize_var = Tk_.BooleanVar(frame, value=True)
        self.minimize = Tk_.Checkbutton(
            self.pi_one_options,
            variable=self.minimize_var,
            text='minimal number of generators',
            bg=GroupBG,
            command=self.compute_pi_one)
        self.minimize.pack(anchor=Tk_.W)
        self.gens_change_var = Tk_.BooleanVar(frame, value=True)
        self.gens_change = Tk_.Checkbutton(
            self.pi_one_options,
            variable=self.gens_change_var,
            text='fillings may affect generators',
            bg=GroupBG,
            command=self.compute_pi_one)
        self.gens_change.pack(anchor=Tk_.W)
        self.pi_one_options.grid(row=3, column=1, padx=30, sticky=Tk_.W)
        self.length_spectrum = NBLabelframe(frame,
                                            text='Length Spectrum')
        self.length_spectrum.grid_columnconfigure(1, weight=1)
        ttk.Label(self.length_spectrum, text='Length Cutoff:'
                  ).grid(row=0, column=0, sticky=Tk_.E, padx=5, pady=5)
        self.length_cutoff = 1.0
        self.cutoff_var=Tk_.StringVar(self.window, self.length_cutoff)
        cutoff_entry = ttk.Entry(
            self.length_spectrum, 
            width=6,
            textvariable=self.cutoff_var,
            validate='focusout',
            validatecommand=(self.window.register(self.validate_cutoff),'%P')
            )
        cutoff_entry.bind('<Return>', lambda event : self.window.focus_set())
        cutoff_entry.grid(row=0, column=1, sticky=Tk_.W, pady=5)
        self.geodesics = geodesics = ttk.Treeview(
            self.length_spectrum,
            columns=['mult', 'length', 'topology', 'parity'],
            show='headings')
        geodesics.heading('mult', text='Mult.')
        geodesics.column('mult', stretch=False, width=40)
        geodesics.heading('length', text='Length')
        geodesics.column('length', stretch=True)
        geodesics.heading('topology', text='Topology')
        geodesics.column('topology', stretch=False, width=80)
        geodesics.heading('parity', text='Parity')
        geodesics.column('parity', stretch=False, width=80)
        geodesics.grid(row=1, columnspan=2, sticky=Tk_.EW, padx=5, pady=5)
        self.length_spectrum.grid(row=4, columnspan=2, padx=10, pady=10,
                                  sticky=Tk_.EW)
        self.notebook.add(frame, text='Invariants', padding=[0])

    def build_side_panel(self):
        window = self.window
        self.side_panel = side_panel = Tk_.Frame(window, bg=WindowBG)
        self.side_panel.grid_rowconfigure(1, weight=1)
        ttk.Button(side_panel, text='Retriangulate',
                   command=self.retriangulate).grid(
                       row=0, column=0, padx=20, pady=10, sticky=Tk_.EW)

        filling = ttk.Labelframe(side_panel, text='Dehn Filling')
        self.filling_vars=[]
        for n in range(self.manifold.num_cusps()):
            cusp = ttk.Labelframe(filling, text='Cusp %d'%n)
            mer_var = Tk_.StringVar(window)
            long_var = Tk_.StringVar(window)
            self.filling_vars.append((mer_var, long_var))
            ttk.Label(cusp, text='Meridian: ').grid(
                row=0, column=0, sticky=Tk_.E)
            ttk.Label(cusp, text='Longitude: ').grid(
                row=1, column=0, sticky=Tk_.E)
            meridian = ttk.Entry(cusp, width=4,
                textvariable=mer_var,
                name=':%s:0'%n,            
                validate='focusout',
                validatecommand=(window.register(self.validate_coeff),'%P','%W')
                )
            meridian.grid(row=0, column=1, sticky=Tk_.W, padx=3, pady=3)
            longitude = ttk.Entry(cusp, width=4,
                textvariable=long_var,
                name=':%s:1'%n,
                validate='focusout',
                validatecommand=(window.register(self.validate_coeff),'%P','%W')
                )
            longitude.grid(row=1, column=1, sticky=Tk_.W, padx=3, pady=3)
            cusp.grid(row=n, pady=8, padx=5)
        ttk.Button(filling, text='Fill',
                   command=self.do_filling).grid(
                       row=n+1, columnspan=2, padx=20, pady=10, sticky=Tk_.EW)
        filling.grid(row=1, column=0, sticky=Tk_.N, pady=10, padx=5)
        return side_panel

    def do_filling(self):
        filling_spec = [( float(x[0].get()), float(x[1].get()) )
                         for x in self.filling_vars]
        self.manifold.dehn_fill(filling_spec)
        self.status.set('%s tetrahedra; %s'%(
            self.manifold.num_tetrahedra(),
            self.manifold.solution_type())
            )
        current_fillings = [c.filling for c in self.manifold.cusp_info()]
        for n, coeffs in enumerate(current_fillings):
            for m in (0,1):
                self.filling_vars[n][m].set('%g'%coeffs[m])
        self.update_current_tab()

    def retriangulate(self):
        self.manifold.randomize()
        self.update_current_tab()

    def update_current_tab(self, event=None):
        self.window.update_idletasks()
        self.update_panel()
        self.window.update_idletasks()
        tab_name = self.notebook.tab(self.notebook.select(), 'text')
        if tab_name == 'Invariants':
            self.window.after_idle(self.update_invariants)
        if tab_name == 'Cusp Nbhds':
            self.window.after_idle(self.update_cusps)
        elif tab_name == 'Dirichlet':
            self.window.after_idle(self.update_dirichlet)

    def update_panel(self):
        self.status.set('%s tetrahedra; %s'%(
            self.manifold.num_tetrahedra(),
            self.manifold.solution_type())
            )
        current_fillings = [c.filling for c in self.manifold.cusp_info()]
        for n, coeffs in enumerate(current_fillings):
            for m in (0,1):
                self.filling_vars[n][m].set('%g'%coeffs[m])

    def update_invariants(self):
        self.volume.set(repr(self.manifold.volume()))
        try:
            self.cs.set(repr(self.manifold.chern_simons()))
        except ValueError:
            self.cs.set('')
        orblty = ('orientable' if self.manifold.is_orientable()
                  else 'non-orientable')
        self.orblty.set(orblty)
        self.homology.set(repr(self.manifold.homology()))
        self.compute_pi_one()
        self.update_length_spectrum()

    def update_length_spectrum(self):
        spectrum = self.manifold.length_spectrum(self.length_cutoff)
        self.geodesics.delete(*self.geodesics.get_children())
        for geodesic in spectrum:
            parity = '+' if geodesic['parity'].endswith('preserving') else '-'
            self.geodesics.insert('', 'end', values=(
                    geodesic['multiplicity'],
                    geodesic['length'],
                    geodesic['topology'],
                    parity))
        
    def update_dirichlet(self):
        try:
            faces = self.manifold.dirichlet_domain().face_list()
        except RuntimeError:
            faces = []
        self.dirichlet_viewer.new_polyhedron(faces)
        # shouldn't be necessary, but it helps ...
        self.window.after(100, self.dirichlet_viewer.widget.tkRedraw()) 

    def update_cusps(self):
        try:
            nbhd = self.manifold.cusp_neighborhood()
        except RuntimeError:
            nbhd = None
        self.horoball_viewer.new_scene(nbhd)
        self.horoball_viewer.configure_sliders()
        # shouldn't be necessary, but it helps ...
        self.window.after(100, self.horoball_viewer.configure_sliders) 

    def compute_pi_one(self):
        fun_gp = self.manifold.fundamental_group(
            simplify_presentation=self.simplify_var.get(),
            minimize_number_of_generators=self.minimize_var.get(),
            fillings_may_affect_generators=self.gens_change_var.get())
        self.pi_one.set(repr(fun_gp))
        
    def close(self):
        self.window.destroy()
        
if __name__ == '__main__':
    from snappy import *
    M = Manifold('m125')
    root = Tk_.Tk()
    root.withdraw()
    browser = Browser(root, M)
    root.wait_window(browser)
    
    
