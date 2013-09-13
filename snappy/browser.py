# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import sys, os
try:
    import Tkinter as Tk_
    import ttk
    from SimpleDialog import SimpleDialog
except ImportError:
    import tkinter as Tk_
    from tkinter import ttk
    from tkinter.simpledialog import SimpleDialog
from snappy.polyviewer import PolyhedronViewer
from snappy.horoviewer import HoroballViewer, GetColor
from snappy.app_menus import dirichlet_menus, horoball_menus, browser_menus
from snappy.app_menus import togl_save_image
from snappy.SnapPy import SnapPeaFatalError

# The ttk.LabelFrame is designed to go in a standard window.
# If placed in a ttk.Notebook it will have the wrong background
# color, since the notebook has a darker background than a
# standard window.  This hack fixes that, by overlaying a label with
# the correct background.

ttk_style = None
window_master = None

def init_style():
    """
    Initialize ttk style attributes.  Must be called after creating
    a Tk root window.
    """
    global ST_args, SM_args, GroupBG, WindowBG
    ttk_style = ttk.Style()
    if sys.platform == 'darwin':
        WindowBG = 'SystemDialogBackgroundActive'
        GroupBG = '#e0e0e0'
    else:
        WindowBG = GroupBG = ttk_style.lookup('TLabelframe', 'background')
    ST_args = {
        'selectborderwidth' : 0,
        'highlightbackground' : WindowBG,
        'highlightcolor' : WindowBG,
        'readonlybackground' : WindowBG,
        'relief' : Tk_.FLAT,
        'takefocus' : False,
        'state' : 'readonly'}

    SM_args = {
        'background' : WindowBG,
        'selectborderwidth' : 0,
        'highlightbackground' : WindowBG,
        'highlightcolor' : WindowBG,
        'takefocus': False,
        'state': Tk_.DISABLED,
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
        self.value = value = Tk_.Entry(self, textvariable=self.var, **ST_args)
        self.value.pack(padx=2, pady=2)
        
    def set(self, value):
        self.var.set(value)
        self.value.selection_clear()

    def get(self):
        return self.var.get()

class AutoScrollbar(Tk_.Scrollbar):
    # Frederick Lundh's scrollbar that hides itself if it's not needed.
    def set(self, lo, hi):
        if float(lo) <= 0.0 and float(hi) >= 1.0:
            self.grid_remove()
        else:
            self.grid()
        Tk_.Scrollbar.set(self, lo, hi)

class SelectableMessage(NBLabelframe):
    def __init__(self, master, labeltext=''):
        NBLabelframe.__init__(self, master, text=labeltext)
        self.scrollbar = AutoScrollbar(self, orient=Tk_.VERTICAL)
        self.value = Tk_.Text(self, width=40, height=10,
                              yscrollcommand=self.scrollbar.set,
                              **SM_args)
        self.value.bind('<KeyPress>', lambda event: 'break')
        self.value.bind('<<Paste>>', lambda event: 'break')
        self.value.bind('<<Copy>>', self.copy)
        self.scrollbar.config(command=self.value.yview)
        self.grid_columnconfigure(0, weight=1)
        self.value.grid(row=0, column=0, sticky=Tk_.NSEW)
        self.scrollbar.grid(row=0, column=1, sticky=Tk_.NS)

    def set(self, value):
        self.value.config(state=Tk_.NORMAL)
        self.value.delete('0.1', Tk_.END)
        self.value.insert(Tk_.INSERT, value)
        self.value.config(state=Tk_.DISABLED)
        self.value.selection_clear()

    def get(self):
        return self.value.get('0.1', Tk_.END)

    def copy(self, event):
        self.value.selection_get(selection='CLIPBOARD')

class DirichletTab(PolyhedronViewer):
    def __init__(self, facedicts, root=None, title='Polyhedron Tab',
                 container=None):
        self.focus_var = Tk_.IntVar()
        self.window_master = window_master
        PolyhedronViewer.__init__(self, facedicts, root=root,
                                  title=title, container=container,
                                  bgcolor=WindowBG)
    def add_help(self):
        pass
    
    build_menus = dirichlet_menus

    save_image = togl_save_image

    def close(self):
        pass

class CuspNeighborhoodTab(HoroballViewer):
    def __init__(self, nbhd, root=None, title='Polyhedron Tab',
                 container=None):
        self.focus_var = Tk_.IntVar()
        self.window_master = window_master
        HoroballViewer.__init__(self, nbhd, root=root,
                                title=title, container=container,
                                bgcolor=WindowBG)
    def add_help(self):
        pass

    def view_check(self):
        if self.horo_var.get():
            self.widget.set_background(0.3, 0.3, 0.4)
        else:
            self.widget.set_background(1.0, 1.0, 1.0)
        self.widget.tkRedraw()

    build_menus = horoball_menus

    save_image = togl_save_image

    def close(self):
        pass

class Browser:
    def __init__(self, manifold, root=None):
        if root is None:
            if Tk_._default_root is None:
                root = Tk_.Tk(className='snappy')
                root.iconify()
            else:
                root = Tk_._default_root
        self.root = root
        self.manifold = manifold
        self.symmetry_group = None
        self.length_spectrum = []
        self.window = window = Tk_.Toplevel(root, class_='snappy')
        window.title(manifold.name())
        if ttk_style == None:
            init_style()
        window.config(bg=GroupBG)
        window.protocol("WM_DELETE_WINDOW", self.close)
        self.window_master = window_master
        self.notebook = notebook = ttk.Notebook(window)
        self.notebook.bind('<<NotebookTabChanged>>', self.update_current_tab)
        self.build_invariants()
        self.dirichlet_frame = Tk_.Frame(window)
        self.dirichlet_viewer = DirichletTab(
            facedicts=[],
            root=window,
            container=self.dirichlet_frame)
        notebook.add(self.dirichlet_frame, text='Dirichlet')
        self.update_dirichlet()
        self.horoball_frame = Tk_.Frame(window)
        self.horoball_viewer = CuspNeighborhoodTab(
            nbhd=None,
            root=window,
            container=self.horoball_frame)
        notebook.add(self.horoball_frame, text='Cusp Nbhds')
        window.grid_columnconfigure(1, weight=1)
        window.grid_rowconfigure(0, weight=1)
        self.side_panel = self.build_side_panel()
        self.bottombar = Tk_.Frame(window, height=20, bg='white')
        self.modeline = Tk_.Text(self.bottombar, height=1,
                                 relief=Tk_.FLAT,
                                 background='white',
                                 selectborderwidth=0,
                                 highlightbackground='white',
                                 highlightcolor='white',
                                 highlightthickness=0,
                                 takefocus=False,
                                 state=Tk_.DISABLED)
        self.modeline.tag_config('alert', foreground='red')
        self.modeline.pack(fill=Tk_.BOTH, expand=True, padx=30)
        self.side_panel.grid(row=0, column=0, sticky=Tk_.NSEW, padx=0, pady=0)
        notebook.grid(row=0, column=1, sticky=Tk_.NSEW, padx=0, pady=0)
        self.bottombar.grid(row=1, columnspan=2, sticky=Tk_.NSEW)
        self.build_menus()
        self.window.config(menu=self.menubar)
        self.update_invariants()

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
            cutoff = float(P)
            if self.length_cutoff != cutoff:
                self.length_cutoff = cutoff
                self.update_length_spectrum()
        except ValueError:
            self.window.after_idle( 
                self.cutoff_var.set, str(self.length_cutoff))
            return False
        return True
    
    def save(self):
        # should save a .tri file
        pass

    build_menus = browser_menus

    def build_invariants(self):
        self.invariant_frame = frame = Tk_.Frame(self.window, bg=GroupBG)
        #frame.columnconfigure(0, weight=1)
        frame.columnconfigure(1, weight=1)
        self.volume = SelectableText(frame, labeltext='Volume')
        self.volume.grid(row=0, column=0, padx=30, pady=5, sticky=Tk_.E)
        self.cs = SelectableText(frame, labeltext='Chern-Simons Invariant')
        self.cs.grid(row=1, column=0, padx=30, pady=5, sticky=Tk_.E)
        self.homology = SelectableText(frame, labeltext='First Homology')
        self.homology.grid(row=2, column=0, padx=30, pady=5, sticky=Tk_.E)
        self.symmetry = SelectableText(frame, labeltext='Symmetry Group')
        self.symmetry.grid(row=3, column=0, padx=30, pady=5, sticky=Tk_.E)
        self.pi_one = SelectableMessage(frame, labeltext='Fundamental Group')
        self.pi_one.grid(row=0, column=1, rowspan=3,
                         padx=30, pady=5, sticky=Tk_.NSEW)
        self.pi_one_options = Tk_.Frame(frame, bg=GroupBG)
        self.simplify_var = Tk_.BooleanVar(frame, value=True)
        self.simplify = Tk_.Checkbutton(
            self.pi_one_options,
            variable=self.simplify_var,
            text='simplified presentation',
            bg=GroupBG, borderwidth=0, highlightthickness=0,
            command=self.compute_pi_one)
        self.simplify.pack(anchor=Tk_.W)
        self.minimize_var = Tk_.BooleanVar(frame, value=True)
        self.minimize = Tk_.Checkbutton(
            self.pi_one_options,
            variable=self.minimize_var,
            text='minimal number of generators',
            bg=GroupBG, borderwidth=0, highlightthickness=0,
            command=self.compute_pi_one)
        self.minimize.pack(anchor=Tk_.W)
        self.gens_change_var = Tk_.BooleanVar(frame, value=True)
        self.gens_change = Tk_.Checkbutton(
            self.pi_one_options,
            variable=self.gens_change_var,
            text='fillings may affect generators',
            bg=GroupBG, borderwidth=0, highlightthickness=0,
            command=self.compute_pi_one)
        self.gens_change.pack(anchor=Tk_.W)
        self.pi_one_options.grid(row=3, column=1, padx=30, sticky=Tk_.EW)
        self.length_spectrum = NBLabelframe(frame,
                                            text='Length Spectrum')
        self.length_spectrum.grid_columnconfigure(1, weight=1)
        ttk.Label(self.length_spectrum, text='Length Cutoff:'
                  ).grid(row=0, column=0, sticky=Tk_.E, padx=5, pady=5)
        self.length_cutoff = 1.0
        self.cutoff_var=Tk_.StringVar(self.window, self.length_cutoff)
        self.cutoff_entry = cutoff_entry = ttk.Entry(
            self.length_spectrum,
            takefocus=False,
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
        self.side_panel.grid_rowconfigure(5, weight=1)
        filling = ttk.Labelframe(side_panel, text='Dehn Filling')
        self.filling_vars=[]
        for n in range(self.manifold.num_cusps()):
            R, G, B, A = GetColor(n)
            color = '#%.3x%.3x%.3x'%(int(R*4095), int(G*4095), int(B*4095))
            cusp = ttk.Labelframe(
                filling,
                labelwidget=ttk.Label(filling,
                                      foreground=color,
                                      text='Cusp %d'%n))
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
        filling.grid(row=0, column=0, sticky=Tk_.N, pady=10, padx=5)
        ttk.Button(side_panel, text='Drill ...',
                   command=self.drill).grid(
                       row=1, column=0, padx=10, pady=10, sticky=Tk_.EW)
        ttk.Button(side_panel, text='Cover ...',
                   command=self.cover).grid(
                       row=2, column=0, padx=10, pady=10, sticky=Tk_.EW)
        ttk.Button(side_panel, text='Identify ...',
                   command=self.cover).grid(
                       row=3, column=0, padx=10, pady=10, sticky=Tk_.EW)
        ttk.Button(side_panel, text='Retriangulate',
                   command=self.retriangulate).grid(
                       row=4, column=0, padx=10, pady=10, sticky=Tk_.EW)
        return side_panel

    def do_filling(self):
        filling_spec = [( float(x[0].get()), float(x[1].get()) )
                         for x in self.filling_vars]
        self.window.config(cursor='watch')
        self.clear_invariants()
        self.window.update_idletasks()
        self.manifold.dehn_fill(filling_spec)
        current_fillings = [c.filling for c in self.manifold.cusp_info()]
        for n, coeffs in enumerate(current_fillings):
            for m in (0,1):
                self.filling_vars[n][m].set('%g'%coeffs[m])
        self.update_current_tab()
        self.window.config(cursor='')

    def drill(self):
        dialog = Driller(self.window, self.manifold)
        dialog.go()
        for n in dialog.result:
            self.manifold.drill(n).browse()

    def cover(self):
        pass

    def retriangulate(self):
        self.manifold.randomize()
        self.update_current_tab()

    def update_modeline(self):
        modeline = self.modeline
        modeline.config(state=Tk_.NORMAL)
        modeline.delete(1.0, Tk_.END)
        if not self.manifold.is_orientable():
            modeline.insert(Tk_.END, 'Non-orientable; ')
        modeline.insert(Tk_.END, '%s tetrahedra; %s'%(
            self.manifold.num_tetrahedra(),
            self.manifold.solution_type())
                        )
        if len(self.dirichlet) == 0:
             modeline.insert(Tk_.END,
                             '  * failed to compute Dirichlet domain!',
                             'alert')
        modeline.config(state=Tk_.DISABLED)

    def update_current_tab(self, event=None):
        self.window.update_idletasks()
        self.update_dirichlet()
        self.update_modeline()
        self.update_panel()
        tab_name = self.notebook.tab(self.notebook.select(), 'text')
        if tab_name == 'Invariants':
            self.window.config(menu=self.menubar)
            self.update_invariants()
        if tab_name == 'Cusp Nbhds':
            self.window.config(menu=self.horoball_viewer.menubar)
            self.update_cusps()
        elif tab_name == 'Dirichlet':
            self.window.config(menu=self.dirichlet_viewer.menubar)
            self.dirichlet_viewer.new_polyhedron(self.dirichlet)
        self.window.update_idletasks()

    def update_panel(self):
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
        try:
            self.symmetry_group = self.manifold.symmetry_group()
        except (ValueError, SnapPeaFatalError):
            self.symmetry_group = str('unknown')
        self.symmetry.set(str(self.symmetry_group))
        self.homology.set(repr(self.manifold.homology()))
        self.compute_pi_one()
        self.update_length_spectrum()

    def clear_invariants(self):
        self.volume.set('')
        self.cs.set('')
        self.symmetry.set('')
        self.homology.set('')
        self.pi_one.set('')
        self.geodesics.delete(*self.geodesics.get_children())

    def update_length_spectrum(self):
        try:
            self.length_spectrum = self.manifold.length_spectrum(
                self.length_cutoff)
        except RuntimeError:
            self.length_spectrum = []
        self.geodesics.delete(*self.geodesics.get_children())
        for geodesic in self.length_spectrum:
            parity = '+' if geodesic['parity'].endswith('preserving') else '-'
            self.geodesics.insert('', 'end', values=(
                    geodesic['multiplicity'],
                    geodesic['length'],
                    geodesic['topology'],
                    parity))
        self.cutoff_entry.selection_clear()
        
    def update_dirichlet(self):
        try:
            self.dirichlet = self.manifold.dirichlet_domain().face_list()
        except RuntimeError:
            self.dirichlet = []

    def update_cusps(self):
        try:
            nbhd = self.manifold.cusp_neighborhood()
        except RuntimeError:
            nbhd = None
        self.horoball_viewer.new_scene(nbhd)
        self.window.after(100,
                          self.horoball_viewer.cutoff_entry.selection_clear)
        self.window.focus_set()
        
    def compute_pi_one(self):
        fun_gp = self.manifold.fundamental_group(
            simplify_presentation=self.simplify_var.get(),
            minimize_number_of_generators=self.minimize_var.get(),
            fillings_may_affect_generators=self.gens_change_var.get())
        self.pi_one.set(repr(fun_gp))
        
    def close(self):
        self.window.destroy()

class Driller(SimpleDialog):
    def __init__(self, master, manifold):
        self.manifold = manifold
        self.num = self.result = None
        self.root = root = Tk_.Toplevel(master, class_='SnapPy')
        title = 'Drill'
        root.title(title)
        root.iconname(title)
        root.bind('<Return>', self.return_event)
        frame = Tk_.Frame(self.root)
        frame.grid_columnconfigure(0, weight=1)
        frame.grid_rowconfigure(1, weight=1)
        msg = Tk_.Label(frame, text='Choose curves to drill out:')
        msg.grid(row=0, column=0, pady=6)
        curve_list = [(x['parity'], x['filled_length'])
                      for x in self.manifold.dual_curves()]
        self.curves = curves = ttk.Treeview(
            frame,
            selectmode='extended',
            columns=['parity', 'length'],
            show='headings')
        curves.heading('parity', text='Parity')
        curves.column('parity', stretch=False, width=80)
        curves.heading('length', text='Length')
        curves.column('length', stretch=True)
        for curve in self.manifold.dual_curves():
            parity = '+' if curve['parity'] == 1 else '-'
            length = curve['filled_length']
            curves.insert( '', 'end', values=(parity, length) )
        self.curves.grid(row=1, column=0, padx=6, pady=6, sticky=Tk_.NSEW)
        frame.pack(fill=Tk_.BOTH, expand=1) 
        button_frame = Tk_.Frame(self.root)
        button = ttk.Button(button_frame, text='Drill', command=self.drill,
                            default='active')
        button.pack(side=Tk_.LEFT, padx=6)
        button = ttk.Button(button_frame, text='Cancel', command=self.cancel)
        button.pack(side=Tk_.LEFT, padx=6)
        button_frame.pack()
        self.root.protocol('WM_DELETE_WINDOW', self.wm_delete_window)
        self._set_transient(master)

    def drill(self):
        self.result = [self.curves.index(x) for x in self.curves.selection()]
        self.root.quit()
    
    def cancel(self):
        self.root.quit()

if __name__ == '__main__':
    from snappy import Manifold
    M = Manifold('m125')
    root = Tk_.Tk()
    root.withdraw()
    browser = Browser(M, root)
    root.wait_window(browser.window)
