# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from builtins import range
import sys, os
from .gui import *
from .polyviewer import PolyhedronViewer
from .horoviewer import HoroballViewer
from .CyOpenGL import GetColor
from .app_menus import browser_menus
from . import app_menus
from .number import Number
from . import database
from .exceptions import SnapPeaFatalError
from plink import LinkViewer, LinkEditor
from plink.ipython_tools import IPythonTkRoot
from spherogram.links.orthogonal import OrthogonalLinkDiagram


main_window = None
#TODO:  Find these dimensions by introspection.
cusp_box_height = 105
if sys.platform == 'darwin':
    cusp_box_width = 180
elif sys.platform == 'win32':
    cusp_box_width = 150
else:
    cusp_box_width = 210


class SelectableText(ttk.Frame):
    """
    A Label and a disabled Entry widget which is disguised as a label
    but, unlike a Label, allows selecting the text.  On the Mac, the
    background color matches the background of a depth 2 LabelFrame
    by default.
    """
    def __init__(self, master, labeltext='', width=18, depth=2):
        ttk.Frame.__init__(self, master)
        self.var = Tk_.StringVar(master)
        style = SnapPyStyle()
        bg_color = style.groupBG if depth == 1 else style.subgroupBG
        self.label = label = ttk.Label(self, text=labeltext)
        self.value = value = Tk_.Entry(self, textvariable=self.var,
            state='readonly', borderwidth=0, readonlybackground=bg_color,
            highlightbackground=bg_color, highlightcolor=bg_color,
            highlightthickness=0, takefocus=False)
        if width:
            value.config(width=width)
        label.pack(side=Tk_.LEFT)
        value.pack(side=Tk_.LEFT, padx=2)

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

class SelectableMessage(ttk.Frame):
    """
    A disabled Text widget which allows selection of text.  On the mac
    the selection does not highlight correctly unless the Text widget has
    focus and does not clear correctly unless the state is NORMAL.
    """
    def __init__(self, master):
        self.master = master
        ttk.Frame.__init__(self, master)
        self.scrollbar = AutoScrollbar(self, orient=Tk_.VERTICAL)
        self.text = text = Tk_.Text(self, width=60, height=12,
            highlightthickness=0, relief=Tk_.FLAT,
            yscrollcommand=self.scrollbar.set)
        self.scrollbar.config(command=text.yview)
        self.grid_columnconfigure(0, weight=1)
        text.grid(row=0, column=0, sticky=Tk_.NSEW)
        text.bind('<<Copy>>', self.copy)
        text.bind('<Button-1>', lambda *args: self.text.focus_set())
        text.bind('<FocusOut>', self.focus_out)
        self.scrollbar.grid(row=0, column=1, sticky=Tk_.NS)
        self.text.config(state=Tk_.DISABLED)

    def focus_out(self, *args):
        self.text.config(state=Tk_.NORMAL)
        self.text.selection_clear()
        self.text.config(state=Tk_.DISABLED)

    def set(self, message):
        self.text.config(state=Tk_.NORMAL)
        self.text.delete('0.1', Tk_.END)
        self.text.selection_clear()
        self.text.insert(Tk_.INSERT, message)
        self.text.config(state=Tk_.DISABLED)

    def get(self):
        return self.text.get('0.1', Tk_.END)

    def copy(self, event):
        self.text.selection_get(selection='CLIPBOARD')

class DirichletTab(PolyhedronViewer):
    def __init__(self, master, facedicts=[], title='Polyhedron Tab', parent=None):
        self.main_window = main_window
        self.style = style = SnapPyStyle()
        PolyhedronViewer.__init__(self, master, facedicts=facedicts,
                                  title=title, bgcolor=style.groupBG)
    def update_menus(self, menubar):
        menubar.children['help'].activate(
            [ app_menus.help_polyhedron_viewer_label,
              app_menus.help_report_bugs_label ])

class CuspNeighborhoodTab(HoroballViewer):
    def __init__(self, nbhd, title='Polyhedron Tab'):
        self.main_window = main_window
        style = SnapPyStyle()
        if main_window:
            HoroballViewer.__init__(self, nbhd, title=title,
                                    bgcolor=style.groupBG, prefs=main_window.prefs)
        else:
            HoroballViewer.__init__(self, nbhd, title=title,
                                    bgcolor=style.groupBG)

    def update_menus(self, menubar):
        menubar.children['help'].activate(
            [ app_menus.help_horoball_viewer_label,
              app_menus.help_report_bugs_label])

class LinkTab(LinkViewer):
    def __init__(self, data, window):
        self.style = style = SnapPyStyle()
        self.canvas = canvas = Tk_.Canvas(window, bg=style.groupBG,
            highlightthickness=0, highlightcolor=style.groupBG)
        LinkViewer.__init__(self, canvas, data)
        canvas.bind("<Configure>", lambda event : self.draw())

    def close(self, event=None):
        pass

class Browser(Tk_.Toplevel):
    """
    Toplevel window displaying a Dehn filling control panel and
    a notebook with tabs that show numerical or graphical information
    about the manifold.
    """

    def __init__(self, manifold, root=None, main_window=None):
        if root is None:
            if Tk_._default_root:
                self.root = root = Tk_._default_root
            else:
                self.root = root = IPythonTkRoot(window_type='Browser')
                root.withdraw()
        else:
            self.root = root
        Tk_.Toplevel.__init__(self, root, class_='snappy')
        self.update_idletasks()
        if isinstance(root, IPythonTkRoot):
            self.withdraw()
            # Avoid showing an empty root window on the screen.
            self.root.after(100, self.deiconify)
        if manifold.num_tetrahedra() == 0:
            raise ValueError('The empty Manifold cannot be browsed.')
        self.manifold = manifold
        self.aka_after_id = None
        self.main_window = main_window
        self.symmetry_group = None
        self.dirichlet = []
        self.cusp_nbhd = None
        self.length_spectrum = []
        self.recompute_invariants = True
        self.style = style = SnapPyStyle()
        self.title(manifold.name())
        self.config(bg=style.groupBG)
        self.protocol("WM_DELETE_WINDOW", self.close)
        if sys.platform == 'darwin':
            self.bind_all('<Command-Key-w>', self.close)
        elif sys.platform == 'linux2' or sys.platform == 'linux':
            self.bind_all('<Alt-Key-F4>', self.close)

        self.side_panel = side_panel = self.build_side_panel()

        self.notebook = notebook = ttk.Notebook(self)
        self.invariants_tab = invariants_tab = self.build_invariants()
        self.dirichlet_viewer = DirichletTab(self)
        self.horoball_viewer = CuspNeighborhoodTab(self)

        self.fillings_changed_callback = None

        self.bottombar = bottombar = ttk.Frame(self, height=20)
        bg = self.style.ttk_style.lookup('TLable', 'background')
        fg = self.style.ttk_style.lookup('TLable', 'foreground')
        self.modeline = modeline = Tk_.Text(
            self.bottombar, height=1,
            relief=Tk_.FLAT,
            background=bg,
            foreground=fg,
            selectborderwidth=0,
            highlightbackground=bg,
            highlightcolor=bg,
            highlightthickness=0,
            takefocus=False,
            state=Tk_.DISABLED)
        self.modeline.tag_config('alert', foreground='red')

        self.inside_view = None

        self.symmetry_tab = symmetry_tab = self.build_symmetry()
        self.link_tab = link_tab = self.build_link()
        notebook.add(invariants_tab, text='Invariants', padding=[0])
        notebook.add(self.dirichlet_viewer, text='Dirichlet')
        notebook.add(self.horoball_viewer, text='Cusp Nbhds')
        notebook.add(self.build_inside_view(), text = 'Inside view')
        notebook.add(symmetry_tab, text='Symmetry', padding=[0])
        if link_tab:
            notebook.add(link_tab.canvas, text='Link')
        notebook.bind('<<NotebookTabChanged>>', self.update_current_tab)

        self.build_menus()
        self.config(menu=self.menubar)
        self.grid_columnconfigure(1, weight=1)
        self.grid_rowconfigure(0, weight=1)
        side_panel.grid(row=0, column=0, sticky=Tk_.NSEW, padx=0, pady=0)
        notebook.grid(row=0, column=1, sticky=Tk_.NSEW, padx=0, pady=0, ipady=5)
        bottombar.grid(row=1, columnspan=2, sticky=Tk_.NSEW)
        self.modeline.pack(fill=Tk_.BOTH, expand=True, padx=30)
        self.update_modeline()
        self.update_idletasks()

    def __repr__(self):
        return 'Browser window for %s\n'%self.manifold

    build_menus = browser_menus

    def build_side_panel(self):
        side_panel = ttk.Frame(self)
        side_panel.grid_rowconfigure(5, weight=1)
        filling = ttk.Labelframe(side_panel, text='Filling Curves',
            padding=(10, 10))
        self.filling_canvas = None
        num_cusps = self.manifold.num_cusps()
        if num_cusps > 3:
            filling.grid_propagate(False)
            filling.configure(width=cusp_box_width + 20,
                              height=int(3.5*cusp_box_height))
            filling.grid_columnconfigure(0, weight=1)
            filling.grid_rowconfigure(0, weight=1)
            filling_scrollbar = ttk.Scrollbar(filling)
            filling_scrollbar.grid(row=0, column=1, sticky=Tk_.NS)
            self.filling_canvas = canvas = Tk_.Canvas(filling, bd=0,
                highlightthickness=0,
                yscrollcommand=filling_scrollbar.set,
                scrollregion=(0, 0, cusp_box_width + 20,
                                  cusp_box_height*num_cusps + 10)
            )
            if sys.platform == 'darwin':
                canvas.configure(background=self.style.groupBG)
            canvas.grid(row=0, column=0, sticky=Tk_.NSEW)
            filling_scrollbar.config(command=canvas.yview)
        self.filling_vars=[]

        # Embedded windows in a canvas are clipped to their parent, not to
        # the canvas.
        cusp_parent = self.filling_canvas if self.filling_canvas else filling
        for n in range(num_cusps):
            R, G, B, A = GetColor(n)
            color = '#%.3x%.3x%.3x'%(int(R*4095), int(G*4095), int(B*4095))
            cusp = ttk.Labelframe(cusp_parent, text='Cusp %d'%n)
            mer_var = Tk_.StringVar(self, value='0')
            long_var = Tk_.StringVar(self, value='0')
            self.filling_vars.append((mer_var, long_var))
            Tk_.Label(cusp, width=0, background=color, bd=1,).grid(
                row=0, column=0, rowspan=2, sticky=Tk_.NS, padx=4, pady=8)
            ttk.Label(cusp, text='Meridian:').grid(
                row=0, column=1, sticky=Tk_.E)
            ttk.Label(cusp, text='Longitude:').grid(
                row=1, column=1, sticky=Tk_.E)
            meridian = Spinbox(cusp, width=4, textvariable=mer_var,
                from_=-1000, to=1000, increment=1,
                name=':%s:0'%n, validate='focusout',
                validatecommand=(self.register(self.validate_coeff),'%P','%W')
                )
            meridian.bind('<Return>', self.do_filling)
            meridian.grid(row=0, column=2, sticky=Tk_.W, padx=0, pady=3)
            longitude = Spinbox(cusp, width=4, textvariable=long_var,
                from_=-1000, to=1000, increment=1,
                name=':%s:1'%n, validate='focusout',
                validatecommand=(self.register(self.validate_coeff),'%P','%W')
                )
            longitude.bind('<Return>', self.do_filling)
            longitude.grid(row=1, column=2, sticky=Tk_.W, padx=0, pady=3)
            if self.filling_canvas:
                self.filling_canvas.create_window(0, n*cusp_box_height,
                    anchor=Tk_.NW, window=cusp)
            else:
                cusp.grid(row=n, pady=8)
        filling.grid(row=0, column=0, padx=10, pady=10)
        modify = ttk.Labelframe(side_panel, text='Modify', padding=(10, 10))
        self.fill_button = ttk.Button(modify, text='Fill', command=self.do_filling
                ).pack(padx=10, pady=5, expand=True, fill=Tk_.X)
        ttk.Button(modify, text='Retriangulate', command=self.retriangulate
                ).pack(padx=10, pady=5, expand=True, fill=Tk_.X)
        modify.grid(row=1, column=0, padx=10, pady=10, sticky=Tk_.EW)
        create = ttk.Labelframe(side_panel, text='Create', padding=(10, 10))
        ttk.Button(create, text='Drill ...', command=self.drill
                ).pack(padx=10, pady=5, expand=True, fill=Tk_.X)
        ttk.Button(create, text='Cover ...', command=self.cover
                ).pack(padx=10, pady=5, expand=True, fill=Tk_.X)
        create.grid(row=2, column=0, padx=10, pady=10, sticky=Tk_.EW)
        return side_panel

    def build_invariants(self):
        style = self.style
        frame = ttk.Frame(self)
        frame.grid_columnconfigure(1, weight=1)
        self.basic = basic = ttk.LabelFrame(frame, text="Basic Invariants",
                                            padding=(10, 10))
        for n in range(4):
            basic.grid_rowconfigure(n, weight=1)
        self.volume = SelectableText(basic, labeltext='Volume:')
        self.volume.grid(row=0, column=0, padx=5, pady=0, sticky=Tk_.E)
        self.cs = SelectableText(basic, labeltext='Chern-Simons:')
        self.cs.grid(row=1, column=0, padx=5, pady=0, sticky=Tk_.E)
        self.homology = SelectableText(basic, labeltext='H\u2081:')
        self.homology.grid(row=2, column=0, padx=5, pady=0, sticky=Tk_.E)
        self.orientability = SelectableText(basic, labeltext='Orientable:')
        self.orientability.grid(row=3, column=0, padx=5, pady=0, sticky=Tk_.E)
        basic.grid(row=0, column=0, sticky=Tk_.NSEW, padx=10, pady=10)

        self.aka = ttk.LabelFrame(frame, text='Also Known As', padding=(10, 10))
        self.aka_viewer = aka_viewer = ttk.Treeview(
            self.aka,
            selectmode='none',
            height=4,
            columns=['manifold', 'as_link'],
            show='headings')
        aka_viewer.heading('manifold', text='Manifold')
        aka_viewer.column('manifold', stretch=True, width=100)
        aka_viewer.heading('as_link', text='Same Link')
        aka_viewer.column('as_link', stretch=True, width=100)
        aka_viewer.pack(expand=True, fill=Tk_.BOTH)
        self.aka.grid(row=1, column=0, padx=10, pady=10, sticky=Tk_.NSEW)

        self.fundamental = fundamental = ttk.LabelFrame(frame,
            text="Fundamental Group", padding=(10, 10))
        self.pi_one = SelectableMessage(fundamental)
        self.pi_one.grid(row=0, column=0, padx=5, pady=5, sticky=Tk_.W)
        self.pi_one_options = ttk.Frame(fundamental)
        self.simplify_var = Tk_.BooleanVar(fundamental, value=True)
        self.simplify = ttk.Checkbutton(self.pi_one_options,
            variable=self.simplify_var,
            text='simplified presentation',
            command=self.compute_pi_one)
        self.simplify.grid(row=0, column=0, sticky=Tk_.W)
        self.minimize_var = Tk_.BooleanVar(fundamental, value=True)
        self.minimize = ttk.Checkbutton(self.pi_one_options,
            variable=self.minimize_var,
            text='minimal number of generators',
            command=self.compute_pi_one)
        self.minimize.grid(row=1, column=0, sticky=Tk_.W)
        self.gens_change_var = Tk_.BooleanVar(fundamental, value=True)
        self.gens_change = ttk.Checkbutton(self.pi_one_options,
            variable=self.gens_change_var,
            text='fillings may affect generators',
            command=self.compute_pi_one)
        self.gens_change.grid(row=2, column=0, sticky=Tk_.W)
        self.pi_one_options.grid(row=1, column=0,
            padx=10, pady=10, sticky=Tk_.W)
        fundamental.grid(row=0, column=1, rowspan=2, padx=10, pady=10,
                             sticky=Tk_.W)

        self.length_spectrum_frame = ttk.LabelFrame(frame,
            text='Length Spectrum', padding=(10, 10))
        self.length_spectrum_frame.grid_columnconfigure(1, weight=1)
        ttk.Label(self.length_spectrum_frame, text='Cutoff:').grid(
            row=0, column=0, sticky=Tk_.E, pady=5)
        self.length_cutoff = 1.0
        self.cutoff_var=Tk_.StringVar(self, self.length_cutoff)
        self.cutoff_entry = cutoff_entry = ttk.Entry(
            self.length_spectrum_frame,
            takefocus=False,
            width=6,
            textvariable=self.cutoff_var,
            validate='focusout',
            validatecommand=(self.register(self.validate_cutoff),'%P')
            )
        cutoff_entry.bind('<Return>', lambda event : self.focus_set())
        cutoff_entry.grid(row=0, column=1, sticky=Tk_.W, pady=5)
        self.geodesics = geodesics = ttk.Treeview(self.length_spectrum_frame,
            height=6, show='headings', selectmode='none',
            columns=['mult', 'length', 'topology', 'parity'],
            )
        #TODO: compute column widths by measuring text.
        geodesics.heading('mult', text='Mult.')
        geodesics.column('mult', stretch=False, width=60)
        geodesics.heading('length', text='Complex Length')
        geodesics.column('length', stretch=True, minwidth=360)
        geodesics.heading('topology', text='Topology')
        geodesics.column('topology', stretch=True, width=100)
        geodesics.heading('parity', text='Parity')
        geodesics.column('parity', stretch=True, width=80)
        geodesics.grid(row=1, column=0, columnspan=2, sticky=Tk_.EW, padx=10,
                           pady=10)
        self.length_spectrum_frame.grid(row=4, columnspan=2, padx=10, pady=10,
                                            sticky=Tk_.NSEW)
        return frame

    def build_symmetry(self):
        style = self.style
        frame = ttk.Frame(self)
        frame.grid_columnconfigure(0, weight=1)
        self.symmetry = SelectableText(frame, labeltext='Symmetry Group:',
                                       width=30, depth=1)
        self.symmetry.grid(row=0, column=0, pady=20)
        message1 = ttk.Label(frame,
            text='Future releases of SnapPy will show more information here.')
        message2 = ttk.Label(frame,
            text = 'Type SymmetryGroup.<tab> in the command shell to see '
            'what is available.')
        message1.grid(row=1, column=0, pady=(40, 10))
        message2.grid(row=2, column=0)
        return frame

    def build_link(self):
        if self.manifold.DT_code() is not None:
            data = None
            try:
                data = OrthogonalLinkDiagram(self.manifold.link()).plink_data()
            except:
                if self.manifold.LE:
                    data = self.manifold.LE.pickle()
            if data:
                return LinkTab(data, self)

    def build_inside_view(self):
        if not self.manifold.is_orientable():
            text = ("Inside view for non-orientable manifolds such as %s "
                    "is not supported yet.") % self.manifold.name()
            return ttk.Label(self, text = text)

        try:
            # delayed import to avoid cycle
            from .raytracing.inside_viewer import InsideViewer
            self.inside_view = InsideViewer(self, self.manifold,
                fillings_changed_callback=self.update_modeline_and_side_panel)
            self.fillings_changed_callback = self.inside_view.pull_fillings_from_manifold
            return self.inside_view
        except Exception:
            import traceback
            text = ("Could not instantiate inside view. "
                    "Error was:\n\n%s" % traceback.format_exc())
            return ttk.Label(self, text = text)

    def update_menus(self, menubar):
        """Default menus used by the Invariants, Symmetry and Link tabs."""
        menubar.children['help'].activate([])

    def update_modeline(self):
        modeline = self.modeline
        modeline.config(state=Tk_.NORMAL)
        modeline.delete(1.0, Tk_.END)
        if not self.manifold.is_orientable():
            modeline.insert(Tk_.END, 'Non-orientable; ')
        modeline.insert(Tk_.END, '%s tetrahedra; %s'%(
            self.manifold.num_tetrahedra(), self.manifold.solution_type()))
        if len(self.dirichlet) == 0:
             modeline.insert(Tk_.END,
                             '  * failed to compute Dirichlet domain!',
                             'alert')
        modeline.config(state=Tk_.DISABLED)

    def update_modeline_and_side_panel(self):
        self.update_modeline()
        self.update_side_panel()

    def update_current_tab(self, event=None):
        self.update_modeline_and_side_panel()
        tab_name = self.notebook.tab(self.notebook.select(), 'text')
        if tab_name == 'Invariants':
            self.update_menus(self.menubar)
            self.update_invariants()
        elif tab_name == 'Cusp Nbhds':
            self.horoball_viewer.update_menus(self.menubar)
            self.horoball_viewer.view_menu.focus_set()
            if self.horoball_viewer.empty:
                self.update_cusps()
            else:
                self.horoball_viewer.redraw()
        elif tab_name == 'Dirichlet':
            self.dirichlet_viewer.update_menus(self.menubar)
        elif tab_name == 'Link':
            self.update_menus(self.menubar)
            self.link_tab.draw()
        elif tab_name == 'Symmetry':
            self.update_menus(self.menubar)
            self.update_symmetry()
        else:
            self.update_menus(self.menubar)
        self.update_idletasks()

    def update_side_panel(self):
        current_fillings = [c.filling for c in self.manifold.cusp_info()]
        for n, coeffs in enumerate(current_fillings):
            for m in (0,1):
                value = '%g'%coeffs[m]
                value = '0' if value == '-0' else value
                self.filling_vars[n][m].set(value)

    def update_invariants(self):
        manifold = self.manifold
        if not self.recompute_invariants:
            return
        self.orientability.set('Yes' if manifold.is_orientable() else 'No')
        try:
            self.volume.set(repr(manifold.volume()))
        except ValueError:
            self.volume.set('')
        try:
            self.cs.set(repr(manifold.chern_simons()))
        except ValueError:
            self.cs.set('')
        try:
            self.homology.set(repr(manifold.homology()))
        except ValueError:
            self.homology.set('')
        self.compute_pi_one()
        self.update_length_spectrum()
        self.update_dirichlet()
        self.update_aka()
        self.recompute_invariants = False

    def clear_invariants(self):
        self.volume.set('')
        self.cs.set('')
        self.homology.set('')
        self.pi_one.set('')
        self.geodesics.delete(*self.geodesics.get_children())
        self.symmetry.set('')
        self.recompute_invariants = True

    def update_length_spectrum(self):
        OK = True
        try:
            self.length_spectrum = self.manifold.length_spectrum(
                self.length_cutoff)
        except RuntimeError:
            OK = False
            self.length_spectrum = []
        self.geodesics.delete(*self.geodesics.get_children())
        self.cutoff_entry.selection_clear()
        if OK:
            for geodesic in self.length_spectrum:
                parity = '+' if geodesic['parity'].endswith('preserving') else '-'
                self.geodesics.insert('', 'end', values=(
                    geodesic['multiplicity'],
                    Number(geodesic['length'], accuracy=25),
                    geodesic['topology'],
                    parity))
        else:
            self.geodesics.insert('', 'end', values=(
                '?', 'Not computable without a Dirichlet domain', '?', '?'))

    def update_aka(self):
        self._write_aka_info()
        #self.identifier.identify(self.manifold)
        #self.aka_after_id = self.after(100, self._aka_callback)
        M = self.manifold.copy()
        def format_name(N):
            if all(N.cusp_info('is_complete')):
                return N.name()
            return repr(N)
        try:
            mflds = {
                'weak' : [format_name(N) for N in M.identify()],
                'strong' : [format_name(N) for N in M.identify(True)] }
        except ValueError:
            mflds = {
                'weak' : [],
                'strong' : [] }
        self._write_aka_info(mflds)

    def _aka_callback(self):
        # Not used, since the Identifier does not work in a Mac GUI
        self.after_cancel(self.aka_after_id)
        self.aka_after_id = None
        if self.identifier.state == 'finished':
            self._write_aka_info(self.identifier.get())
        else:
            self.aka_after_id = self.after(1000, self._aka_callback)

    def _write_aka_info(self, mflds=None):
        aka_viewer = self.aka_viewer
        try:
            all_items = aka_viewer.get_children()
        except Tk_.TclError: # the widget has been destroyed
            return
        if all_items:
            aka_viewer.delete(*all_items)
        if mflds is None:
            aka_viewer.insert('', 'end', values=('Working ...', ''))
        else:
            strong = set(mflds['strong'])
            weak = set(mflds['weak']) - strong
            for mfld in strong:
                aka_viewer.insert('', 'end', values=(mfld,'Yes'))
            for mfld in weak:
                aka_viewer.insert('', 'end', values=(mfld,'No'))

    def update_dirichlet(self):
        try:
            self.dirichlet = self.manifold.dirichlet_domain().face_list()
        except RuntimeError:
            self.dirichlet = []
        self.update_modeline()
        self.dirichlet_viewer.new_polyhedron(self.dirichlet)

    def update_cusps(self):
        try:
            self.cusp_nbhd = self.manifold.cusp_neighborhood()
        except RuntimeError:
            self.cusp_nbhd = None
        self.horoball_viewer.new_scene(self.cusp_nbhd)
        self.after(100,
                          self.horoball_viewer.cutoff_entry.selection_clear)

    def update_symmetry(self):
        'update_symmetry'
        try:
            self.symmetry_group = self.manifold.symmetry_group()
        except (ValueError, SnapPeaFatalError):
            self.symmetry_group = str('unknown')
        self.symmetry.set(str(self.symmetry_group))

    def validate_coeff(self, P, W):
        tkname, cusp, curve = W.split(':')
        cusp, curve = int(cusp), int(curve)
        try:
            float(P)
        except ValueError:
            var = self.filling_vars[cusp][curve]
            if P == '':
                var.set('0')
            else:
                value = '%g'%self.manifold.cusp_info()[cusp].filling[curve]
                value = '0' if value == '-0' else value
                var.set(value)
            return False
        return True

    def validate_cutoff(self, P):
        try:
            cutoff = float(P)
            if self.length_cutoff != cutoff:
                self.length_cutoff = cutoff
                self.update_length_spectrum()
        except ValueError:
            self.after_idle(
                self.cutoff_var.set, str(self.length_cutoff))
            return False
        return True

    def do_filling(self, event=None):
        filling_spec = [( float(x[0].get() if x[0].get() else 0),
                          float(x[1].get() if x[1].get() else 0) )
                         for x in self.filling_vars]
        self.config(cursor='watch')
        self.clear_invariants()
        self.manifold.dehn_fill(filling_spec)
        current_fillings = [c.filling for c in self.manifold.cusp_info()]
        for n, coeffs in enumerate(current_fillings):
            for m in (0,1):
                value = '%g'%coeffs[m]
                value = '0' if value == '-0' else value
                self.filling_vars[n][m].set(value)
        self.update_cusps()
        self.update_current_tab()
        self.config(cursor='')

        if self.fillings_changed_callback:
            self.fillings_changed_callback()

    def drill(self):
        dialog = Driller(self, self.manifold)
        dialog.go()
        for n in dialog.result:
            self.manifold.drill(n).browse()

    def cover(self):
        dialog = Coverer(self, self.manifold)
        dialog.go()
        for manifold in dialog.result:
            manifold.browse()

    def retriangulate(self):
        self.manifold.randomize()
        self.clear_invariants()
        self.update_current_tab()

    def compute_pi_one(self):
        fun_gp = self.manifold.fundamental_group(
            simplify_presentation=self.simplify_var.get(),
            minimize_number_of_generators=self.minimize_var.get(),
            fillings_may_affect_generators=self.gens_change_var.get())
        self.pi_one.set(repr(fun_gp))

    def save(self, event=None):
        self.manifold.save()

    def close(self, event=None):
        self.destroy()

    def edit_actions(self):
        tab_name = self.notebook.tab(self.notebook.select(), 'text')
        if tab_name in ('Invariants', 'Link', 'Symmetry'):
            try:
                selected =  self.selection_get()
            except:
                selected = False
            if selected:
                return {'Copy' : self.edit_copy}
        return {}

    def edit_copy(self):
        try:
            self.clipboard_clear()
            self.clipboard_append(self.selection_get())
            self.selection_clear()
        except:
            pass

    def test(self):
        self.update_idletasks()
        print('Testing browser')
        self.after(1000, self.notebook.select, self.dirichlet_viewer)
        self.after(2500, self.notebook.select, self.horoball_viewer)
        self.after(4000, self.notebook.select, self.inside_view)
        if self.link_tab:
            self.after(5500, self.notebook.select, self.link_tab.canvas)
            self.after(7000, self.close)
        else:
            self.after(5500, self.close)
        self.wait_window(self)

class Driller(SimpleDialog):
    def __init__(self, master, manifold):
        self.manifold = manifold
        self.num = 0 # make the superclass happy
        self.max_segments = 6
        self.result = []
        style = SnapPyStyle()
        self.root = root = Tk_.Toplevel(master, class_='SnapPy',
                                        bg=style.windowBG)
        title = 'Drill'
        root.title(title)
        root.iconname(title)
        root.bind('<Return>', self.handle_return)
        top_frame = ttk.Frame(self.root)
        top_frame.grid_columnconfigure(0, weight=1)
        top_frame.grid_rowconfigure(2, weight=1)
        msg_font = Font(family=style.font_info['family'],
                        weight='bold',
                        size=int(style.font_info['size']*1.2))
        msg = ttk.Label(top_frame, font=msg_font,
                        text='Choose which curves to drill out:')
        msg.grid(row=0, column=0, pady=10)
        segment_frame = ttk.Frame(top_frame)
        self.segment_var = segment_var = Tk_.StringVar(root)
        segment_var.set(str(self.max_segments))
        ttk.Label(segment_frame, text='Max segments: ').pack(
            side=Tk_.LEFT, padx=4)
        self.segment_entry = segment_entry = ttk.Entry(
            segment_frame,
            takefocus=False,
            width=2,
            textvariable=segment_var,
            validate='focusout',
            validatecommand=(root.register(self.validate_segments),'%P')
            )
        segment_entry.pack(side=Tk_.LEFT)
        segment_frame.grid(row=1, column=0, pady=2)
        self.curves = curves = ttk.Treeview(
            top_frame,
            selectmode='extended',
            columns=['index', 'parity', 'length'],
            show='headings')
        curves.heading('index', text='#')
        curves.column('index', stretch=False, width=20)
        curves.heading('parity', text='Parity')
        curves.column('parity', stretch=False, width=80)
        curves.heading('length', text='Length')
        curves.column('length', stretch=True, width=460)
        curves.bind('<Double-Button-1>', self.drill)
        self.curves.grid(row=2, column=0, padx=6, pady=6, sticky=Tk_.NSEW)
        self.show_curves()
        top_frame.pack(fill=Tk_.BOTH, expand=1)
        button_frame = ttk.Frame(self.root)
        button = ttk.Button(button_frame, text='Drill', command=self.drill,
                            default='active')
        button.pack(side=Tk_.LEFT, padx=6)
        button = ttk.Button(button_frame, text='Cancel', command=self.cancel)
        button.pack(side=Tk_.LEFT, padx=6)
        button_frame.pack(pady=6)
        self.root.protocol('WM_DELETE_WINDOW', self.wm_delete_window)
        self._set_transient(master)

    def show_curves(self):
        self.curves.delete(*self.curves.get_children())
        for curve in self.manifold.dual_curves(max_segments=self.max_segments):
            n = curve['index']
            parity = '+' if curve['parity'] == 1 else '-'
            length = Number(curve['filled_length'], precision=25)
            self.curves.insert( '', 'end', values=(n, parity, length) )

    def handle_return(self, event):
        if event.widget != self.segment_entry:
            self.drill()
        else:
            self.curves.focus_set()

    def drill(self, event=None):
        self.result = [self.curves.index(x) for x in self.curves.selection()]
        self.root.quit()

    def cancel(self):
        self.root.quit()

    def validate_segments(self, P):
        try:
            new_max = int(P)
            if self.max_segments != new_max:
                self.max_segments = new_max
                self.segment_var.set(str(self.max_segments))
                self.show_curves()
        except ValueError:
            self.root.after_idle(
                self.segment_var.set, str(self.max_segments))
            return False
        return True

class Coverer(SimpleDialog):
    def __init__(self, master, manifold):
        self.manifold = manifold.copy()
        self.num = 0 # make the superclass happy
        self.result = []
        style = SnapPyStyle()
        self.root = root = Tk_.Toplevel(master, class_='SnapPy', bg=style.windowBG)
        title = 'Cover'
        root.title(title)
        root.iconname(title)
        root.bind('<Return>', self.handle_return)
        top_frame = ttk.Frame(root)
        top_frame.grid_rowconfigure(2, weight=1)
        top_frame.grid_columnconfigure(0, weight=1)
        top_frame.grid_columnconfigure(1, weight=1)
        msg_font = Font(family=style.font_info['family'],
                        weight='bold',
                        size=int(style.font_info['size']*1.2))
        msg = ttk.Label(top_frame, font=msg_font,
                        text='Choose covering spaces to browse:')
        msg.grid(row=0, column=0, columnspan=3, pady=10)
        degree_frame = ttk.Frame(top_frame)
        degree_frame.grid_columnconfigure(1, weight=1)
        self.degree_var = degree_var = Tk_.StringVar()
        degree_var.trace('w', self.show_covers)
        ttk.Label(degree_frame, text='Degree: ').grid(
            row=0, column=0, sticky=Tk_.E)
        self.degree_option = degree_option = ttk.OptionMenu(
            degree_frame,
            degree_var,
            None,
            *range(2,9)
            )
        degree_option.grid(row=0, column=1)
        self.cyclic_var = cyclic_var = Tk_.BooleanVar()
        cyclic_var.trace('w', self.show_covers)
        cyclic_or_not = ttk.Checkbutton(degree_frame,
                                        variable=cyclic_var,
                                        text='cyclic covers only',
                                       )
        cyclic_or_not.grid(row=0, column=2, padx=6, sticky=Tk_.W)
        self.action = action = ttk.Button(degree_frame, text='Recompute',
                                          command=self.show_covers)
        action.grid(row=0, column=3, padx=8, sticky=Tk_.W)
        degree_frame.grid(row=1, column=0, pady=2, padx=6, sticky=Tk_.EW)
        self.covers = covers = ttk.Treeview(
            top_frame,
            selectmode='extended',
            columns=['index', 'cover_type', 'num_cusps', 'homology'],
            show='headings')
        covers.heading('index', text='')
        covers.column('index', stretch=False, width=40, minwidth=40)
        covers.heading('cover_type', text='Type')
        covers.column('cover_type', stretch=False, width=100)
        covers.heading('num_cusps', text='# Cusps')
        covers.column('num_cusps', stretch=False, width=100, anchor=Tk_.CENTER)
        covers.heading('homology', text='Homology')
        covers.column('homology', stretch=True, width=300)
        covers.bind('<Double-Button-1>', self.choose)
        self.covers.grid(row=2, column=0, columnspan=2, padx=6, pady=6,
                         sticky=Tk_.NSEW)
        top_frame.pack(fill=Tk_.BOTH, expand=1)
        button_frame = ttk.Frame(self.root)
        button_frame.grid_columnconfigure(0, weight=1)
        button_frame.grid_columnconfigure(1, weight=1)
        self.browse = ttk.Button( button_frame, text='Browse', command=self.choose,
            default='active')
        self.browse.grid(row=0, column=0, sticky=Tk_.E, padx=6)
        button = ttk.Button(button_frame, text='Cancel', command=self.cancel)
        button.grid(row=0, column=1, sticky=Tk_.W, padx=6)
        button_frame.pack(pady=6, fill=Tk_.BOTH, expand=1)
        self.root.protocol('WM_DELETE_WINDOW', self.cancel)
        self._set_transient(master)
        degree_var.set('2')
        cyclic_var.set(True)
        self.show_covers()

    def clear_list(self, *args):
        self.covers.delete(*self.covers.get_children())
        self.browse.config(default='normal')
        self.action.config(default='active')
        self.state = 'not ready'

    def show_covers(self, *args):
        self.state = 'ready'
        self.browse.config(default='active')
        self.action.config(default='normal')
        self.covers.delete(*self.covers.get_children())
        degree = int(self.degree_var.get())
        if self.cyclic_var.get():
            self.cover_list = self.manifold.covers(
                degree, cover_type='cyclic')
        else:
            self.cover_list = self.manifold.covers(degree)
        for n, N in enumerate(self.cover_list):
            cusps = repr(N.num_cusps())
            homology = repr(N.homology())
            name = N.name()
            cover_type = N.cover_info()['type']
            self.covers.insert( '', 'end',
                                values=(n, cover_type, cusps, homology))

    def handle_return(self, event):
        if self.state == 'ready':
            self.choose()
        else:
            self.show_covers()

    def choose(self, event=None):
        self.result = [self.cover_list[self.covers.index(x)]
                       for x in self.covers.selection()]
        self.root.destroy()

    def cancel(self):
        self.result = []
        self.root.destroy()

    def go(self):
        self.root.grab_set()
        self.root.wait_window()

if __name__ == '__main__':
    from snappy import Manifold
    root = Tk_.Tk()
    root.withdraw()
    browser1 = Browser(Manifold('m125'), root)
    browser2 = Browser(Manifold('12n345'))
    root.wait_window(browser1.window)
    root.wait_window(browser2.window)
