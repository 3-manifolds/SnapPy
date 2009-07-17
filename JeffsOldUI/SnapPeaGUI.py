#    SnapPeaGUI.py
#
#    A graphical user interface for SnapPea, using the Python wrapper
#    for the SnapPea kernel found in SnapPea.py and SnapPeaC.c.

#    WARNING #1:  This GUI uses Python MegaWidgets.
#    MegaWidgets are great -- but leak memory, in the sense that
#    a Toplevel containing a MegaWidget is never freed.
#    [The guy maintaining PMW plans to fix this.]

#    WARNING #2:  Given that PMW leaks memory already, I've
#    succumbed to the temptation to write memory-leaking code
#    myself.  A SnapPeaGUITriangulation keeps pointers to its
#    constituent SnapPeaGUITriangulationDataPanes, and the
#    SnapPeaGUITriangulationDataPanes keep pointers to the
#    SnapPeaGUITriangulation, so neither will ever be freed.
#    It's embarassing to intentionally write code like this,
#    but for a rapid prototyping tool I guess it's OK.
#    Serious work should take place in the command line interpreter.
#    [Perhaps I could eventually fix this by manually deleting
#    the references.  Ugh.  But until PMW is fixed, there's no point,
#    because I can't verify that everything gets deleted OK.]
#    [See Nathan's comments on this.]


from Tkinter import *
import Pmw

from snappy import *

lighter_color = "blue"

#    A simple dialog to request a file name.
class SaveAsDialog(Toplevel):
    def __init__(self, name = 'untitled'):
        Toplevel.__init__(self)
        self.file_name = ''        # by convention the empty name '' means "don't save"
        self.title('Save As')
        self.entry = Entry(self, relief=SUNKEN)
        self.entry.insert(0, name)
        self.entry.selection_range(0, END)
        self.entry.pack(fill='x', expand=1)
        self.entry.bind('<Key-Return>', self.dismiss)
        self.grab_set()
        self.entry.focus_set()
        self.wait_window()
    def dismiss(self, event):
        self.file_name = self.entry.get()
        self.destroy()


#    The various SnapPeaGUITriangulation<whatever>Panes display
#    information about a SnapPeaGUITriangulation.
#    SnapPeaGUITriangulationDataPane is their common base class.

class SnapPeaGUITriangulationDataPane(Pmw.Group):

    def __init__(self, owner, label):
    
        Pmw.Group.__init__(    self,
                            owner.options_frame.interior(),
                            tag_pyclass=Checkbutton,
                            tag_text=label)
        self.component('tag').select()        
        self.component('tag').configure(command = self.DeleteSelf, selectcolor = '#C00000')
        
        self.owner = owner
    
    def DeleteSelf(self):
        self.owner.DeleteDataPane(self)

class SnapPeaGUITriangulationDehnFillingPane(SnapPeaGUITriangulationDataPane):

    def __init__(self, owner):
    
        SnapPeaGUITriangulationDataPane.__init__(self, owner, 'Dehn filling coefficients')

        theNumCusps = owner.triangulation.num_cusps()
        
        for i in range(theNumCusps):
            theFillButton = Checkbutton(self.interior())
            theFillButton.grid(row=i, column=0)
            theFillButton.select()
            theFillButton['command']        = lambda o=owner, ii=i: o.FillCusp(ii)
            theFillButton['selectcolor']    = '#C00000'

            m = Entry(self.interior(), textvariable=owner.coef[i][0], width=16) 
            m.grid(row=i, column=1)
            m.bind('<Key-Return>', owner.Recompute)

            if owner.triangulation.cusp_info(i)['topology'] == 'torus cusp':
                l = Entry(self.interior(), textvariable=owner.coef[i][1], width=16) 
                l.grid(row=i, column=2)
                l.bind('<Key-Return>', owner.Recompute)    

        Button(    self.interior(),
                text='Recompute',
                font='Helvetica 10',
                command=owner.Recompute).grid(row=theNumCusps, column=1, pady=2)
        Button(    self.interior(),
                text='Reset',
                font='Helvetica 10',
                command=owner.Reset).grid(row=theNumCusps, column=2, pady=2)
        
    def Update(self):
        M = self.owner.triangulation
        for i in range(M.num_cusps()):
            if M.cusp_info(i)['complete?']:
                self.owner.coef[i][0].set('')
                self.owner.coef[i][1].set('')
            else:
                filling = M.cusp_info(i)['filling']
                self.owner.coef[i][0].set(filling[0])
                self.owner.coef[i][1].set(filling[1])

class SnapPeaGUITriangulationVolumePane(SnapPeaGUITriangulationDataPane):

    def __init__(self, owner):
    
        SnapPeaGUITriangulationDataPane.__init__(self, owner, 'volume')
        self.volume = DoubleVar()
        Entry(    self.interior(),
                textvariable=self.volume,
                state=DISABLED,
                justify=CENTER,
                relief=FLAT).pack(pady=2, fill='x', expand=1)
        self.Update()
    
    def Update(self):
        self.volume.set(self.owner.triangulation.volume())

class SnapPeaGUITriangulationHomologyPane(SnapPeaGUITriangulationDataPane):

    def __init__(self, owner):
    
        SnapPeaGUITriangulationDataPane.__init__(self, owner, 'homology')
        self.homology = StringVar()
        Entry(    self.interior(),
                textvariable=self.homology,
                state=DISABLED,
                justify=CENTER,
                relief=FLAT).pack(pady=2, fill='x', expand=1)
        self.Update()
    
    def Update(self):
        self.homology.set(self.owner.triangulation.homology())

class SnapPeaGUITriangulationOrientabilityPane(SnapPeaGUITriangulationDataPane):

    def __init__(self, owner):
    
        SnapPeaGUITriangulationDataPane.__init__(self, owner, 'orientability')
        self.orientability = StringVar()
        Entry(    self.interior(),
                textvariable=self.orientability,
                state=DISABLED,
                justify=CENTER,
                relief=FLAT).pack(pady=2, fill='x', expand=1)
        self.Update()
    
    def Update(self):
        if self.owner.triangulation.is_orientable():
            self.orientability.set('oriented')
        else:
            self.orientability.set('nonorientable')

class SnapPeaGUITriangulationSolutionTypePane(SnapPeaGUITriangulationDataPane):

    def __init__(self, owner):
    
        SnapPeaGUITriangulationDataPane.__init__(self, owner, 'solution type')
        self.solution_type = StringVar()
        Entry(    self.interior(),
                textvariable=self.solution_type,
                state=DISABLED,
                justify=CENTER,
                relief=FLAT).pack(pady=2, fill='x', expand=1)
        self.Update()
    
    def Update(self):
        self.solution_type.set(self.owner.triangulation.solution_type())

class SnapPeaGUITriangulationDrillingPane(SnapPeaGUITriangulationDataPane):

    def __init__(self, owner):
    
        SnapPeaGUITriangulationDataPane.__init__(self, owner, 'drilling')
        
        self.box = Pmw.ScrolledListBox(
                    self.interior(),
                    listbox_height = 6,
                    dblclickcommand = self.DrillCurve)
        self.box.pack(fill='x', expand=1, pady=2)
        self.box.component('listbox')['font'] = 'Courier 12'
    
        self.max_segments = IntVar()
        self.max_segments.set(6)
        self.counter = Pmw.Counter(
                    self.interior(),
                    labelpos = 'w',
                    label_text = 'max. segments',
                    entry_width = 2,
                    entryfield_value = 6,
                    entryfield_validate = {'validator':'integer', 'min':1, 'max':12})
        self.counter.pack()
        self.counter.component('label')['font'] = 'Helvetica 10'
        self.counter.component('entryfield').component('entry')['font'] = 'Helvetica 10'
        self.counter.component('entryfield').component('entry')['textvariable'] = self.max_segments
        self.counter.component('entryfield').configure(modifiedcommand = self.Update)
        
        self.Update()
    
    def Update(self):
        #    The max_segments could temporarily be empty,
        #    e.g. after the user has erased an old value,
        #    but before he/she has entered a new one.
        if self.counter.component('entryfield').valid() != 1:
            self.box.setlist('')
            return
        
        theNumericalData = self.owner.triangulation.dual_curves(self.max_segments.get())
        theTextData = []
        for theCurve in theNumericalData:
            filled = theCurve['filled length']
            complete = theCurve['complete length']
            if not theCurve['complete length']:    # orientation revsersing
                theLine = '%7.5lf           (%7.5lf         )'%(filled.real, complete.real)
                    
                if (filled.imag != 0.0 or complete.imag != 0.0):
                    raise RuntimeError
            else:            # orientation preserving
                theLine = '%7.5lf %8.5lf  (%7.5lf %8.5lf)'%(
                    filled.real, filled.imag, complete.real, complete.imag)
            theTextData.append(theLine)
        self.box.setlist(theTextData)
    
    def DrillCurve(self):

        theSelection = self.box.curselection()

        #    Due to a bug in one version of Python, this list
        #    could contain strings instead of integers.
        #    Here's the workaround from the Tkinter Intro book.
        try:
            theSelection = map(int, theSelection)
        except ValueError:
            pass

        if len(theSelection) > 0:
            self.owner.DrillCurve(theSelection[0], self.max_segments.get())

class SnapPeaGUITriangulationSplittingsPane(SnapPeaGUITriangulationDataPane):

    def __init__(self, owner):
    
        SnapPeaGUITriangulationDataPane.__init__(self, owner, 'splitting surfaces')
        
        self.box = Pmw.ScrolledListBox(
                    self.interior(),
                    listbox_height = 3,
                    dblclickcommand = self.Split)
        self.box.pack(fill='x', expand=1, pady=2)
#        self.box.component('listbox')['font'] = 'Courier 12'
        
        self.Update()
    
    def Update(self):
        self.box.setlist(self.owner.triangulation.normal_surfaces())
    
    def Split(self):

        theSelection = self.box.curselection()

        #    Due to a bug in one version of Python, this list
        #    could contain strings instead of integers.
        #    Here's the workaround from the Tkinter Intro book.
        try:
            theSelection = map(int, theSelection)
        except ValueError:
            pass

        if len(theSelection) > 0:
            thePieces = self.owner.triangulation.split_on_normal_surface(theSelection[0])
            for thePiece in thePieces:
                SnapPeaGUITriangulation(thePiece)

class SnapPeaGUITriangulationFundamentalGroupPane(SnapPeaGUITriangulationDataPane):

    def __init__(self, owner):
    
        SnapPeaGUITriangulationDataPane.__init__(self, owner, 'fundamental group')

        #    Set up the main display for the group presentation.
        
        self.generators = Text(self.interior(),
            state=DISABLED, relief=GROOVE, width=1, height=1)
        self.generators.pack(fill='x', expand=1, padx=2, pady=2)

        self.relations = Pmw.ScrolledText(self.interior(),
            usehullsize=1, hull_width=16, hull_height=80,
            hscrollmode='none', vscrollmode='dynamic')
        self.relations.component('text').configure(state=DISABLED, relief=GROOVE)
        self.relations.pack(fill='x', expand=1, padx=2, pady=2)
        
        self.option_simplify = BooleanVar()
        self.option_simplify.set(1)
        
        self.option_fillings_affect_generators = BooleanVar()
        self.option_fillings_affect_generators.set(1)
        
        self.option_minimize_num_generators = BooleanVar()
        self.option_minimize_num_generators.set(0)
        
        #    Set up the options sub-pane.
        
        self.show_options = BooleanVar()
        self.show_options.set(0)

        self.options = Pmw.Group(    self.interior(),
                                    tag_pyclass=Checkbutton,
                                    tag_text='options')
        self.options.component('tag').configure(
                                    command = self.ToggleOptions,
                                    variable = self.show_options,
                                    selectcolor = '#F08000')
        self.options.pack(fill='x', padx=6, pady=2)

        #    Set up the representation sub-pane.
        
        self.show_representation = BooleanVar()
        self.show_representation.set(0)
        self.word = StringVar()
        self.word.set('a')
        self.matrix_text = None
        self.use_O31 = BooleanVar()
        self.use_O31.set(1)

        self.representation = Pmw.Group(self.interior(),
                                        tag_pyclass=Checkbutton,
                                        tag_text='representation')
        self.representation.component('tag').configure(
                                        command = self.ToggleRepresentation,
                                        variable = self.show_representation,
                                        selectcolor = '#F08000')
        self.representation.pack(fill='x', padx=6, pady=2)
        
        #    Set up the peripheral curves sub-pane.
        
        self.show_peripheral_curves = BooleanVar()
        self.show_peripheral_curves.set(0)

        self.peripheral_curves = Pmw.Group(self.interior(),
                                        tag_pyclass=Checkbutton,
                                        tag_text='peripheral curves')
        self.peripheral_curves.component('tag').configure(
                                        command = self.TogglePeripheralCurves,
                                        variable = self.show_peripheral_curves,
                                        selectcolor = '#F08000')
        self.peripheral_curves.pack(fill='x', padx=6, pady=2)
        
        #    Update the display.
        self.Update()
    
    def destroy(self):

        #    Given the memory management problems currently plagueing
        #    Python Mega Widgets, let's be safe and explicitly cut
        #    loose any structure that keeps a pointer to kernel memory,
        #    just in case this pane isn't properly released.
        self.presentation = None
        
        #    Call the default destroy().
        SnapPeaGUITriangulationDataPane.destroy(self)
    
    def Update(self):
        
        self.UpdatePresentation()

        #    The options display doesn't need explicit updates
        #    because it changes only under direct user control.
        
        if self.show_representation.get() == 1:
            self.UpdateRepresentation()

        if self.show_peripheral_curves.get() == 1:
            self.UpdatePeripheralCurves()
    
    def UpdatePresentation(self):

        self.presentation    = self.owner.triangulation.fundamental_group(
                                self.option_simplify.get(),
                                self.option_fillings_affect_generators.get(),
                                self.option_minimize_num_generators.get())
        
        self.generators.configure(state=NORMAL)
        self.generators.delete(1.0, END)
        self.generators.insert(END, self.presentation.generators())
        self.generators.configure(state=DISABLED)
        
        self.relations.component('text').configure(state=NORMAL)
        self.relations.delete(1.0, END)
        self.relations.insert(END, self.presentation.relators())
        self.relations.component('text').configure(state=DISABLED)

    def ToggleOptions(self):

        if self.show_options.get() == 1:
            self.CreateOptions()
        else:
            self.options_frame.pack_forget()
            self.options_frame.destroy()
            self.options_frame = None
            self.options.interior().configure(width = 1, height = 1)
    
    def CreateOptions(self):
        
        self.options_frame = Frame(self.options.interior())
        self.options_frame.pack()
                        
        thePair1 = Frame(self.options_frame)
        thePair1.pack(fill='x')
        Radiobutton(thePair1,
                    text = 'simplified presentation',
                    variable = self.option_simplify,
                    value = 1,
                    selectcolor = lighter_color,
                    command = self.Update).grid(sticky=W)
        Radiobutton(thePair1,
                    text = 'geometric presentation',
                    variable = self.option_simplify,
                    value = 0,
                    selectcolor = lighter_color,
                    command = self.Update).grid(sticky=W)
        
        thePair2 = Frame(self.options_frame)
        thePair2.pack(fill='x')
        Radiobutton(thePair2,
                    text = 'generators may depend on fillings',
                    variable = self.option_fillings_affect_generators,
                    value = 1,
                    selectcolor = lighter_color,
                    command = self.Update).grid(sticky=W)
        Radiobutton(thePair2,
                    text = 'same generators for all fillings',
                    variable = self.option_fillings_affect_generators,
                    value = 0,
                    selectcolor = lighter_color,
                    command = self.Update).grid(sticky=W)
        
        thePair3 = Frame(self.options_frame)
        thePair3.pack(fill='x')
        Radiobutton(thePair3,
                    text = 'minimize number of generators',
                    variable = self.option_minimize_num_generators,
                    value = 1,
                    selectcolor = lighter_color,
                    command = self.Update).grid(sticky=W)
        Radiobutton(thePair3,
                    text = 'minimize total length of relations',
                    variable = self.option_minimize_num_generators,
                    value = 0,
                    selectcolor = lighter_color,
                    command = self.Update).grid(sticky=W)
    
    def ToggleRepresentation(self):

        if self.show_representation.get() == 1:
            self.CreateRepresentation()
        else:
            self.matrix_text = None
            self.representation_frame.pack_forget()
            self.representation_frame.destroy()
            self.representation_frame = None
            self.representation.interior().configure(width = 1, height = 1)
    
    def CreateRepresentation(self):
        
        self.representation_frame = Frame(self.representation.interior())
        self.representation_frame.pack(fill='x', expand=1)
        self.representation_frame.columnconfigure(0, weight=1)
        self.representation_frame.columnconfigure(1, weight=1)
        
        Pmw.EntryField(    self.representation_frame,
                        entry_textvariable = self.word,
                        validate = {'validator':'alphabetic'},
                        modifiedcommand=self.UpdateRepresentation).grid(
            row=0, column=0, columnspan=2, sticky=W+E, padx=6, pady=6)

        self.matrix_text = Text(self.representation_frame,
                                width=52,
                                height=4,
                                font='Courier 12',
                                state=DISABLED,
                                relief=GROOVE)
        self.matrix_text.grid(row=1, column=0, columnspan=2, sticky=W+E, padx=16, pady=6)

        Radiobutton(self.representation_frame,
                    text = 'O(1,3)',
                    variable = self.use_O31,
                    value = 1,
                    selectcolor = lighter_color,
                    command = self.UpdateRepresentation).grid(row=2, column=0)
        Radiobutton(self.representation_frame,
                    text = 'SL(2,C)',
                    variable = self.use_O31,
                    value = 0,
                    selectcolor = lighter_color,
                    command = self.UpdateRepresentation).grid(row=2, column=1)
        
        self.UpdateRepresentation()
    
    def UpdateRepresentation(self):


        if self.use_O31.get() == 1:
            theText =  repr(self.presentation.O31(self.word.get()))
        else:
            theText =  repr(self.presentation.SL2C(self.word.get()))

  
        self.matrix_text.config(state=NORMAL)
        self.matrix_text.delete(1.0, END)
        self.matrix_text.insert(END, theText)
        self.matrix_text.config(state=DISABLED)
    
    def TogglePeripheralCurves(self):

        if self.show_peripheral_curves.get() == 1:
            self.CreatePeripheralCurves()
        else:
            self.peripheral_curve_entries = None
            self.peripheral_curves_frame.pack_forget()
            self.peripheral_curves_frame.destroy()
            self.peripheral_curves_frame = None
            self.peripheral_curves.interior().configure(width = 1, height = 1)
    
    def CreatePeripheralCurves(self):
        
        self.peripheral_curves_frame = Frame(self.peripheral_curves.interior())
        self.peripheral_curves_frame.pack(fill='x', expand=1)

        self.peripheral_curves_frame.columnconfigure(0, weight=1)
        self.peripheral_curves_frame.columnconfigure(1, weight=1)
        
        self.peripheral_curve_entries = []
        
        for i in range(self.owner.triangulation.num_cusps()):

            theMeridian = Entry(self.peripheral_curves_frame, width=1, relief=GROOVE)
            theMeridian.grid (row = i, column = 0, sticky=W+E)
            
            theLongitude = Entry(self.peripheral_curves_frame, width=1, relief=GROOVE)
            theLongitude.grid(row = i, column = 1, sticky=W+E)
            
            self.peripheral_curve_entries.append([theMeridian, theLongitude])
        
        self.UpdatePeripheralCurves()
    
    def UpdatePeripheralCurves(self):
    
        thePeripheralCurves = self.presentation.peripheral_curves()
        
        for i in range(self.owner.triangulation.num_cusps()):
            
            for j in range(2):    # meridian, longitude
            
                self.peripheral_curve_entries[i][j].config(state=NORMAL)
                self.peripheral_curve_entries[i][j].delete(0, END)
                self.peripheral_curve_entries[i][j].insert(END, thePeripheralCurves[i][j])
                self.peripheral_curve_entries[i][j].config(state=DISABLED)

class SnapPeaGUITriangulationCoreGeodesicsPane(SnapPeaGUITriangulationDataPane):

    def __init__(self, owner):
    
        SnapPeaGUITriangulationDataPane.__init__(self, owner, 'core geodesics')

        self.text = Pmw.ScrolledText(self.interior(),
            usehullsize=1, hull_width=16, hull_height=80,
            hscrollmode='none', vscrollmode='dynamic')
        self.text.component('text').configure(    state=DISABLED,
                                                relief=GROOVE,
                                                font='Courier 12')
        self.text.pack(fill='x', expand=1, padx=2, pady=2)
        
        self.option_menu = Pmw.OptionMenu(    self.interior(),
                                            items=[    'complex length',
                                                    'holonomy',
                                                    'trace in PSL(2,C)',
                                                    'trace squared in PSL(2,C)',
                                                    'eigenvalue in PSL(2,C)'],
                                            initialitem='complex length',
                                            command=self.UpdateFormat)
        self.option_menu.component('menu'      ).configure(font='Helvetica 10')
        self.option_menu.component('menubutton').configure(font='Helvetica 10')
        self.option_menu.pack()

        self.Update()
    
    def Update(self):
        self.UpdateFormat(self.option_menu.getvalue())
    
    def UpdateFormat(self, format):

        self.text.component('text').configure(state=NORMAL)

        self.text.delete(1.0, END)

        for i in range(self.owner.triangulation.num_cusps()):

            self.text.insert(END, '%2d '%i)

            theCore = self.owner.triangulation.core_geodesic(i)

            #    singular index is zero => the Cusp is unfilled
            #    or the Dehn filling coefficients are not integers

            if theCore['singularity index'] > 0:

                if theCore['singularity index'] > 1:
                    self.text.insert(END, '%2d '%theCore['singularity index'])
                else:
                    self.text.insert(END, '   ')

                if format == 'complex length':
                    theValue = theCore['complex length']
                if format == 'holonomy':
                    theValue = theCore['holonomy']
                if format == 'trace in PSL(2,C)':
                    theValue = theCore['trace']
                if format == 'trace squared in PSL(2,C)':
                    theValue = theCore['trace squared']
                if format == 'eigenvalue in PSL(2,C)':
                    theValue = theCore['eigenvalue']
                
                thePrecision = theCore['precision']
                theSpacing = '              '
                if self.owner.triangulation.cusp_is_orientable(i):
                    self.text.insert(END, '%*.*f%s %*.*f%s i'
                        %(    thePrecision + 3, thePrecision, theValue.real, theSpacing[thePrecision:],
                            thePrecision + 3, thePrecision, theValue.imag, theSpacing[thePrecision:]))
                else:
                    self.text.insert(END, '%*.*f%s'
                        %(    thePrecision + 3, thePrecision, theValue.real, theSpacing[thePrecision:]))
                    if theValue.imag != 0.0:
                        raise RuntimeError, 'nonorientable geodesic has nonzero torison'

            self.text.insert(END, '\n')

        self.text.component('text').configure(state=DISABLED)


class SnapPeaGUITriangulationChangePeripheralCurvesPane(SnapPeaGUITriangulationDataPane):

    def __init__(self, owner):
    
        SnapPeaGUITriangulationDataPane.__init__(self, owner, 'change peripheral curves')

        Button(    self.interior(),
                text='shortest curves become meridians',
                command=self.owner.ShortestCurvesBecomeMeridians).pack(padx=2, pady=2)

        Button(    self.interior(),
                text='current Dehn fillings become meridians',
                command=self.owner.CurrentFillingsBecomeMeridians).pack(padx=2, pady=2)

    def Update(self):
        pass

class SnapPeaGUITriangulationRetriangulationPane(SnapPeaGUITriangulationDataPane):

    def __init__(self, owner):
    
        SnapPeaGUITriangulationDataPane.__init__(self, owner, 'retriangulation')

        Label(    self.interior(),
                text='%d tetrahedra'%owner.triangulation.num_tetrahedra()
                ).pack(padx=2, pady=2)

        Button(    self.interior(),
                text='simplify',
                command=self.owner.Simplify).pack(padx=2, pady=2)

        Button(    self.interior(),
                text='randomize',
                command=self.owner.Randomize).pack(padx=2, pady=2)

        Button(    self.interior(),
                text='reverse orientation',
                command=self.owner.Reflect).pack(padx=2, pady=2)

        Button(    self.interior(),
                text='canonize',
                command=self.owner.Canonize).pack(padx=2, pady=2)

    def Update(self):
        pass


class SnapPeaGUITriangulationSymmetryGroupPane(SnapPeaGUITriangulationDataPane):

    def __init__(self, owner):
    
        SnapPeaGUITriangulationDataPane.__init__(self, owner, 'symmetry group')
        
        Label(    self.interior(),
                text='symmetry group of manifold',
                font='Helvetica 12').pack()
        self.manifold_group = StringVar()
        Entry(    self.interior(),
                textvariable=self.manifold_group,
                state=DISABLED,
                justify=CENTER,
                relief=FLAT).pack(pady=2, fill='x', expand=1)
        self.manifold_chirality = StringVar()
        Entry(    self.interior(),
                textvariable=self.manifold_chirality,
                state=DISABLED,
                justify=CENTER,
                relief=FLAT).pack(pady=2, fill='x', expand=1)

        Label(    self.interior(),
                text='symmetry group of associated link',
                font='Helvetica 12').pack()
        self.link_group = StringVar()
        Entry(    self.interior(),
                textvariable=self.link_group,
                state=DISABLED,
                justify=CENTER,
                relief=FLAT).pack(pady=2, fill='x', expand=1)
        self.link_chirality = StringVar()
        Entry(    self.interior(),
                textvariable=self.link_chirality,
                state=DISABLED,
                justify=CENTER,
                relief=FLAT).pack(pady=2, fill='x', expand=1)
        self.knot_invertibility = StringVar()
        Entry(    self.interior(),
                textvariable=self.knot_invertibility,
                state=DISABLED,
                justify=CENTER,
                relief=FLAT).pack(pady=2, fill='x', expand=1)
        
        #    Set up the details sub-pane.
        
        self.show_details = BooleanVar()
        self.show_details.set(0)

##        self.details = Pmw.Group(    self.interior(),
##                                    tag_pyclass=Checkbutton,
##                                    tag_text='group structure')
##        self.details.component('tag').configure(
##                                    command = self.ToggleDetails,
##                                    variable = self.show_details,
##                                    selectcolor = '#F08000')
##        self.details.pack(fill='x', padx=6, pady=2)

        #    Update.
        
        self.Update()
    
    def destroy(self):

        #    Given the memory management problems currently plagueing
        #    Python Mega Widgets, let's be safe and explicitly cut
        #    loose any structure that keeps a pointer to kernel memory,
        #    just in case this pane isn't properly released.
        self.symmetry_group = None
        
        #    Call the default destroy().
        SnapPeaGUITriangulationDataPane.destroy(self)
    
    def Update(self):
    
        self.symmetry_group = {}
        self.symmetry_group['manifold'] = self.owner.triangulation.symmetry_group()
        self.symmetry_group['link'] = self.owner.triangulation.symmetry_group(of_link=True)
        

        self.UpdateMainDisplay()
        
        if self.show_details.get() == 1:
            self.UpdateDetails()    
    
    def UpdateMainDisplay(self):

        if self.symmetry_group['manifold'] != None:
            self.manifold_group.set(repr(self.symmetry_group['manifold']))
            if self.owner.triangulation.is_orientable() == 1:
                if self.symmetry_group['manifold'].is_amphicheiral() == 1:
                    self.manifold_chirality.set('amphicheiral')
                else:
                    self.manifold_chirality.set('chiral')
            else:
                self.manifold_chirality.set('nonorientable')
        else:
            self.manifold_group.set('-')
            self.manifold_chirality.set('')

        if self.symmetry_group['link'] != None:
            self.link_group.set(repr(self.symmetry_group['link']))
            if self.owner.triangulation.is_orientable() == 1:
                if self.symmetry_group['link'].is_amphicheiral() == 1:
                    self.link_chirality.set('amphicheiral')
                else:
                    self.link_chirality.set('chiral')
                if self.owner.triangulation.num_cusps() == 1:
                    if self.symmetry_group['link'].is_invertible_knot() == 1:
                        self.knot_invertibility.set('invertible knot')
                    else:
                        self.knot_invertibility.set('noninvertible knot')
                else:
                    self.knot_invertibility.set('')
            else:
                self.link_chirality.set('nonorientable')
                self.knot_invertibility.set('')
        else:
            self.link_group.set('-')
            self.link_chirality.set('')
            self.knot_invertibility.set('')

    def ToggleDetails(self):

        if self.show_details.get() == 1:
            self.CreateDetails()
        else:
            self.details_frame.pack_forget()
            self.details_frame.destroy()
            self.details_frame = None
            self.details.interior().configure(width = 1, height = 1)
    
    def CreateDetails(self):
        
        self.details_frame = Frame(self.details.interior())
        self.details_frame.pack(fill='x', expand=1)

        Label(    self.details_frame,
                text='commutator',
                font='Helvetica 12').pack()
        self.commutator_text_manifold = StringVar()
        Entry(    self.details_frame,
                textvariable=self.commutator_text_manifold,
                state=DISABLED,
                justify=CENTER,
                relief=FLAT).pack(pady=2, fill='x', expand=1)
        self.commutator_text_link = StringVar()
        Entry(    self.details_frame,
                textvariable=self.commutator_text_link,
                state=DISABLED,
                justify=CENTER,
                relief=FLAT).pack(pady=2, fill='x', expand=1)

        Label(    self.details_frame,
                text='abelianization',
                font='Helvetica 12').pack()
        self.abelianization_text_manifold = StringVar()
        Entry(    self.details_frame,
                textvariable=self.abelianization_text_manifold,
                state=DISABLED,
                justify=CENTER,
                relief=FLAT).pack(pady=2, fill='x', expand=1)
        self.abelianization_text_link = StringVar()
        Entry(    self.details_frame,
                textvariable=self.abelianization_text_link,
                state=DISABLED,
                justify=CENTER,
                relief=FLAT).pack(pady=2, fill='x', expand=1)

        Label(    self.details_frame,
                text='center',
                font='Helvetica 12').pack()
        self.center_text_manifold = StringVar()
        Entry(    self.details_frame,
                textvariable=self.center_text_manifold,
                state=DISABLED,
                justify=CENTER,
                relief=FLAT).pack(pady=2, fill='x', expand=1)
        self.center_text_link = StringVar()
        Entry(    self.details_frame,
                textvariable=self.center_text_link,
                state=DISABLED,
                justify=CENTER,
                relief=FLAT).pack(pady=2, fill='x', expand=1)

        Label(    self.details_frame,
                text='presentation',
                font='Helvetica 12').pack()
        self.presentation_text_manifold = StringVar()
        Entry(    self.details_frame,
                textvariable=self.presentation_text_manifold,
                state=DISABLED,
                justify=CENTER,
                relief=FLAT).pack(pady=2, fill='x', expand=1)
        self.presentation_text_link = StringVar()
        Entry(    self.details_frame,
                textvariable=self.presentation_text_link,
                state=DISABLED,
                justify=CENTER,
                relief=FLAT).pack(pady=2, fill='x', expand=1)
        
        self.UpdateDetails()

    def UpdateDetails(self):
    
        if (self.symmetry_group['manifold'] != None and
            self.symmetry_group['manifold'].is_full_group):
            
            self.commutator_text_manifold.set(
                self.symmetry_group['manifold'].commutator_subgroup().__repr__())
                
            self.abelianization_text_manifold.set(
                self.symmetry_group['manifold'].abelianization().__repr__())
                
            self.center_text_manifold.set(
                self.symmetry_group['manifold'].center().__repr__())
                
            self.presentation_text_manifold.set(
                self.symmetry_group['manifold'].presentation_text())

        else:
            self.commutator_text_manifold.set('-')
            self.abelianization_text_manifold.set('-')
            self.center_text_manifold.set('-')
            self.presentation_text_manifold.set('-')
    
    
        if (self.symmetry_group['link'] != None and
            self.symmetry_group['link'].is_full_group):
            
            self.commutator_text_link.set(
                self.symmetry_group['link'].commutator_subgroup().__repr__())
                
            self.abelianization_text_link.set(
                self.symmetry_group['link'].abelianization().__repr__())
                
            self.center_text_link.set(
                self.symmetry_group['link'].center().__repr__())
                
            self.presentation_text_link.set(
                self.symmetry_group['link'].presentation_text())

        else:
            self.commutator_text_link.set('-')
            self.abelianization_text_link.set('-')
            self.center_text_link.set('-')
            self.presentation_text_link.set('-')

class SnapPeaGUITriangulationTetShapesPane(SnapPeaGUITriangulationDataPane):

    def __init__(self, owner):
    
        SnapPeaGUITriangulationDataPane.__init__(self, owner, 'tetrahedron shapes')

        self.text = Pmw.ScrolledText(self.interior(),
            usehullsize=1, hull_width=360, hull_height=80,
            hscrollmode='none', vscrollmode='dynamic')
        self.text.component('text').configure(    state=DISABLED,
                                                relief=GROOVE,
                                                font='Courier 12')
        self.text.grid(row=0, column=0, columnspan=2, sticky=W+N+E+S, padx=2, pady=2)

        self.format_menu = Pmw.OptionMenu(    self.interior(),
                                            items=[    'rectangular',
                                                    'logarithmic'],
                                            initialitem='rectangular',
                                            command=self.UpdateFromMenu)
        self.format_menu.component('menu'      ).configure(font='Helvetica 10')
        self.format_menu.component('menubutton').configure(font='Helvetica 10')
        self.format_menu.grid(row=1, column=0)

        self.coord_menu = Pmw.OptionMenu(    self.interior(),
                                            items=[    'nicest edge parameters',
                                                    'fixed edge parameters'],
                                            initialitem='nicest edge parameters',
                                            command=self.UpdateFromMenu)
        self.coord_menu.component('menu'      ).configure(font='Helvetica 10')
        self.coord_menu.component('menubutton').configure(font='Helvetica 10')
        self.coord_menu.grid(row=1, column=1)

        self.Update()
    
    def UpdateFromMenu(self, menu_selection):
        self.Update()
    
    def Update(self):

        self.text.component('text').configure(state=NORMAL)

        self.text.delete(1.0, END)

        alignment =  self.coord_menu.getvalue() == 'fixed edge parameters'
        theTetShapes = self.owner.triangulation.tetrahedra_shapes(fixed_alignment=alignment)
                
        theFormat = self.format_menu.getvalue()
        
        for i in range(self.owner.triangulation.num_tetrahedra()):
            
                
            self.text.insert(END, '%2d '%i)
            
            if theFormat == 'rectangular':
            
                self.text.insert(END, '%*.*f%s %*.*f'%(
                                    
                                    6 + theTetShapes[i]['precisions'][0],
                                    theTetShapes[i]['precisions'][0],
                                    theTetShapes[i]['rect'].real,
                                    '                '[:(16 - theTetShapes[i]['precisions'][0])],
                                    
                                    6 + theTetShapes[i]['precisions'][1],
                                    theTetShapes[i]['precisions'][1],
                                    theTetShapes[i]['rect'].imag))


            elif theFormat == 'logarithmic':
            
                self.text.insert(END, '%*.*f%s %*.*f'%(

                                      6 + theTetShapes[i]['precisions'][2],
                                    theTetShapes[i]['precisions'][2],
                                    theTetShapes[i]['log'].real,
                                    '                '[:(16 - theTetShapes[i]['precisions'][2])],
                                    
                                    6 + theTetShapes[i]['precisions'][3],
                                    theTetShapes[i]['precisions'][3],
                                    theTetShapes[i]['log'].imag))

            else:
                raise RuntimeError, 'bad tet shape format'

            self.text.insert(END, '\n')

        self.text.component('text').configure(state=DISABLED)            

class SnapPeaGUITriangulationDirichletPane(SnapPeaGUITriangulationDataPane):

    def __init__(self, owner):
    
        SnapPeaGUITriangulationDataPane.__init__(self, owner, 'Dirichlet domain (.off)')

        self.text = Pmw.ScrolledText(self.interior(),
            usehullsize=1, hull_width=16, hull_height=80,
            hscrollmode='none', vscrollmode='dynamic')
        self.text.component('text').configure(state=DISABLED, relief=GROOVE, font='Courier 12')
        self.text.pack(fill='x', expand=1, padx=2, pady=2)

        self.basepoint_menu = Pmw.OptionMenu(    self.interior(),
                                            items=[    'basepoint on edge',
                                                    'basepoint at centroid'],
                                            initialitem='basepoint at centroid',
                                            command=self.Update)
        self.basepoint_menu.component('menu'      ).configure(font='Helvetica 10')
        self.basepoint_menu.component('menubutton').configure(font='Helvetica 10')
        self.basepoint_menu.pack()
        
        self.maximize_inj_rad = IntVar()
        self.maximize_inj_rad.set(0)
        Checkbutton(    self.interior(),
                        text='maximize injectivity radius',
                        var=self.maximize_inj_rad,
                        command=self.Update,
                        selectcolor=lighter_color).pack()

        #    Set up the face_pairings sub-pane.
        
        self.show_face_pairings = BooleanVar()
        self.show_face_pairings.set(0)

        self.face_pairings = Pmw.Group(    self.interior(),
                                    tag_pyclass=Checkbutton,
                                    tag_text='face pairings')
        self.face_pairings.component('tag').configure(
                                    command = self.ToggleFacePairings,
                                    variable = self.show_face_pairings,
                                    selectcolor = '#F08000')
        self.face_pairings.pack(fill='x', padx=6, pady=2)
        self.face_pairings_text = None
        
        #    Set up the basepoint displacement sub-pane.
        
        self.dx = DoubleVar()
        self.dx.set(0.0)
        self.dy = DoubleVar()
        self.dy.set(0.0)
        self.dz = DoubleVar()
        self.dz.set(0.0)
        
        self.show_displacement = BooleanVar()
        self.show_displacement.set(0)

        self.displacement = Pmw.Group(    self.interior(),
                                    tag_pyclass=Checkbutton,
                                    tag_text='basepoint displacement')
        self.displacement.component('tag').configure(
                                    command = self.ToggleDisplacement,
                                    variable = self.show_displacement,
                                    selectcolor = '#F08000')
        self.displacement.pack(fill='x', padx=6, pady=2)

        self.Update()
    
    def Update(self, menu_choice=''):
    
        theDirichletDomain = self.owner.triangulation.Dirichlet(
            self.basepoint_menu.get() == 'basepoint at centroid',
            self.maximize_inj_rad.get(),
            (self.dx.get(), self.dy.get(), self.dz.get()))
        
        self.text.component('text').configure(state=NORMAL)
        self.text.delete(1.0, END)
        self.text.insert(END, theDirichletDomain.off())
        self.text.component('text').configure(state=DISABLED)
        
        if self.face_pairings_text != None:
            self.face_pairings_text.component('text').configure(state=NORMAL)
            self.face_pairings_text.delete(1.0, END)
            theFacePairings = theDirichletDomain.face_pairings()
            for i in range(len(theFacePairings)):
                for j in range(4):
                    for k in range(4):
                        self.face_pairings_text.insert(END, ' %11.6f'%theFacePairings[i][j][k])
                    self.face_pairings_text.insert(END, '\n')
                self.face_pairings_text.insert(END, '\n')
            self.face_pairings_text.component('text').configure(state=DISABLED)
        
        #    self.displacement requires no update.

    def ToggleFacePairings(self):
        if self.show_face_pairings.get() == 1:
            self.face_pairings_text = Pmw.ScrolledText(
                self.face_pairings.interior(),
                usehullsize=1, hull_width=16, hull_height=120,
                hscrollmode='none', vscrollmode='dynamic')
            self.face_pairings_text.component('text').configure(state=DISABLED, relief=GROOVE, font='Courier 12')
            self.face_pairings_text.pack(fill='x', expand=1, padx=2, pady=2)
            self.face_pairings_text.pack()
            self.Update()
        else:
            self.face_pairings_text.pack_forget()
            self.face_pairings_text.destroy()
            self.face_pairings_text = None
            self.face_pairings.interior().configure(width = 1, height = 1)
    
    def ToggleDisplacement(self):
        if self.show_displacement.get() == 1:
            self.displacement_frame = Frame(self.displacement.interior())
            self.displacement_frame.pack(fill='x', expand=1)

            theEntryDx = Entry(self.displacement_frame, textvariable=self.dx, width=8)
            theEntryDx.grid(row=0, column=0, sticky=W+E, padx=6)
            theEntryDx.bind('<Key-Return>', self.Update)

            theEntryDy = Entry(self.displacement_frame, textvariable=self.dy, width=8)
            theEntryDy.grid(row=0, column=1, sticky=W+E, padx=6)
            theEntryDy.bind('<Key-Return>', self.Update)

            theEntryDz = Entry(self.displacement_frame, textvariable=self.dz, width=8)
            theEntryDz.grid(row=0, column=2, sticky=W+E, padx=6)
            theEntryDz.bind('<Key-Return>', self.Update)
        else:
            self.displacement_frame.pack_forget()
            self.displacement_frame.destroy()
            self.displacement_frame = None
            self.displacement.interior().configure(width = 1, height = 1)
        
        

#    SnapPeaGUITriangulation represents a triangulation.

class SnapPeaGUITriangulation(Pmw.MegaToplevel):
    
    def __init__(self, triangulation):

        Pmw.MegaToplevel.__init__(self)
        self.title(triangulation.name())
        self.triangulation = triangulation
        
        self.InitDehnCoefficients()

        self.menubar = Frame(self.component('hull'), relief=RAISED, bd=2)
        self.menubar.pack(side=TOP, fill='x')

        theFileButton = Menubutton(self.menubar, text='File', underline=0)
        theFileButton.pack(side=LEFT)
        theFileMenu = Menu(theFileButton)
        theFileMenu.add_command(    label='Save As...',    command=self.SaveAs)
        theFileMenu.add_command(    label='Clone',        command=self.Clone)
        theFileMenu.add_command(    label='Close',        command=self.destroy)
        theFileButton['menu'] = theFileMenu

        theViewButton = Menubutton(self.menubar, text='View', underline=0)
        theViewButton.pack(side=LEFT)
        theViewMenu = Menu(theViewButton)
        theViewMenu.add_command(
            label='Dehn filling coefficients',
            command=lambda s=self: s.MakeDataPane(SnapPeaGUITriangulationDehnFillingPane))
        theViewMenu.add_command(
            label='volume',
            command=lambda s=self: s.MakeDataPane(SnapPeaGUITriangulationVolumePane))
        theViewMenu.add_command(
            label='homology',
            command=lambda s=self: s.MakeDataPane(SnapPeaGUITriangulationHomologyPane))
        theViewMenu.add_command(
            label='orientability',
            command=lambda s=self: s.MakeDataPane(SnapPeaGUITriangulationOrientabilityPane))
        theViewMenu.add_command(
            label='solution type',
            command=lambda s=self: s.MakeDataPane(SnapPeaGUITriangulationSolutionTypePane))
        theViewMenu.add_command(
            label='drilling',
            command=lambda s=self: s.MakeDataPane(SnapPeaGUITriangulationDrillingPane))
#        theViewMenu.add_command(
#            label='splittings',
#            command=lambda s=self: s.MakeDataPane(SnapPeaGUITriangulationSplittingsPane))
        theViewMenu.add_command(
            label='fundamental group',
            command=lambda s=self: s.MakeDataPane(SnapPeaGUITriangulationFundamentalGroupPane))
#        theViewMenu.add_command(
#            label='core geodesics',
#            command=lambda s=self: s.MakeDataPane(SnapPeaGUITriangulationCoreGeodesicsPane))
        theViewMenu.add_command(
            label='change peripheral curves',
            command=lambda s=self: s.MakeDataPane(SnapPeaGUITriangulationChangePeripheralCurvesPane))
        theViewMenu.add_command(
            label='retriangulation',
            command=lambda s=self: s.MakeDataPane(SnapPeaGUITriangulationRetriangulationPane))
        theViewMenu.add_command(
            label='symmetry group',
            command=lambda s=self: s.MakeDataPane(SnapPeaGUITriangulationSymmetryGroupPane))
        theViewMenu.add_command(
            label='tetrahedron shapes',
            command=lambda s=self: s.MakeDataPane(SnapPeaGUITriangulationTetShapesPane))
#        theViewMenu.add_command(
#            label='Dirichlet domain (.off format)',
#            command=lambda s=self: s.MakeDataPane(SnapPeaGUITriangulationDirichletPane))
        theViewButton['menu'] = theViewMenu
        
        self.options_frame = Pmw.ScrolledFrame(self.component('hull'),
            usehullsize=1, hull_width=400, hull_height=450,
            hscrollmode='none', vscrollmode='dynamic',
            horizflex='expand', vertflex='fixed')
        self.options_frame.pack(fill='both', expand=1, padx=6, pady=6)

        self.displays = []
        self.MakeDataPane(SnapPeaGUITriangulationDehnFillingPane)
        self.MakeDataPane(SnapPeaGUITriangulationSolutionTypePane)
        self.MakeDataPane(SnapPeaGUITriangulationVolumePane)
        self.MakeDataPane(SnapPeaGUITriangulationHomologyPane)
        
    #    temporary __del__ method to investigate weird behavior
    #    of embedded Python MegaWidgets.
#    def __del__(self):
#        print 'deleting manifold window'

    def destroy(self):
        #    Due to a bug in Python Mega Widgets, a Pmw.MegaToplevel
        #    is never freed.  (Most likely the problem is circular
        #    reference chains.)  So if we are about to be destroyed,
        #    we should at least free our Triangulation.
        #    (After PMW fixes their bug, I'll probably have to free
        #    my own panes manually to get this SnapPeaGUITriangulation
        #    to be freed.  At that point it'll no longer be necessary
        #    to manually free the Triangulation.)
        self.triangulation = None
        Pmw.MegaToplevel.destroy(self)
        
    def InitDehnCoefficients(self):
        #    Set up shadow variables for the Dehn filling coefficients.
        self.coef = []
        for i in range(self.triangulation.num_cusps()):
            self.coef.append([DoubleVar(), DoubleVar()])
            self.coef[i][0].set(self.triangulation.cusp_info(i)['filling'][0])
            self.coef[i][0].set(self.triangulation.cusp_info(i)['filling'][1])
            if self.coef[i][0].get() == 0.0 and self.coef[i][1].get() == 0.0:
                self.coef[i][0].set('')
                self.coef[i][1].set('')
    
    def MakeDataPane(self, DisplayType):
        theDisplay = DisplayType(self)
        theDisplay.pack(fill='x', padx=6, pady=6)
        self.displays.append(theDisplay)

    def DeleteDataPane(self, pane):
        self.displays.remove(pane)
        pane.destroy()
        
    def Recompute(self, event=None):
    
        for i in range(self.triangulation.num_cusps()):
        
            #    Convert nonnumeric entries (including '') to (0,0).
            #    (I wish I knew how to do this more elegantly.)
            for j in [0,1]:
                try:
                    self.coef[i][j].get()
                except:
                    self.coef[i][j].set(0.0)
            
            if (
#                    self.coef[i][0].get() == '' or self.coef[i][1].get() == '' or
                    (self.coef[i][0].get() == 0.0 and self.coef[i][1].get() == 0.0)):
                self.triangulation.dehn_fill( (0,0), i)
                self.coef[i][0].set('')
                self.coef[i][1].set('')
            else:
                self.triangulation.dehn_fill( (self.coef[i][0].get(), self.coef[i][1].get()), i)
                self.coef[i][0].set(self.coef[i][0].get())
                self.coef[i][1].set(self.coef[i][1].get())

        self.Update()
    
    def Reset(self, event=None):
        M = self.triangulation
        n = M.num_cusps()
        for i in range(n):
            for j in [0,1]:
                self.coef[i][j].set('')
        
        M.dehn_fill(n*[ (0,0) ])

        self.Update()
    
    #    Update() handles mild changes to the manifold,
    #    such as new Dehn filling coefficients and a new solution.
    def Update(self):
        for theDisplay in self.displays:
            theDisplay.Update()
    
    #    Overhaul() handles severe changes to the manifold,
    #    such as changing the combinatorics of the triangulation.
    def Overhaul(self):

        #    Reinitialize the shadow variables
        #    for the Dehn filling coefficients.
        self.InitDehnCoefficients()
        
        #    Trash the old display panes, and
        #    create new ones of the same types.
        theDisplayClasses = []
        for theDisplay in self.displays:
            theDisplayClasses.append(theDisplay.__class__)
            theDisplay.destroy()
        self.displays = []
        for theClass in theDisplayClasses:
            self.MakeDataPane(theClass)
    
    def SaveAs(self):
        theFileName = SaveAsDialog(self.triangulation.name()).file_name
        if (theFileName != ''):
            self.triangulation.save(theFileName)
            self.title(theFileName)
    
    def Clone(self):

        theCopy = SnapPeaGUITriangulation(self.triangulation.clone())

        theDisplayClasses = []
        for theDisplay in self.displays:
            theDisplayClasses.append(theDisplay.__class__)
        theCopy.SetDisplayClasses(theDisplayClasses)

    def SetDisplayClasses(self, new_display_classes):
        for theDisplay in self.displays:
            theDisplay.destroy()
        self.displays = []
        for theClass in new_display_classes:
            self.MakeDataPane(theClass)        
    
    def FillCusp(self, i):
        self.Recompute()
        self.triangulation.fill_cusp(i)
        self.Overhaul()

    def DrillCurve(self, i, max_segments):
        self.triangulation = self.triangulation.drill(i, max_segments)
        self.Overhaul()
    
    def ShortestCurvesBecomeMeridians(self):
        self.triangulation.set_peripheral_curves('shortest')
        self.Update()
    
    def CurrentFillingsBecomeMeridians(self):
        self.triangulation.set_peripheral_curves('fillings')
        self.Update()
    
    def Simplify(self):
        self.triangulation.simplify()
        self.Overhaul()
    
    def Randomize(self):
        self.triangulation.randomize()
        self.Overhaul()
    
    def Reflect(self):
        self.triangulation.reverse_orientation()
        self.Overhaul()
    
    def Canonize(self):
        self.triangulation.canonize()
        self.Overhaul()


def ManifoldFromFile(event=None):
    if len(theFileName.get()) > 0:
        SnapPeaGUITriangulation(Manifold(theFileName.get()))
    else:
        raise ValueError, 'The file name is empty.'

def ManifoldFromCensus(event=None):
    theCensusNumber            = [5, 6, 6, 7, 7]
    theCensusOrientability    = [1, 1, 0, 1, 0]
    SnapPeaGUITriangulation(Triangulation(    theCensusNumber[theCuspedCensus.get()],
                                            theCensusOrientability[theCuspedCensus.get()],
                                            theCensusIndex.get()))

if __name__ == '__main__':
    root = Tk()
    Pmw.initialise(root, fontScheme='pmw1')
    root.title('SnapPea (experimental version)')

    theFileGroup = Pmw.Group(root, tag_text='manifold from descriptor')
    theFileGroup.pack(fill='x', expand=0, padx=6, pady=6)
    theFileName = StringVar()
    theFileName.set('m125')
    theFileNameEntry = Entry(    theFileGroup.interior(),
                                textvariable=theFileName,
                                width=24)
    theFileNameEntry.pack(fill='x', expand=0, padx=6, pady=6)
    theFileNameEntry.bind('<Key-Return>', ManifoldFromFile)
    Button(    theFileGroup.interior(),
            text='Create',
            command=ManifoldFromFile,
            font='Helvetica 10').pack()

    Button(root, text='Quit SnapPea', command=root.quit, font='Helvetica 12').pack()

    root.mainloop()
