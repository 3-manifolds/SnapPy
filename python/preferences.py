import os, sys, re, time
try:
    import plistlib
except ImportError:
    from . import plistlib
from .gui import *

class Preferences:
    def __init__(self, text_widget):
        self.text_widget = text_widget
        self.prefs_dict = {
            'autocall' : False,
            'automagic' : False,
            'font' : self.current_font_tuple(),
            'cusp_horoballs' : True,
            'cusp_triangulation' : True,
            'cusp_ford_domain' : True,
            'cusp_parallelogram' : True,
            'cusp_labels' : True,
            'cusp_cutoff' : '0.1000',
            'keyboard' : 'QWERTY'}
        self.cache = {}
        self.cache_prefs()
        self.find_prefs()
        self.read_prefs()
        
    def __getitem__(self, x):
        return self.prefs_dict[x]
    
    def __setitem__(self, x, y):
        self.prefs_dict[x] = y

    def __repr__(self):
        return str(self.prefs_dict)

    def items(self):
        return self.prefs_dict.items()

    def get(self, key, default):
        return self.prefs_dict.get(key, default)

    def current_font_dict(self):
        font_string = self.text_widget.cget('font')
        return Font(font=font_string).actual()

    def current_font_tuple(self):
        font = self.current_font_dict()
        style = '%s %s'%(font['weight'], font['slant'])
        return (font['family'], font['size'], style) 

    def find_prefs(self):
        if sys.platform == 'darwin':
            home = os.environ['HOME']
            self.prefs_file = os.path.join(home, 'Library',
                                           'Preferences',
                                           'edu.t3m.SnapPy.plist')
        elif sys.platform == 'linux2' or sys.platform == 'linux':
            home = os.environ['HOME']
            self.prefs_file = os.path.join(home, '.SnapPy.plist')
        elif sys.platform == 'win32':
            home = os.environ['USERPROFILE']
            self.prefs_file = os.path.join(home, '.SnapPy.plist')
        else:
            self.prefs_file = None

    def read_prefs(self):
        if self.prefs_file:
            try:
                if hasattr(plistlib, 'load'):
                    with open(self.prefs_file, 'rb') as file:
                        self.prefs_dict.update(plistlib.load(file))
                else:
                    self.prefs_dict.update(plistlib.readPlist(self.prefs_file))
                # plistlib screws up tuples
                self.prefs_dict['font'] = tuple(self.prefs_dict['font'])
            except IOError:
                pass

    def write_prefs(self):
        if self.prefs_file:
            plistlib.writePlist(self.prefs_dict, self.prefs_file)

    def cache_prefs(self):
        self.cache.update(self.prefs_dict)

    def revert_prefs(self):
        self.prefs_dict.update(self.cache)
        self.apply_prefs()

    def changed(self):
        return [key for key in self.cache.keys() if 
                self.cache[key] != self.prefs_dict[key]]

    # Override this in a subclass.
    def apply_prefs(self):
        self.text_widget.config(font=self.prefs_dict['font'])
        print(self.prefs_dict)

class PreferenceDialog(Dialog):
    def __init__(self, parent, prefs, title='SnapPy Preferences'):
        self.parent = parent
        self.style = SnapPyStyle()
        self.prefs = prefs
        self.prefs.cache_prefs()
        self.okay = False
        Tk_.Toplevel.__init__(self, master=parent, class_='snappy')
        frame = ttk.Frame(self, padding=(0, 10, 0, 0))
        self.title(title)
        self.notebook = notebook = ttk.Notebook(frame)
        self.build_font_pane(notebook)
        self.build_shell_pane(notebook)
        self.build_cusp_pane(notebook)
        self.build_inside_pane(notebook)
        notebook.add(self.font_frame, text='Font')
        notebook.add(self.shell_frame, text='Shell')
        notebook.add(self.cusp_frame, text='Cusps')
        notebook.add(self.inside_frame, text='Inside')
        notebook.grid(row=0, column=0)
        self.buttonbox()
        notebook.pack(padx=0, pady=0)
        frame.pack()
        self.button_frame.pack(fill=Tk_.X)
        self.protocol('WM_DELETE_WINDOW', self.cancel)
        self.wait_window(self)

    def validate(self):
        cutoff = self.cutoff.get()
        try:
            float(cutoff)
            self.prefs['cusp_cutoff'] = cutoff 
        except ValueError:
            showerror('Invalid input',
                      'Please enter a number for the cutoff.')
            return False
        self.prefs.apply_prefs()
        self.okay = True
        return True

    def revert(self):
        self.prefs.revert_prefs()
        self.cancel()

    def buttonbox(self):
        self.button_frame = box = ttk.Frame(self, padding=(0, 0, 0, 20))
        box.grid_columnconfigure(0, weight=1)
        box.grid_columnconfigure(2, weight=1)
        OK = ttk.Button(box, text="OK", width=10, command=self.ok,
                        default=Tk_.ACTIVE)
        OK.grid(row=0, column=0, sticky=Tk_.NE, padx=5)
        Apply = ttk.Button(box, text="Apply", width=10,
                           command=self.prefs.apply_prefs)
        Apply.grid(row=0, column=1, sticky=Tk_.N, padx=5)
        Cancel = ttk.Button(box, text="Cancel", width=10, command=self.revert)
        Cancel.grid(row=0, column=2, sticky=Tk_.NW, padx=5)
        self.bind("<Return>", lambda event : OK.focus_set())
        self.bind("<Escape>", self.cancel)

    def show_body(self):
        self.body_frame.grid(row=1, column=0, padx=5, pady=5,
                             sticky=Tk_.N + Tk_.S + Tk_.W + Tk_.E)
        self.body_frame.focus_set()

    def build_font_pane(self, master):
        current_font = self.prefs.current_font_dict()
        groupBG = self.style.groupBG
        self.font_frame = font_frame = ttk.Frame(master)
        font_frame.columnconfigure(2, weight=1)
        font_frame.columnconfigure(0, weight=1)
        self.list_frame = list_frame = ttk.Frame(font_frame)
        self.font_list = font_list = Tk_.Listbox(list_frame,
                                                 selectmode=Tk_.SINGLE)
        self.families = families = [f for f in font_families()
                                    if f[0] != '@'] # omit vertical fonts
        families.sort()
        for family in families:
            font_list.insert(Tk_.END, family)
        self.current_family = families.index(current_font['family'])
        font_list.selection_set(self.current_family)
        font_list.see(self.current_family)
        font_list.grid(row=0, column=0, pady=(20,30))
        font_scroller = ttk.Scrollbar(list_frame,
                                      command=font_list.yview)
        font_list.config(yscrollcommand=font_scroller.set)
        font_scroller.grid(row=0, column=1, sticky=Tk_.N + Tk_.S, pady=(20,30))
        list_frame.grid(rowspan=6, column=0)
        
        label = ttk.Label(self.font_frame, text='Size: ')
        label.grid(row=0, column=1, sticky=Tk_.E, pady=(20,0))
        self.font_sizer = sizer = Spinbox(font_frame, from_=10, to_=36, width=4,
                                             command=self.set_font_sample,
                                             )
        sizer.delete(0,2)
        sizer.insert(0, current_font['size'] )
        sizer.bind('<Return>', self.set_font_sample)
        sizer.bind('<Tab>', self.set_font_sample)
        self.current_size =  current_font['size']
        sizer.grid(row=0, column=2, sticky=Tk_.W, pady=(20,0))
        label = ttk.Label(font_frame, text='Weight: ')
        label.grid(row=1, column=1, sticky=Tk_.E)
        self.font_weight = weight = Tk_.StringVar(value=current_font['weight'])
        radio = ttk.Radiobutton(font_frame, text='normal', variable=weight,
                                    value='normal', command=self.set_font_sample)
        radio.grid(row=1, column=2, sticky=Tk_.W)
        radio = ttk.Radiobutton(font_frame, text='bold', variable=weight,
                                    value='bold', command=self.set_font_sample)
        radio.grid(row=2, column=2, sticky=Tk_.W)
        label = ttk.Label(font_frame, text='Slant: ')
        label.grid(row=3, column=1, sticky=Tk_.E)
        self.font_slant = slant = Tk_.StringVar(value=current_font['slant'])
        radio = ttk.Radiobutton(font_frame, text='roman', variable=slant,
                                 value='roman', command=self.set_font_sample)
        radio.grid(row=3, column=2, sticky=Tk_.W)
        radio = ttk.Radiobutton(font_frame, text='italic', variable=slant,
                                    value='italic', command=self.set_font_sample)
        radio.grid(row=4, column=2, sticky=Tk_.W)
        self.sample = sample = Tk_.Text(self.font_frame,
                                        width=40, height=6,
                                        highlightthickness=0,
                                        relief=Tk_.FLAT,
                                        font=self.prefs['font'])
        self.sample.bind('<Button-1>', lambda event: 'break')
        self.sample.insert(Tk_.INSERT, '\nABCDEFGHIJKLMNOPQRSTUVWXYZ\n'\
                                           'abcdefghijklmnopqrstuvwxyz')
        sample.tag_add('all', '1.0', Tk_.END)
        self.set_font_sample()
        font_list.bind('<ButtonRelease-1>', self.set_font_sample)
        self.sample.grid(row=6, columnspan=4, padx=10, pady=10, sticky=Tk_.E+Tk_.W)

    def get_font(self):
        selection = self.font_list.curselection()
        if selection:
            index = int(selection[0])
            self.current_family = index
        else:
            index = self.current_family
        family = self.families[index]
        try:
            size = int(self.font_sizer.get())
        except ValueError:
            size = self.current_size
            self.font_sizer.delete(0, Tk_.END)
            self.font_sizer.insert(0, str(size))
        self.font_list.selection_set(self.current_family)
        style = '%s %s'%(self.font_weight.get(), self.font_slant.get())
        return (family, size, style.strip()) 

    def set_font_sample(self, event=None):
        new_font = self.get_font()
        self.prefs['font'] = new_font
        self.sample.tag_config('all', justify=Tk_.CENTER,
                               font=new_font) 

    def build_shell_pane(self, master):
        groupBG = self.style.groupBG
        self.autocall = Tk_.BooleanVar(value=self.prefs['autocall'])
        self.automagic = Tk_.BooleanVar(value=self.prefs['automagic'])
        self.update_idletasks()
        self.shell_frame = shell_frame = ttk.Frame(master)
        shell_frame.rowconfigure(3, weight=1)
        shell_frame.columnconfigure(0, weight=1)
        shell_frame.columnconfigure(3, weight=1)
        # Keep the height the same as the height of the font pane.
        strut = ttk.Frame(shell_frame, width=1)
        strut.grid(rowspan=5, column=0)
        next_label = ttk.Label(shell_frame, anchor=Tk_.W,
            text='Which IPython features would you like to enable?')
        next_label.grid(row=0, column=1, columnspan=2, sticky=Tk_.W, pady=(20,0))
        next_check = ttk.Checkbutton(shell_frame, variable = self.autocall,
                                    text='IPython autocall',
                                    command=self.set_autocall)
        next_check.grid(row=1, column=1, sticky=Tk_.W, pady=(10,0))
        next_check = ttk.Checkbutton(shell_frame, variable = self.automagic,
                                    text='IPython automagic',
                                    command=self.set_automagic)
        next_check.grid(row=2, column=1, sticky=Tk_.W, pady=(5,0))

    def set_autocall(self):
        self.prefs['autocall'] = self.autocall.get()

    def set_automagic(self):
        self.prefs['automagic'] = self.automagic.get()

    def build_cusp_pane(self, master):
        groupBG = self.style.groupBG
        self.cusp_frame = cusp_frame = ttk.Frame(master)
        self.horoballs = Tk_.BooleanVar(value=self.prefs['cusp_horoballs'])
        self.triangulation = Tk_.BooleanVar(value=self.prefs['cusp_triangulation'])
        self.ford = Tk_.BooleanVar(value=self.prefs['cusp_ford_domain'])
        self.labels = Tk_.BooleanVar(value=self.prefs['cusp_labels'])
        self.parallelogram = Tk_.BooleanVar(value=self.prefs['cusp_parallelogram'])
        self.cutoff = Tk_.StringVar(value=self.prefs['cusp_cutoff'])
        self.update_idletasks()
        cusp_frame.rowconfigure(8, weight=1)
        cusp_frame.columnconfigure(0, weight=1)
        cusp_frame.columnconfigure(3, weight=1)
        strut = ttk.Frame(cusp_frame, width=1)
        strut.grid(rowspan=8, column=0)
        next_label = ttk.Label(cusp_frame,
                               text='Which elements should be visible when you first '
                               'view the cusp neighborhood?')
        next_label.grid(row=0, column=1, columnspan=2, sticky=Tk_.W, pady=(20,10))
        next_check = ttk.Checkbutton(cusp_frame, variable = self.horoballs,
                                     text='Horoballs',
                                     command=self.set_horoballs)
        next_check.grid(row=1, column=1, sticky=Tk_.W, padx=(30, 0))
        next_check = ttk.Checkbutton(cusp_frame, variable = self.triangulation,
                                     text='Triangulation',
                                     command=self.set_triangulation)
        next_check.grid(row=2, column=1, sticky=Tk_.W, padx=(30, 0))
        next_check = ttk.Checkbutton(cusp_frame, variable = self.ford,
                                     text='Ford domain',
                                     command=self.set_ford)
        next_check.grid(row=3, column=1, sticky=Tk_.W, padx=(30, 0))
        next_check = ttk.Checkbutton(cusp_frame, variable = self.labels,
                                     text='Labels',
                                     command=self.set_labels)
        next_check.grid(row=4, column=1, sticky=Tk_.W, padx=(30, 0))
        next_check = ttk.Checkbutton(cusp_frame, variable = self.parallelogram,
                                     text='Parallelogram',
                                     command=self.set_parallelogram)
        next_check.grid(row=5, column=1, sticky=Tk_.W, padx=(30, 0))
        next_label = ttk.Label(cusp_frame,
                          text='What should the initial cutoff be?')
        next_label.grid(row=6, column=1, columnspan=2, pady=(20,10), sticky=Tk_.W)
        cutoff_entry = ttk.Entry(cusp_frame, textvariable=self.cutoff, width=15)
        cutoff_entry.grid(row=7, column=1, columnspan=2, sticky=Tk_.W,
                          pady=(0,10), padx=(30,0))

    def set_horoballs(self):
        self.prefs['cusp_horoballs'] = self.horoballs.get()

    def set_triangulation(self):
        self.prefs['cusp_triangulation'] = self.triangulation.get()

    def set_ford(self):
        self.prefs['cusp_ford_domain'] = self.ford.get()

    def set_labels(self):
        self.prefs['cusp_labels'] = self.labels.get()

    def set_parallelogram(self):
        self.prefs['cusp_parallelogram'] = self.parallelogram.get()

    def build_inside_pane(self, master):
        groupBG = self.style.groupBG
        self.keyboard = Tk_.Variable(value=self.prefs['keyboard'])
        self.inside_frame = frame = ttk.Frame(master)
        frame.rowconfigure(3, weight=1)
        frame.columnconfigure(0, weight=1)
        frame.columnconfigure(3, weight=1)
        # Keep the height the same as the height of the font pane.
        strut = ttk.Frame(frame, width=1)
        strut.grid(rowspan=5, column=0)
        keyboard_label = ttk.Label(frame, anchor=Tk_.W,
            text='Which keyboard layout are you using?')
        keyboard_label.grid(row=0, column=1, columnspan=2, sticky=Tk_.W, pady=(20,0))
        keyboard_button = ttk.OptionMenu(frame, self.keyboard, self.keyboard.get(),
            'QWERTY', 'AZERTY', 'QWERTZ', command=self.set_keyboard)
        keyboard_button.grid(row=1, column=1, columnspan=2, sticky=Tk_.W, pady=(10, 0))

    def set_keyboard(self, value):
        self.prefs['keyboard'] = value

if __name__ == '__main__':
    parent = Tk_.Tk(className='snappy')
    text_widget = Tk_.Text(parent)
    text_widget.insert(Tk_.INSERT, """
Lorem ipsum dolor sit amet, consectetur adipisicing
elit, sed do eiusmod tempor incididunt ut labore et
dolore magna aliqua. Ut enim ad minim veniam, quis
nostrud exercitation ullamco laboris nisi ut aliquip
ex ea commodo consequat. Duis aute irure dolor in
reprehenderit in voluptate velit esse cillum dolore
eu fugiat nulla pariatur. Excepteur sint occaecat
cupidatat non proident, sunt in culpa qui officia
deserunt mollit anim id est laborum.""")
    text_widget.tag_add('all', '1.0', Tk_.END)
    text_widget.tag_config('all', lmargin1=20, lmargin2=20)
    text_widget.pack()
    prefs = Preferences(text_widget)
    prefs.apply_prefs()
    PreferenceDialog(parent, prefs)
    parent.mainloop()
