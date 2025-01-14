import os
import sys
import plistlib
from .gui import *
from .app_menus import ListedWindow
from .tkterminal import FontChoice, default_terminal_font

def clip_font_size(size):
    try:
        value = int(size)
    except:
        return 13 if sys.platform == 'darwin' else 11
    if value < 10:
        return 10
    if value > 36:
        return 36
    return size

class Settings:
    def __init__(self):
        self.setting_dict = {
            'autocall' : False,
            'font' : self.default_font(),
            'cusp_horoballs' : True,
            'cusp_triangulation' : True,
            'cusp_ford_domain' : True,
            'cusp_parallelogram' : True,
            'cusp_labels' : False,
            'cusp_cutoff' : '0.1000',
            'keyboard' : 'QWERTY'}
        self.cache = {}
        self.cache_settings()
        self.find_settings()
        self.read_settings()

    def __getitem__(self, x):
        return self.setting_dict[x]

    def __setitem__(self, x, y):
        self.setting_dict[x] = y

    def __repr__(self):
        return str(self.setting_dict)

    def items(self):
        return self.setting_dict.items()

    def get(self, key, default):
        return self.setting_dict.get(key, default)

    def default_font(self):
        return default_terminal_font()

    def find_settings(self):
        if sys.platform == 'darwin':
            home = os.environ['HOME']
            # File is 'Preferences' rather than 'Settings' for backwards
            # compatibility.
            self.setting_file = os.path.join(home, 'Library',
                                           'Preferences',
                                           'org.computop.SnapPy.plist')
        elif sys.platform == 'linux2' or sys.platform == 'linux':
            home = os.environ['HOME']
            self.setting_file = os.path.join(home, '.SnapPy.plist')
        elif sys.platform == 'win32':
            home = os.environ['USERPROFILE']
            self.setting_file = os.path.join(home, '.SnapPy.plist')
        else:
            self.setting_file = None
        
    def read_settings(self):
        if self.setting_file:
            try:
                if hasattr(plistlib, 'load'):
                    with open(self.setting_file, 'rb') as setting_file:
                        self.setting_dict.update(plistlib.load(setting_file))
                else:
                    self.setting_dict.update(plistlib.readPlist(self.setting_file))
                family, size, info = self.setting_dict['font']
                # Guard against crazy values in the settings plist file.
                size = clip_font_size(size)
                weight, slant = info.split()
                self.setting_dict['font'] = FontChoice(family, size, weight, slant)
            except OSError:
                pass

    def write_settings(self):
        if self.setting_file:
            setting_dict = self.setting_dict.copy()
            font = self.setting_dict['font']
            setting_dict['font'] = (font.family, font.size, font.rest)
            if hasattr(plistlib, 'dump'):
                with open(self.setting_file, 'wb') as setting_file:
                    plistlib.dump(setting_dict, setting_file)
            else:
                plistlib.writePlist(setting_dict, self.setting_file)

    def cache_settings(self):
        self.cache.update(self.setting_dict)

    def revert_settings(self):
        self.setting_dict.update(self.cache)
        self.apply_settings()

    def changed(self):
        return [key for key in self.cache.keys() if
                self.cache[key] != self.setting_dict[key]]

    # Override this in a subclass.
    def apply_settings(self):
        print(self.setting_dict)


class SettingsDialog(Dialog):
    def __init__(self, parent, settings, title='SnapPy Settings'):
        self.parent = parent
        self.style = SnapPyStyle()
        self.settings = settings
        self.settings.cache_settings()
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
        self.attributes('-topmost', True)

    def run(self):
        self.wait_window(self)

    def validate(self):
        cutoff = self.cutoff.get()
        try:
            float(cutoff)
            self.settings['cusp_cutoff'] = cutoff
        except ValueError:
            showerror('Invalid input',
                      'Please enter a number for the cutoff.')
            return False
        self.settings.apply_settings()
        self.okay = True
        return True

    def revert(self):
        self.settings.revert_settings()
        self.cancel()

    def buttonbox(self):
        self.button_frame = box = ttk.Frame(self, padding=(0, 0, 0, 20))
        box.grid_columnconfigure(0, weight=1)
        box.grid_columnconfigure(2, weight=1)
        OK = ttk.Button(box, text="OK", width=10, command=self.ok,
                        default=Tk_.ACTIVE)
        OK.grid(row=0, column=0, sticky=Tk_.NE, padx=5)
        Apply = ttk.Button(box, text="Apply", width=10,
                           command=self.settings.apply_settings)
        Apply.grid(row=0, column=1, sticky=Tk_.N, padx=5)
        Cancel = ttk.Button(box, text="Cancel", width=10, command=self.revert)
        Cancel.grid(row=0, column=2, sticky=Tk_.NW, padx=5)
        self.bind("<Return>", lambda event : OK.focus_set())
        self.bind("<Escape>", self.cancel)

    def show_body(self):
        self.body_frame.grid(row=1, column=0, padx=5, pady=5,
                             sticky=Tk_.N + Tk_.S + Tk_.W + Tk_.E)
        self.body_frame.focus_set()

    def build_font_pane(self, parent):
        current_font = self.settings['font']
        groupBG = self.style.groupBG
        self.font_frame = font_frame = ttk.Frame(parent)
        font_frame.columnconfigure(2, weight=1)
        font_frame.columnconfigure(0, weight=1)
        self.list_frame = list_frame = ttk.Frame(font_frame)
        self.font_list = font_list = Tk_.Listbox(list_frame,
                                                 selectmode=Tk_.SINGLE)
        font_list.grid(row=0, column=0, pady=(20,30))
        font_scroller = ttk.Scrollbar(list_frame,
                                      command=font_list.yview)
        font_list.config(yscrollcommand=font_scroller.set)
        font_scroller.grid(row=0, column=1, sticky=Tk_.N + Tk_.S, pady=(20,30))
        list_frame.grid(rowspan=6, column=0)

        label = ttk.Label(self.font_frame, text='Size: ')
        label.grid(row=0, column=1, sticky=Tk_.E, pady=(20,0))
        self.font_sizer = sizer = Spinbox(font_frame, from_=10, to_=36, width=4,
                                          command=self.set_font_sample)
        def check_font_size(size):
            try:
                value = int(size)
            except:
                sizer.set(13 if sys.platform == 'darwin' else 11)
                return False
            if value < 10:
                sizer.set(10)
                return False
            if value > 36:
                sizer.set(36)
                return False
            return True
        validator = parent.register(check_font_size)
        sizer.config(validate ="key", validatecommand =(validator, '%P'))
        sizer.delete(0,2)
        sizer.insert(0, clip_font_size(current_font.size))
        sizer.bind('<Return>', self.set_font_sample)
        sizer.bind('<Tab>', self.set_font_sample)
        self.current_size = current_font.size
        sizer.grid(row=0, column=2, sticky=Tk_.W, pady=(20,0))
        label = ttk.Label(font_frame, text='Weight: ')
        label.grid(row=1, column=1, sticky=Tk_.E)
        self.font_weight = weight = Tk_.StringVar(value=current_font.weight)
        radio = ttk.Radiobutton(font_frame, text='normal', variable=weight,
                                    value='normal', command=self.set_font_sample)
        radio.grid(row=1, column=2, sticky=Tk_.W)
        radio = ttk.Radiobutton(font_frame, text='bold', variable=weight,
                                    value='bold', command=self.set_font_sample)
        radio.grid(row=2, column=2, sticky=Tk_.W)
        label = ttk.Label(font_frame, text='Slant: ')
        label.grid(row=3, column=1, sticky=Tk_.E)
        self.font_slant = slant = Tk_.StringVar(value=current_font.slant)
        radio = ttk.Radiobutton(font_frame, text='roman', variable=slant,
                                 value='roman', command=self.set_font_sample)
        radio.grid(row=3, column=2, sticky=Tk_.W)
        radio = ttk.Radiobutton(font_frame, text='italic', variable=slant,
                                    value='italic', command=self.set_font_sample)
        radio.grid(row=4, column=2, sticky=Tk_.W)

        self.fixed_only =  Tk_.BooleanVar(value=True)
        self.check_fixed = ttk.Checkbutton(font_frame, variable=self.fixed_only,
                                           text='Fixed-width fonts only',
                                           command=self.set_font_families)
        self.check_fixed.grid(row=6, column=0, pady=(0, 10), sticky=Tk_.N)
        self.sample = sample = Tk_.Text(font_frame,
                                        width=40, height=6,
                                        highlightthickness=0,
                                        relief=Tk_.FLAT,
                                        font=self.settings['font'].as_tuple())
        self.sample.bind('<Button-1>', lambda event: 'break')
        self.sample.insert(Tk_.INSERT, '\nABCDEFGHIJKLMNOPQRSTUVWXYZ\n'
                           'abcdefghijklmnopqrstuvwxyz')
        sample.tag_add('all', '1.0', Tk_.END)
        font_list.bind('<ButtonRelease-1>', self.set_font_sample)
        self.sample.grid(row=7, columnspan=4, padx=10, pady=10, sticky=Tk_.E+Tk_.W)
        self.set_font_families()
        self.set_font_sample()

    def set_font_families(self):
        families = {f for f in font_families()
                    if f[0] != '@' and f.lower().find('emoji') < 0} # omit vertical fonts, emojis
        if self.fixed_only.get():
            families = {f for f in families if Font(family=f).metrics('fixed')}
        families.add(self.settings['font'].family)
        self.families = families = sorted(families)
        if self.font_list.size() > 0:
            self.font_list.delete(0, self.font_list.size() - 1)
        for family in families:
            self.font_list.insert(Tk_.END, family)
        self.current_family = families.index(self.settings['font'].family)
        self.font_list.selection_set(self.current_family)
        self.font_list.see(self.current_family)

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
        return FontChoice(family, size,
                          self.font_weight.get(), self.font_slant.get())

    def set_font_sample(self, event=None):
        new_font = self.get_font()
        self.settings['font'] = new_font
        self.sample.tag_config('all', justify=Tk_.CENTER,
                               font=new_font.as_tuple())

    def build_shell_pane(self, parent):
        groupBG = self.style.groupBG
        self.autocall = Tk_.BooleanVar(value=self.settings['autocall'])
        self.update_idletasks()
        self.shell_frame = shell_frame = ttk.Frame(parent)
        shell_frame.rowconfigure(3, weight=1)
        shell_frame.columnconfigure(0, weight=1)
        shell_frame.columnconfigure(3, weight=1)
        # Keep the height the same as the height of the font pane.
        strut = ttk.Frame(shell_frame, width=1)
        strut.grid(rowspan=5, column=0)
        next_label = ttk.Label(shell_frame, anchor=Tk_.W,
            text='Which IPython features would you like to enable?')
        next_label.grid(row=0, column=1, columnspan=2, sticky=Tk_.W, pady=(20,0))
        next_check = ttk.Checkbutton(shell_frame, variable=self.autocall,
                                    text='IPython autocall',
                                    command=self.set_autocall)
        next_check.grid(row=1, column=1, sticky=Tk_.W, pady=(10,0))
        next_check.grid(row=2, column=1, sticky=Tk_.W, pady=(5,0))

    def set_autocall(self):
        self.settings['autocall'] = self.autocall.get()

    def build_cusp_pane(self, parent):
        groupBG = self.style.groupBG
        self.cusp_frame = cusp_frame = ttk.Frame(parent)
        self.horoballs = Tk_.BooleanVar(value=self.settings['cusp_horoballs'])
        self.triangulation = Tk_.BooleanVar(value=self.settings['cusp_triangulation'])
        self.ford = Tk_.BooleanVar(value=self.settings['cusp_ford_domain'])
        self.labels = Tk_.BooleanVar(value=self.settings['cusp_labels'])
        self.parallelogram = Tk_.BooleanVar(value=self.settings['cusp_parallelogram'])
        self.cutoff = Tk_.StringVar(value=self.settings['cusp_cutoff'])
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
        next_check = ttk.Checkbutton(cusp_frame, variable=self.horoballs,
                                     text='Horoballs',
                                     command=self.set_horoballs)
        next_check.grid(row=1, column=1, sticky=Tk_.W, padx=(30, 0))
        next_check = ttk.Checkbutton(cusp_frame, variable=self.triangulation,
                                     text='Triangulation',
                                     command=self.set_triangulation)
        next_check.grid(row=2, column=1, sticky=Tk_.W, padx=(30, 0))
        next_check = ttk.Checkbutton(cusp_frame, variable=self.ford,
                                     text='Ford domain',
                                     command=self.set_ford)
        next_check.grid(row=3, column=1, sticky=Tk_.W, padx=(30, 0))
        next_check = ttk.Checkbutton(cusp_frame, variable=self.labels,
                                     text='Labels',
                                     command=self.set_labels)
        next_check.grid(row=4, column=1, sticky=Tk_.W, padx=(30, 0))
        next_check = ttk.Checkbutton(cusp_frame, variable=self.parallelogram,
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
        self.settings['cusp_horoballs'] = self.horoballs.get()

    def set_triangulation(self):
        self.settings['cusp_triangulation'] = self.triangulation.get()

    def set_ford(self):
        self.settings['cusp_ford_domain'] = self.ford.get()

    def set_labels(self):
        self.settings['cusp_labels'] = self.labels.get()

    def set_parallelogram(self):
        self.settings['cusp_parallelogram'] = self.parallelogram.get()

    def build_inside_pane(self, parent):
        groupBG = self.style.groupBG
        self.keyboard = Tk_.Variable(value=self.settings['keyboard'])
        self.inside_frame = frame = ttk.Frame(parent)
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
        self.settings['keyboard'] = value


if __name__ == '__main__':
    parent = Tk_.Tk(className='snappy')
    settings = Settings()
    settings.apply_settings()
    SettingsDialog(parent, settings)
    parent.mainloop()
