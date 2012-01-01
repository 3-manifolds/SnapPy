try:
    import Tkinter as Tk_
    import tkSimpleDialog, tkFont
    from tkMessageBox import showerror
except ImportError:
    import tkinter as Tk_
    import tkinter.simpledialog as tkSimpleDialog
    import tkinter.font as tkFont
    from tkinter.messagebox import showerror

import os, sys, re, time
from string import ascii_letters
try:
    import plistlib
except ImportError:
    from . import plistlib

class Preferences:
    def __init__(self, text_widget):
        self.text_widget = text_widget
        self.prefs_dict = {'autocall' : False,
                           'automagic' : False,
                           'font' : self.current_font_tuple(),
                           'cusp_horoballs' : True,
                           'cusp_triangulation' : True,
                           'cusp_ford_domain' : True,
                           'cusp_parallelogram' : True,
                           'cusp_labels' : True,
                           'cusp_cutoff' : '0.1000'}
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

    def current_font_dict(self):
        font_string = self.text_widget.cget('font')
        return tkFont.Font(font=font_string).actual()

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
        elif sys.platform == 'linux2':
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

class PreferenceDialog(tkSimpleDialog.Dialog):
    def __init__(self, parent, prefs, title='SnapPy Preferences'):
        self.parent = parent
        self.prefs = prefs
        self.prefs.cache_prefs()
        self.okay = False
        Tk_.Toplevel.__init__(self, parent)
#        if parent.winfo_viewable():
#            self.transient(parent)
        self.title(title)
        self.build_font_panel()
        self.body_frame=self.font_frame
        self.build_shell_panel()
        self.build_cusp_panel()
        tabs = [('Font', self.show_font_panel),
                ('Shell', self.show_shell_panel),
                ('Cusps', self.show_cusp_panel)]
        self.build_navbar(width=500, tabs=tabs)
        self.buttonbox()
        self.font_button.invoke()
        self.initial_focus = self.body_frame
        self.initial_focus.focus_set()
#        self.wait_visibility(self)
#        self.grab_set()
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

    def build_navbar(self, width=500, tabs=[]):
        navbox = Tk_.Frame(self)
        navbox.columnconfigure(0, weight=1)
        navbox.columnconfigure(len(tabs)+1, weight=1)
        selectedButton = Tk_.StringVar()
        for n in range(len(tabs)-1, -1, -1):
            tabtext, tabfunc = tabs[n]
            button = Tk_.Radiobutton(navbox, text=tabtext, width=10,
                                     command=tabfunc, variable=selectedButton,
                                     value=tabtext, indicatoron=0)
            button.grid(row=0, column=n, padx=0, pady=5, sticky=Tk_.E)
        # If nothing is packed into a frame, it keeps its initial size.
        strut=Tk_.Frame(navbox, width=width, bg='Black')
        strut.grid(row=1, columnspan=6)
        navbox.grid(row=0, column=0, pady=10)
        self.font_button = button

    def buttonbox(self):
        box = Tk_.Frame(self)
        OK = Tk_.Button(box, text="OK", width=10, command=self.ok,
                        default=Tk_.ACTIVE)
        OK.pack(side=Tk_.LEFT, padx=5, pady=5)
        Apply = Tk_.Button(box, text="Apply", width=10,
                           command=self.prefs.apply_prefs)
        Apply.pack(side=Tk_.LEFT, padx=5, pady=5)
        Cancel = Tk_.Button(box, text="Cancel", width=10, command=self.revert)
        Cancel.pack(side=Tk_.LEFT, padx=5, pady=5)
        self.bind("<Return>", lambda event : OK.focus_set())
        self.bind("<Escape>", self.cancel)
        box.grid(row=2, column=0)

    def show_body(self):
        self.body_frame.grid(row=1, column=0, padx=5, pady=5,
                             sticky=Tk_.N + Tk_.S + Tk_.W + Tk_.E)
        self.body_frame.focus_set()

    def build_font_panel(self):
        current_font = self.prefs.current_font_dict()
        self.font_frame = font_frame = Tk_.Frame(self)
        font_frame.columnconfigure(2, weight=1)
        font_frame.columnconfigure(0, weight=1)
        self.list_frame = list_frame = Tk_.Frame(font_frame)
        self.font_list = font_list = Tk_.Listbox(list_frame,
                                                 selectmode=Tk_.SINGLE)
        self.families = families = [x for x in tkFont.families()
                                    if x[0] in ascii_letters]
        families.sort()
        for family in families:
            font_list.insert(Tk_.END, family)
        current_family = families.index(current_font['family'])
        font_list.selection_set(current_family)
        font_list.see(current_family)
        font_list.grid(row=0, column=0)
        font_scroller = Tk_.Scrollbar(list_frame,
                                      command=font_list.yview)
        font_list.config(yscrollcommand=font_scroller.set)
        font_scroller.grid(row=0, column=1, sticky=Tk_.N + Tk_.S)
        list_frame.grid(rowspan=6, column=0)
        
        label = Tk_.Label(self.font_frame, text='Size: ')
        label.grid(row=0, column=1, sticky=Tk_.E)
        self.font_sizer = sizer = Tk_.Spinbox(font_frame,
                                              from_=10, to_=36,
                                              width=4,
                                              font='Helvetica 14',
                                              command=self.set_font_sample)
        sizer.delete(0,2)
        sizer.insert(0, current_font['size'] )
        sizer.grid(row=0, column=2, sticky=Tk_.W)
        label = Tk_.Label(font_frame, text='Weight: ')
        label.grid(row=1, column=1, sticky=Tk_.E)
        self.font_weight = weight = Tk_.StringVar(value=current_font['weight'])
        radio = Tk_.Radiobutton(font_frame, text='normal',
                                variable=weight, value='normal',
                                command=self.set_font_sample)
        radio.grid(row=1, column=2, sticky=Tk_.W)
        radio = Tk_.Radiobutton(font_frame, text='bold',
                                variable=weight, value='bold',
                                command=self.set_font_sample)
        radio.grid(row=2, column=2, sticky=Tk_.W)
        label = Tk_.Label(font_frame, text='Slant: ')
        label.grid(row=3, column=1, sticky=Tk_.E)
        self.font_slant = slant = Tk_.StringVar(value=current_font['slant'])
        radio = Tk_.Radiobutton(font_frame, text='roman',
                                variable=slant, value='roman',
                                command=self.set_font_sample)
        radio.grid(row=3, column=2, sticky=Tk_.W)
        radio = Tk_.Radiobutton(font_frame, text='italic',
                                variable=slant, value='italic',
                                command=self.set_font_sample)
        radio.grid(row=4, column=2, sticky=Tk_.W)
        self.sample = sample = Tk_.Text(self.font_frame,
                                        width=40, height=6,
                                        highlightthickness=0,
                                        relief=Tk_.RIDGE,
                                        font=self.prefs['font'])
        self.sample.bind('<Button-1>', lambda event: 'break')
        self.sample.insert(Tk_.INSERT, 'ABCDEFGHIJKLMNOPQRSTUVWXYZ\n'\
                                           'abcdefghijklmnopqrstuvwxyz')
        sample.tag_add('all', '1.0', Tk_.END)
        self.set_font_sample()
        font_list.bind('<ButtonRelease-1>', self.set_font_sample)
        self.sample.grid(row=6, columnspan=4, pady=10, sticky=Tk_.E+Tk_.W)

    def get_font(self):
        index = self.font_list.curselection()[0]
        family = self.families[int(index)]
        size = int(self.font_sizer.get())
        style = '%s %s'%(self.font_weight.get(), self.font_slant.get())
        return (family, size, style.strip()) 

    def set_font_sample(self, event=None):
        new_font = self.get_font()
        self.prefs['font'] = new_font
        self.sample.tag_config('all', justify=Tk_.CENTER,
                               font=new_font) 

    def show_font_panel(self):
        self.body_frame.grid_forget()
        self.body_frame = self.font_frame
        self.show_body()

    def build_shell_panel(self):
        self.autocall = Tk_.BooleanVar(value=self.prefs['autocall'])
        self.automagic = Tk_.BooleanVar(value=self.prefs['automagic'])
        self.update_idletasks()
        self.shell_frame = shell_frame = Tk_.Frame(self)
        shell_frame.rowconfigure(0, weight=1)
        shell_frame.rowconfigure(4, weight=1)
        shell_frame.columnconfigure(0, weight=1)
        shell_frame.columnconfigure(3, weight=1)
        # Keep the height the same as the height of the font panel.
        strut = Tk_.Frame(shell_frame, width=1,
                             height=self.font_frame.winfo_reqheight())
        strut.grid(rowspan=5, column=0)
        next_check = Tk_.Checkbutton(shell_frame, variable = self.autocall,
                                     text='IPython autocall',
                                     command=self.set_autocall)
        next_check.grid(row=1, column=1, sticky=Tk_.W, pady=5)
        next_check = Tk_.Checkbutton(shell_frame, variable = self.automagic,
                                     text='IPython automagic',
                                     command=self.set_automagic)
        next_check.grid(row=2, column=1, sticky=Tk_.W, pady=5)
        next_check.grid(row=3, column=1, sticky=Tk_.W, pady=5)

    def set_autocall(self):
        self.prefs['autocall'] = self.autocall.get()

    def set_automagic(self):
        self.prefs['automagic'] = self.automagic.get()

    def show_shell_panel(self):
        self.body_frame.grid_remove()
        self.body_frame = self.shell_frame
        self.show_body()

    def build_cusp_panel(self):
        self.horoballs = Tk_.BooleanVar(value=self.prefs['cusp_horoballs'])
        self.triangulation = Tk_.BooleanVar(value=self.prefs['cusp_triangulation'])
        self.ford = Tk_.BooleanVar(value=self.prefs['cusp_ford_domain'])
        self.labels = Tk_.BooleanVar(value=self.prefs['cusp_labels'])
        self.parallelogram = Tk_.BooleanVar(value=self.prefs['cusp_parallelogram'])
        self.cutoff = Tk_.StringVar(value=self.prefs['cusp_cutoff'])
        self.update_idletasks()
        self.cusp_frame = cusp_frame = Tk_.Frame(self)
        cusp_frame.rowconfigure(8, weight=1)
        cusp_frame.columnconfigure(0, weight=1)
        cusp_frame.columnconfigure(3, weight=1)
        # Keep the height the same as the height of the font panel.
        strut = Tk_.Frame(cusp_frame, width=1,
                             height=self.font_frame.winfo_reqheight())
        strut.grid(rowspan=8, column=0)
        next_label = Tk_.Label(cusp_frame,
                          text='Which elements should be visible when you first\n'
                               'view the cusp neighborhood?')
        next_label.grid(row=0, column=1, columnspan=2, sticky=Tk_.N)
        next_check = Tk_.Checkbutton(cusp_frame, variable = self.horoballs,
                                     text='Horoballs',
                                     command=self.set_horoballs)
        next_check.grid(row=1, column=1, sticky=Tk_.W)
        next_check = Tk_.Checkbutton(cusp_frame, variable = self.triangulation,
                                     text='Triangulation',
                                     command=self.set_triangulation)
        next_check.grid(row=2, column=1, sticky=Tk_.W)
        next_check = Tk_.Checkbutton(cusp_frame, variable = self.ford,
                                     text='Ford domain',
                                     command=self.set_ford)
        next_check.grid(row=3, column=1, sticky=Tk_.W)
        next_check = Tk_.Checkbutton(cusp_frame, variable = self.labels,
                                     text='Labels',
                                     command=self.set_labels)
        next_check.grid(row=4, column=1, sticky=Tk_.W)
        next_check = Tk_.Checkbutton(cusp_frame, variable = self.parallelogram,
                                     text='Parallelogram',
                                     command=self.set_parallelogram)
        next_check.grid(row=5, column=1, sticky=Tk_.W)
        next_label = Tk_.Label(cusp_frame,
                          text='What should the initial cutoff be?')
        next_label.grid(row=6, columnspan=2, pady=(10,0))
        cutoff_entry = Tk_.Entry(cusp_frame, textvariable=self.cutoff)
        cutoff_entry.grid(row=7, column=1, columnspan=2, sticky=Tk_.NW, pady=(0,10))

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

    def show_cusp_panel(self):
        self.body_frame.grid_remove()
        self.body_frame = self.cusp_frame
        self.show_body()


if __name__ == '__main__':
    parent = Tk_.Tk()
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
