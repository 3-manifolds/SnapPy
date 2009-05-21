import Tkinter as Tk_
import tkSimpleDialog, tkFont
from tkMessageBox import askyesno
import os, sys
from string import ascii_letters
try:
    import plistlib
except ImportError:
    import snappy.plistlib as plistlib

class Preferences:
    def __init__(self):
        self.prefs_dict = {'autocall' : False,
                           'automagic' : False,
                           'tracebacks' : False}
        if sys.platform == 'darwin':
            self.prefs_dict['font'] = ('Monaco', 16, 'normal')
        elif sys.platform == 'linux2':
            self.prefs_dict['font'] = ('fixed', 16, 'normal')
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

    def find_prefs(self):
        home = os.environ['HOME']
        if sys.platform == 'darwin':
            self.prefs_file = os.path.join(home, 'Library',
                                           'Preferences',
                                           'edu.t3m.SnapPy.plist')
        elif sys.platform == 'linux2':
            self.prefs_file = os.path.join(home, '.SnapPy.plist')
        else:
            self.prefs_file = None

    def read_prefs(self):
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

    def changed(self):
        return [key for key in self.cache.keys() if 
                self.cache[key] != self.prefs_dict[key]]

    # Override this in a subclass.
    def apply_prefs(self):
        print self.prefs_dict

class PreferenceDialog(tkSimpleDialog.Dialog):
    def __init__(self, parent, prefs, title='SnapPy Preferences'):
        Tk_.Toplevel.__init__(self, parent)
        self.parent = parent
        self.prefs = prefs
        self.prefs.cache_prefs()
        self.title(title)
        self.result = None
        self.build_font_panel()
        self.body_frame=self.font_frame
        self.build_shell_panel()
        tabs = [('Font', self.show_font_panel),
                ('Shell', self.show_shell_panel),
                ('This', self.show_shell_panel),
                ('That', self.show_shell_panel)]
        self.build_navbar(width=500, tabs=tabs)
        self.grab_set()
        self.protocol('WM_DELETE_WINDOW', self.cancel)
        self.buttonbox()
        self.font_button.invoke()
        self.wait_window(self)

    def ok_check(self):
        self.apply()
        answer = askyesno('Save?', 'Do you want to save these settings?')
        if answer:
            self.prefs.write_prefs()
        self.ok()

    def apply(self):
        self.prefs.apply_prefs()

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
        OK = Tk_.Button(box, text="OK", width=10, command=self.ok_check,
                        default=Tk_.ACTIVE)
        OK.pack(side=Tk_.LEFT, padx=5, pady=5)
        Apply = Tk_.Button(box, text="Apply", width=10, command=self.apply)
        Apply.pack(side=Tk_.LEFT, padx=5, pady=5)
        Cancel = Tk_.Button(box, text="Cancel", width=10, command=self.cancel)
        Cancel.pack(side=Tk_.LEFT, padx=5, pady=5)
        self.bind("<Return>", lambda event : OK.focus_set())
        self.bind("<Escape>", self.cancel)
        box.grid(row=2, column=0)

    def show_body(self):
        self.body_frame.grid(row=1, column=0, padx=5, pady=5,
                             sticky=Tk_.N + Tk_.S + Tk_.W + Tk_.E)
        self.body_frame.focus_set()

    def build_font_panel(self):
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
        current_family, current_size, current_style = self.prefs['font']
        current_family = families.index(current_family)
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
        sizer.insert(0, current_size, )
        sizer.grid(row=0, column=2, sticky=Tk_.W)
        label = Tk_.Label(font_frame, text='Weight: ')
        label.grid(row=1, column=1, sticky=Tk_.E)
        self.font_weight = weight = Tk_.StringVar(value='normal')
        if current_style.find('bold') > -1:
            weight.set('bold')
        radio = Tk_.Radiobutton(font_frame, text='normal',
                                variable=weight, value='normal',
                                command=self.set_font_sample)
        radio.grid(row=1, column=2, sticky=Tk_.W)
        radio = Tk_.Radiobutton(font_frame, text='bold',
                                variable=weight, value='bold',
                                command=self.set_font_sample)
        radio.grid(row=2, column=2, sticky=Tk_.W)
        label = Tk_.Label(font_frame, text='Style: ')
        label.grid(row=3, column=1, sticky=Tk_.E)
        self.font_style = style = Tk_.StringVar(value='roman')
        if current_style.find('italic') > -1:
            style.set('italic')
        radio = Tk_.Radiobutton(font_frame, text='roman',
                                variable=style, value='roman',
                                command=self.set_font_sample)
        radio.grid(row=3, column=2, sticky=Tk_.W)
        radio = Tk_.Radiobutton(font_frame, text='italic',
                                variable=style, value='italic',
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
        style = '%s %s'%(self.font_weight.get(), self.font_style.get())
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
        self.tracebacks = Tk_.BooleanVar(value=self.prefs['tracebacks'])
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
        next_check = Tk_.Checkbutton(shell_frame, variable = self.tracebacks,
                                     text='Show long tracebacks',
                                     command=self.set_tracebacks)
        next_check.grid(row=3, column=1, sticky=Tk_.W, pady=5)

    def set_autocall(self):
        self.prefs['autocall'] = self.autocall.get()

    def set_automagic(self):
        self.prefs['automagic'] = self.automagic.get()

    def set_tracebacks(self):
        self.prefs['tracebacks'] = self.tracebacks.get()

    def show_shell_panel(self):
        self.body_frame.grid_remove()
        self.body_frame = self.shell_frame
        self.show_body()


if __name__ == '__main__':
    parent = Tk_.Tk()
    parent.withdraw()
    prefs = Preferences()
    PreferenceDialog(parent, prefs)
    parent.deiconify()
    parent.mainloop()
