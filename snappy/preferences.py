import Tkinter as Tk_
import tkSimpleDialog
import tkFont
from string import ascii_letters

class PreferenceDialog(tkSimpleDialog.Dialog):
    def __init__(self, parent, title='Preferences'):
        Tk_.Toplevel.__init__(self, parent)
        self.pref_dict = {'font' : ('sanserif', 16, 'normal'),
                          'autocall' : True,
                          'automagic' : False,
                          'tracebacks' : False}
        self.title(title)
        self.parent = parent
        self.result = None
        self.build_font_panel()
        self.body_frame=self.font_frame
        self.build_shell_panel()
        tabs = [('Font', self.show_font_panel),
                ('Shell', self.show_shell_panel),
                ('Three', self.show_shell_panel),
                ('Four', self.show_shell_panel)]
        self.build_navbar(width=500, tabs=tabs)
        self.grab_set()
        self.protocol('WM_DELETE_WINDOW', self.cancel)
        self.body_frame = self.font_frame
        self.show_font_panel()
        self.buttonbox()
        self.wait_window(self)
    
    def build_navbar(self, width=500, tabs=[]):
        navbox = Tk_.Frame(self)
        navbox.columnconfigure(0, weight=1)
        navbox.columnconfigure(len(tabs)+1, weight=1)
        var = Tk_.StringVar()
        for n in range(len(tabs)-1, -1, -1):
            tabtext, tabfunc = tabs[n]
            button = Tk_.Radiobutton(navbox, text=tabtext, width=10,
                                     command=tabfunc, variable=var,
                                     value=tabtext, indicatoron=0)
            button.grid(row=0, column=n, padx=0, pady=5, sticky=Tk_.E)
        strut=Tk_.Frame(navbox, width=width, bg='Black')
        strut.grid(row=1, columnspan=6)
        navbox.grid(row=0, column=0, pady=10)
        button.select()  # This does not work.  It is a Tkinter bug
                         # that arises when indicatoron=0.
        tabfunc()

    def buttonbox(self):
        box = Tk_.Frame(self)
        OK = Tk_.Button(box, text="OK", width=10, command=self.ok,
                        default=Tk_.ACTIVE)
        OK.pack(side=Tk_.LEFT, padx=5, pady=5)
        Apply = Tk_.Button(box, text="Apply", width=10, command=self.apply)
        Apply.pack(side=Tk_.LEFT, padx=5, pady=5)
        Cancel = Tk_.Button(box, text="Cancel", width=10, command=self.cancel)
        Cancel.pack(side=Tk_.LEFT, padx=5, pady=5)
        self.bind("<Return>", self.ok)
        self.bind("<Escape>", self.cancel)
        box.grid(row=2, column=0)

    def show_body(self):
        self.body_frame.grid(row=1, column=0, padx=5, pady=5,
                             sticky=Tk_.N + Tk_.S + Tk_.W + Tk_.E)
        self.body_frame.focus_set()

    def build_font_panel(self):
        self.font_frame = font_frame = Tk_.Frame(self)
        self.list_frame = list_frame = Tk_.Frame(font_frame)
        self.font_list = font_list = Tk_.Listbox(list_frame,
                                                 selectmode=Tk_.SINGLE)
        self.families = families = [x for x in tkFont.families()
                                    if x[0] in ascii_letters]
        families.sort()
        for family in families:
            font_list.insert(Tk_.END, family)
        self.font_scroller = font_scroller = Tk_.Scrollbar(list_frame,
                                                           command=font_list.yview)
        font_list.config(yscrollcommand=font_scroller.set)
        self.font_scroller.grid(row=0, column=1, sticky=Tk_.N + Tk_.S)
        self.sample = sample = Tk_.Text(self.font_frame,
                                        width=40, height=4,
                                        highlightthickness=0,
                                        relief=Tk_.RIDGE,
                                        font=self.pref_dict['font'])
        self.sample.bind('<Button-1>', lambda event: 'break')
        self.sample.insert(Tk_.INSERT, 'ABCDEFGHIJKLMNOPQRSTUVWXYZ\n'\
                                           'abcdefghijklmnopqrstuvwxyz')
        sample.tag_add('all', '1.0', Tk_.END)
        sample.tag_config('all', justify=Tk_.CENTER,
                          font=self.pref_dict['font'])
        font_list.bind('<ButtonRelease-1>', self.set_font_sample)
        self.font_list.grid(row=0, column=0)
        self.sample.grid(row=1, column=0, pady=10, sticky=Tk_.E+Tk_.W)
        self.list_frame.grid(row=0, column=0)

    def set_font_sample(self, event):
        index = self.font_list.curselection()[0]
        family = self.families[int(index)]
        size = 18
        weight = 'normal'
        self.sample.tag_config('all', justify=Tk_.CENTER,
                               font=(family, size, weight)) 

    def show_font_panel(self):
        self.body_frame.grid_forget()
        self.body_frame = self.font_frame
        self.show_body()

    def build_shell_panel(self):
        self.autocall = Tk_.BooleanVar(value=self.pref_dict['autocall'])
        self.automagic = Tk_.BooleanVar(value=self.pref_dict['automagic'])
        self.tracebacks = Tk_.BooleanVar(value=self.pref_dict['tracebacks'])
        self.update_idletasks()
        self.shell_frame = shell_frame = Tk_.Frame(self)
        shell_frame.rowconfigure(0, weight=1)
        shell_frame.rowconfigure(4, weight=1)
        shell_frame.columnconfigure(0, weight=1)
        shell_frame.columnconfigure(3, weight=1)
        strut = Tk_.Frame(shell_frame, width=1,
                             height=self.font_frame.winfo_reqheight())
        strut.grid(rowspan=5, column=0)
        next_check = Tk_.Checkbutton(shell_frame, variable = self.autocall,
                                     text='IPython autocall')
        next_check.grid(row=1, column=1, sticky=Tk_.W, pady=5)
        next_check = Tk_.Checkbutton(shell_frame, variable = self.automagic,
                                  text='IPython automagic')
        next_check.grid(row=2, column=1, sticky=Tk_.W, pady=5)
        next_check = Tk_.Checkbutton(shell_frame, variable = self.tracebacks,
                                  text='Show long tracebacks')
        next_check.grid(row=3, column=1, sticky=Tk_.W, pady=5)

    def show_shell_panel(self):
        self.body_frame.grid_remove()
        self.body_frame = self.shell_frame
        self.show_body()

if __name__ == '__main__':
    parent = Tk_.Tk()
    PreferenceDialog(parent)
    parent.mainloop()
