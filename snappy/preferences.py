import Tkinter as Tk_
import tkSimpleDialog
import tkFont
from string import ascii_letters

class PreferenceDialog(tkSimpleDialog.Dialog):
    def __init__(self, parent, title='Preferences',
                 prefs={'font': ('Monaco', 18, 'normal')}):
        Tk_.Toplevel.__init__(self, parent)
        self.title(title)
        self.parent = parent
        self.result = None
        self.prefs = prefs
        self.navbar()
        self.build_font_panel()
        self.grab_set()
        self.protocol('WM_DELETE_WINDOW', self.cancel)
        self.body_frame = self.font_frame
        self.show_font_panel()
        self.buttonbox()
        self.build_shell_panel()
        self.wait_window(self)
    
    def navbar(self):
        box = Tk_.Frame(self)
        var = Tk_.IntVar()
        self.font_button = Tk_.Radiobutton(box, text="Font", width=10,
                                      command=self.show_font_panel,
                                      variable=var, value=0, indicatoron=0)
        self.font_button.pack(side=Tk_.LEFT, padx=1, pady=5)
        self.shell_button = Tk_.Radiobutton(box, text="Shell", width=10,
                                       command=self.show_shell_panel,
                                       variable=var, value=1, indicatoron=0)
        self.shell_button.pack(side=Tk_.LEFT, padx=1, pady=5)
        self.Xbutton = Tk_.Radiobutton(box, text="Three", width=10,
                                       command=self.show_shell_panel,
                                       variable=var, value=2, indicatoron=0)
        self.Xbutton.pack(side=Tk_.LEFT, padx=1, pady=5)
        self.Ybutton = Tk_.Radiobutton(box, text="Four", width=10,
                                       command=self.show_shell_panel,
                                       variable=var, value=3, indicatoron=0)
        self.Ybutton.pack(side=Tk_.LEFT, padx=1, pady=5)
        box.grid(row=0, column=0)

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
        self.body_frame.grid(row=1, column=0, padx=5, pady=5)
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
                                        width=50, height=4,
                                        highlightthickness=0,
                                        relief=Tk_.FLAT,
                                        font=self.prefs['font'])
        self.sample.bind('<Button-1>', lambda event: 'break')
        self.sample.insert(Tk_.INSERT, 'ABCDEFGHIJKLMNOPQRSTUVWXYZ\n'\
                                           'abcdefghijklmnopqrstuvwxyz')
        sample.tag_add('all', '1.0', Tk_.END)
        sample.tag_config('all', justify=Tk_.CENTER,
                          font=self.prefs['font'])
        font_list.bind('<ButtonRelease-1>', self.set_font_sample)
        self.font_list.grid(row=0, column=0)
        self.sample.grid(row=1, column=0, pady=10)
        self.list_frame.grid(row=0, column=0)
        self.font_button.select()

    def show_font_panel(self):
        self.body_frame.grid_forget()
        self.body_frame = self.font_frame
        self.show_body()

    def set_font_sample(self, event):
        index = self.font_list.curselection()[0]
        family = self.families[int(index)]
        size = 18
        weight = 'normal'
        self.sample.tag_config('all', justify=Tk_.CENTER,
                               font=(family, size, weight)) 

    def build_shell_panel(self):
        self.update_idletasks()
        self.shell_frame = shell_frame = Tk_.Frame(self,
                     width=self.font_frame.winfo_reqwidth(),
                     height=self.font_frame.winfo_reqheight())
        
    def show_shell_panel(self):
        self.body_frame.grid_remove()
        self.body_frame = self.shell_frame
        self.show_body()

if __name__ == '__main__':
    parent = Tk_.Tk()
    PreferenceDialog(parent)
    parent.mainloop()
