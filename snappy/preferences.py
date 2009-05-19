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
        self.body_frame = body_frame = Tk_.Frame(self)
        self.initial_focus = self.font_panel(body_frame)
        body_frame.pack(padx=5, pady=5)
        self.buttonbox()
        self.grab_set()
        if not self.initial_focus:
            self.initial_focus = self
        self.protocol('WM_DELETE_WINDOW', self.cancel)
        self.initial_focus.focus_set()
        self.wait_window(self)
    
    def navbar(self):
        box = Tk_.Frame(self)
        Font = Tk_.Button(box, text="Font", width=10,
                          command=self.font_panel, default=Tk_.ACTIVE)
        Font.pack(side=Tk_.LEFT, padx=1, pady=5)
        button = Tk_.Button(box, text="Two", width=10)
        button.pack(side=Tk_.LEFT, padx=1, pady=5)
        button = Tk_.Button(box, text="Three", width=10)
        button.pack(side=Tk_.LEFT, padx=1, pady=5)
        button = Tk_.Button(box, text="Four", width=10)
        button.pack(side=Tk_.LEFT, padx=1, pady=5)
        box.pack()

    def buttonbox(self):
        box = Tk_.Frame(self)
        OK = Tk_.Button(box, text="OK", width=10, command=self.ok, default=Tk_.ACTIVE)
        OK.pack(side=Tk_.LEFT, padx=5, pady=5)
        Apply = Tk_.Button(box, text="Apply", width=10, command=self.apply)
        Apply.pack(side=Tk_.LEFT, padx=5, pady=5)
        Cancel = Tk_.Button(box, text="Cancel", width=10, command=self.cancel)
        Cancel.pack(side=Tk_.LEFT, padx=5, pady=5)
        self.bind("<Return>", self.ok)
        self.bind("<Escape>", self.cancel)
        box.pack(side=Tk_.BOTTOM)

    def font_panel(self, body_frame):
        self.list_frame = list_frame = Tk_.Frame(body_frame)
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
        font_list.grid(row=0, column=0)
        self.font_scroller.grid(row=0, column=1, sticky=Tk_.N + Tk_.S)
        list_frame.grid(row=0, column=0)
        self.sample = sample = Tk_.Text(self.body_frame,
                                        width=50, height=4,
                                        highlightthickness=0,
                                        font=self.prefs['font'])
        self.sample.bind('<Button-1>', lambda event: 'break')
        self.sample.insert(Tk_.INSERT, 'ABCDEFGHIJKLMNOPQRSTUVWXYZ\n'\
                                           'abcdefghijklmnopqrstuvwxyz')
        sample.tag_add('all', '1.0', Tk_.END)
        sample.tag_config('all', justify=Tk_.CENTER,
                          font=self.prefs['font'])
        sample.grid(row=1, column=0, pady=10)
        font_list.bind('<ButtonRelease-1>', self.set_font_sample)
        return font_list

    def set_font_sample(self, event):
        index = self.font_list.curselection()[0]
        family = self.families[int(index)]
        size = 18
        weight = 'normal'
        self.sample.tag_config('all', justify=Tk_.CENTER,
                               font=(family, size, weight)) 
        

if __name__ == '__main__':
    parent = Tk_.Tk()
    PreferenceDialog(parent)
    parent.mainloop()
