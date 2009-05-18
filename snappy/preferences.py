import Tkinter as Tk_
import tkSimpleDialog
import tkFont
from string import ascii_letters

class PreferenceDialog(tkSimpleDialog.Dialog):
    def __init__(self, parent, title='Preferences'):
        Tk_.Toplevel.__init__(self, parent)
        self.title(title)
        self.parent = parent
        self.result = None
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

    def font_panel(self, bodyframe):
        self.font_list = font_list = Tk_.Listbox(bodyframe)
        fonts = [x for x in tkFont.families() if x[0] in ascii_letters]
        fonts.sort()
        for font in fonts:
            font_list.insert(Tk_.END, font)
        self.font_scroller = font_scroller = Tk_.Scrollbar(bodyframe,
                                                           command=font_list.yview)
        font_list.config(yscrollcommand=font_scroller.set)
        font_list.grid(row=0, column=0)
        self.font_scroller.grid(row=0, column=1, sticky=Tk_.N + Tk_.S)
        return font_list

if __name__ == '__main__':
    parent = Tk_.Tk()
    PreferenceDialog(parent)
    parent.mainloop()
