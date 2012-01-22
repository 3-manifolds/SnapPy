# -*- coding: utf-8 -*-
from __future__ import unicode_literals
try:  # Python 2
    import Tkinter as Tk_
    from tkSimpleDialog import Dialog
except ImportError: # Python 3
    import tkinter as Tk_
    from tkinter.simpledialog import Dialog

# Try to find a native button
Button = Tk_.Button # fallback

try:
    from ttk import Button  # Python 2 with ttk
except ImportError:
    pass

try:
    from tkinter.ttk import Button  # Python 3 with ttk
except ImportError:
    pass

import sys, os
from snappy import __file__ as snappy_file
snappy_path = os.path.dirname(snappy_file)
icon_file = os.path.join(snappy_path, 'info_icon.gif')

class InfoDialog(Dialog):
    def __init__(self, master, title='', content=''):
        self.content = content
        # We need to figure how to set the background for Windows
        if sys.platform == 'darwin':
            self.bg = 'systemButtonFace'
        elif sys.platform == 'linux2':
            self.bg = Button(master).cget('background')
        else:
            self.bg = 'lightgray'
        self.image = Tk_.PhotoImage(file=icon_file)
        Dialog.__init__(self, master, title=title)
        
    def body(self, master):
        self.config(bg=self.bg)
        self.resizable(False, False)
        box = Tk_.Frame(self, bg=self.bg)
        icon = Tk_.Label(box, image=self.image, bg=self.bg)
        icon.pack(side=Tk_.LEFT, pady=30, anchor=Tk_.N)
        message = Tk_.Message(box, text=self.content, bg=self.bg)
        message.pack(side=Tk_.LEFT, padx=20, pady=10)
        box.pack()

    def buttonbox(self):
        box = Tk_.Frame(self, bg=self.bg)
        button = Button(box, text="OK", width=5,
                       command=self.ok, default=Tk_.ACTIVE)
        button.pack(side=Tk_.RIGHT, ipadx=10, padx=20, pady=20)
        self.bind("<Return>", self.ok)
        box.pack(anchor=Tk_.E)

if __name__ == '__main__':
    from snappy.app import about_snappy
    root = Tk_.Tk()
    info = InfoDialog(root, title='About SnapPy', content=about_snappy)
