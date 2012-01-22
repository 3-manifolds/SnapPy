# -*- coding: utf-8 -*-
from __future__ import unicode_literals
try:
    import Tkinter as Tk_
    from tkSimpleDialog import Dialog
    import ttk
except ImportError: # Python 3
    import tkinter as Tk_
    from tkinter.simpledialog import Dialog
    from tkinter import ttk
import sys, os
from snappy import __file__ as snappy_file
snappy_path = os.path.dirname(snappy_file)
# Once we switch to Tk 8.5 or higher, get rid of this.
if Tk_.TkVersion == 8.4:
    icon_file = os.path.join(snappy_path, 'info_icon.gif')
else:
    icon_file = os.path.join(snappy_path, 'info_icon.png')

# We need to figure how to set the background for Windows and linux
if sys.platform =='darwin':
    background = 'systemButtonFace'
else:
    background = 'lightgray'

class InfoDialog(Dialog):
    def __init__(self, master, title='', content=''):
        self.content = content
        self.image = Tk_.PhotoImage(file=icon_file)
        Dialog.__init__(self, master, title=title)
        
    def body(self, master):
        self.config(bg=background)
        self.resizable(False, False)
        box = Tk_.Frame(self, bg=background)
        icon = Tk_.Label(box, image=self.image, bg=background)
        icon.pack(side=Tk_.LEFT, pady=30, anchor=Tk_.N)
        message = Tk_.Message(box, text=self.content, bg=background)
        message.pack(side=Tk_.LEFT, padx=20, pady=10)
        box.pack()

    def buttonbox(self):
        box = Tk_.Frame(self, bg=background)
        button = ttk.Button(box, text="OK", width=5,
                       command=self.ok, default=Tk_.ACTIVE)
        button.pack(side=Tk_.RIGHT, ipadx=10, padx=20, pady=20)
        self.bind("<Return>", self.ok)
        box.pack(anchor=Tk_.E)

if __name__ == '__main__':
    from snappy.app import about_snappy
    root = Tk_.Tk()
    info = InfoDialog(root, title='About SnapPy', content=about_snappy)
