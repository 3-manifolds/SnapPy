# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import sys, os, datetime
try:  # Python 2
    import Tkinter as Tk_
    from tkSimpleDialog import Dialog
    import ttk
except ImportError: # Python 3
    import tkinter as Tk_
    from tkinter.simpledialog import Dialog
    from tkinter import ttk
from snappy.version import version as SnapPy_version

snappy_path = os.path.dirname(__file__)
icon_file = os.path.join(snappy_path, 'info_icon.gif')

class InfoDialog(Dialog):
    def __init__(self, master, title='', content=''):
        self.content = content
        self.style = ttk.Style(master)
        if sys.platform == 'darwin':
            self.bg = 'SystemDialogBackgroundActive'
        else:
            self.bg = self.style.lookup('Button', 'background')
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
        button = ttk.Button(box, text="OK", width=5,
                       command=self.ok, default=Tk_.ACTIVE)
        button.pack(side=Tk_.RIGHT, ipadx=10, padx=20, pady=20)
        self.bind("<Return>", self.ok)
        box.pack(anchor=Tk_.E)

about_snappy_text = """
For information on how to use SnapPy, please see the Help menu.

SnapPy is a user interface for the SnapPea kernel, which was originally written by Jeff Weeks.  The kernel has been extended by contributors to the SnapPy project, including Marc Culler, Nathan Dunfield and Matthias Gӧrner.

SnapPy was written by Marc Culler and Nathan Dunfield and is distributed under the GNU Public License, version 2 or later.  Its home page is:
     http://snappy.computop.org/

The release number of this SnapPy is %s.

SnapPy is written in the Python language, using Cython to incorporate the SnapPea kernel code. The graphical interface uses Tcl/Tk, via Python's Tkinter module.

Information, downloads, and source code for the SnapPea kernel and for the user interfaces written by Jeff Weeks are available at:
     http://www.geometrygames.org/SnapPea-old/
     http://www.geometrygames.org/SnapPea/

Development of SnapPy was made possible in part by generous support from the National Science Foundation of the United States.

Copyright © 2009-%d, Marc Culler, Nathan Dunfield, and others.
"""% (SnapPy_version, datetime.datetime.now().year)

def about_snappy(window):
        InfoDialog(window, 'About SnapPy', about_snappy_text)

if __name__ == '__main__':
    root = Tk_.Tk()
    info = InfoDialog(root, title='About SnapPy', content=about_snappy_text)
