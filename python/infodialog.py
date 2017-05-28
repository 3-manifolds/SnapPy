# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import sys, os, datetime


if sys.version_info[0] < 3: 
    import Tkinter as Tk_
    from tkSimpleDialog import Dialog
    import ttk
else:
    import tkinter as Tk_
    from tkinter.simpledialog import Dialog
    from tkinter import ttk as ttk
    
from .version import version as SnapPy_version
from IPython import __version__ as IPython_version

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

SnapPy is a program for studying the topology and geometry of 3-manifolds, with a focus on hyperbolic structures. It was written by Marc Culler, Nathan Dunfield, Matthias Gӧrner, and Jeff Weeks, with additional contributions by many others.  Its homepage is

     http://snappy.computop.org/

This is version %s of SnapPy, running on Python %s using Tk %s and IPython %s.

Development of SnapPy was made possible in part by generous support from the National Science Foundation of the United States.

SnapPy is copyright © 2009-%d by Marc Culler, Nathan Dunfield, Matthias Gӧrner, Jeff Weeks, and others and is distributed under the GNU Public License, version 2 or later.  
"""% (SnapPy_version,
      sys.version.split()[0],
      Tk_.Tcl().eval('info patchlevel'),
      IPython_version,
      datetime.datetime.now().year)

def about_snappy(window):
        InfoDialog(window, 'About SnapPy', about_snappy_text)

if __name__ == '__main__':
    root = Tk_.Tk()
    info = InfoDialog(root, title='About SnapPy', content=about_snappy_text)
