# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import sys, os, datetime
from .gui import *
from .version import version as SnapPy_version
from IPython import __version__ as IPython_version

snappy_path = os.path.dirname(__file__)
icon_file = os.path.join(snappy_path, 'info_icon.gif')

class InfoDialog(Dialog):
    def __init__(self, master, title='', content=''):
        self.content = content
        self.style = ttk.Style(master)
        self.image = Tk_.PhotoImage(file=icon_file)
        Dialog.__init__(self, master, title=title)
        
    def body(self, master):
        self.resizable(False, False)
        box = ttk.Frame(self)
        icon = ttk.Label(box, image=self.image)
        icon.grid(row=0, column=0, pady=30, sticky=Tk_.N)
        if sys.platform == 'darwin':
            message = Tk_.Message(box, text=self.content)
        else:
            bgColor=self.style.lookup('Button', 'background')
            message = Tk_.Message(box, text=self.content, bg=bgColor)
        message.grid(row=0, column=1, padx=20, pady=10)
        box.pack()

    def buttonbox(self):
        box = ttk.Frame(self, padding=(0, 0, 0, 20))
        if sys.platform in ('linux', 'linux2'):
            # Work around a bug in linux Tk by using a clunky Tk button.
            w = Tk_.Button(box, text="OK", width=10, command=self.ok, default=Tk_.ACTIVE)
        else:
            w = ttk.Button(box, text="OK", width=10, command=self.ok, default=Tk_.ACTIVE)
        w.grid(padx=5, pady=5)
        self.bind("<Return>", self.ok)
        self.bind("<Escape>", self.ok)
        box.pack()
        
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
