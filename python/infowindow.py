# -*- coding: utf-8 -*-
import sys
import os
import datetime
from .gui import *
from .version import version as SnapPy_version
from IPython import __version__ as IPython_version

snappy_path = os.path.dirname(__file__)
icon_file = os.path.join(snappy_path, 'info_icon.gif')


class InfoWindow(Tk_.Toplevel):
    def __init__(self, root, title, content):
        self.root = root
        Tk_.Toplevel.__init__(self, self.root, class_='snappy')
        self.title(title)
        self.content = content
        self.image = Tk_.PhotoImage(file=icon_file)
        self.style = SnapPyStyle()
        self.body()

    def body(self):
        self.resizable(False, False)
        box = ttk.Frame(self)
        icon = ttk.Label(box, image=self.image)
        icon.grid(row=0, column=0, pady=30, sticky=Tk_.N)
        message = Tk_.Message(box, text=self.content)
        message.grid(row=0, column=1, padx=20, pady=10)
        box.pack()


about_snappy_text = """
For information on how to use SnapPy, please see the Help menu.

SnapPy is a program for studying the topology and geometry of 3-manifolds, with a focus on hyperbolic structures. It was written by Marc Culler, Nathan Dunfield, Matthias Gӧrner, and Jeff Weeks, with additional contributions by many others.  Its homepage is

     http://snappy.computop.org/

This is version %s of SnapPy, running on Python %s using Tk %s and IPython %s.

Development of SnapPy was made possible in part by generous support from the National Science Foundation of the United States.

SnapPy is copyright © 2009-%d by Marc Culler, Nathan Dunfield, Matthias Gӧrner, Jeff Weeks, and others and is distributed under the GNU Public License, version 2 or later.
""" % (SnapPy_version,
       sys.version.split()[0],
       Tk_.Tcl().eval('info patchlevel'),
       IPython_version,
       datetime.datetime.now().year)


def about_snappy(window):
    InfoWindow(window, title='About SnapPy', content=about_snappy_text)


if __name__ == '__main__':
    root = Tk_.Tk()
    info = InfoWindow(root, title='About SnapPy', content=about_snappy_text)
