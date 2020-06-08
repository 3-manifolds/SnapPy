# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import sys, time
import tkinter as Tk_
from tkinter import ttk as ttk
from tkinter.font import Font, families as font_families
if sys.version_info.major < 3:
    from tkSimpleDialog import Dialog
    from SimpleDialog import SimpleDialog
    from tkMessageBox import askyesno
else:
    from tkinter.simpledialog import Dialog, SimpleDialog
    from tkinter.messagebox import askyesno

from plink.ipython_tools import IPythonTkRoot

if sys.version_info.major < 3 or sys.version_info.minor < 7:
    class Spinbox(ttk.Entry):
        def __init__(self, master=None, **kw):
            ttk.Entry.__init__(self, master, "ttk::spinbox", **kw)

        def set(self, value):
            self.tk.call(self._w, "set", value)
else:
    Spinbox = ttk.Spinbox

class SnapPyStyle:
    def __init__(self):
        self.ttk_style = ttk_style = ttk.Style()
        # The windowBG, groupBG and subgroupBG colors can be used to match Tk objects to
        # Ttk containers.
        if sys.platform == 'darwin':
            try:
                # check if our Tk supports the new semantic colors
                test = Tk_._default_root.winfo_rgb('systemWindowBackgroundColor1')
                self.windowBG = 'systemWindowBackgroundColor'
                self.groupBG = 'systemWindowBackgroundColor1'
                self.subgroupBG = 'systemWindowBackgroundColor2'
            except:
                # These will be wrong in dark mode.
                self.windowBG = '#ededed'
                self.groupBG = '#e5e5e5'
                self.subgroupBG = '#dddddd'
        else:
            self.windowBG = ttk_style.lookup('TLabelframe', 'background')
            self.groupBG = self.subgroupBG = self.windowBG
        self.font = ttk_style.lookup('TLabel', 'font')
        self.font_info = fi = Font(font=self.font).actual()
        fi['size'] = abs(fi['size']) # Why would the size be negative???

class ViewerWindow(Tk_.Toplevel):
    def __init__(self, view_class, *args, **kwargs):
        window_type = kwargs.pop('window_type', 'untyped')
        self.root = kwargs.pop('root', self._get_root(window_type))
        Tk_.Toplevel.__init__(self, master=self.root, class_='snappy')
        self.protocol("WM_DELETE_WINDOW", self.close)
        self.title(string=kwargs.pop('title', ''))
        self.view = view_class(self, *args, **kwargs)
        self.view.pack(expand=True, fill=Tk_.BOTH)
        self.update_idletasks()

    def __repr__(self):
        return 'New window: %s\n'%self.title()

    def _get_root(self, window_type):
        if Tk_._default_root:
            return Tk_._default_root
        root = IPythonTkRoot(window_type = window_type)
        root.withdraw()
        return root

    def close(self):
        if hasattr(self.view, 'dealloc'):
            self.view.dealloc()
        self.view = None
        self.destroy()

    def test(self):
        print('Testing viewer for %s'%self.title())
        time.sleep(0.5)
        if hasattr(self.view, 'test'):
            if not self.view.winfo_ismapped():
                self.view.wait_visibility()
            self.update()
            time.sleep(0.5)
            self.view.focus_force()
            self.view.test()
            time.sleep(0.5)
            self.close()
