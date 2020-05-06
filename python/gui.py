# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import sys
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

class WindowOrFrame:
    def __init__(self, parent = None,
                 root = None, title = '', window_type = 'untyped'):
        if root:
            self.root = root
        else:
            self.root = _get_root(window_type)

        if parent:
            self.container = ttk.Frame(parent)
            self.window = parent
        else:
            self.container = Tk_.Toplevel(
                master = _get_root(window_type),
                class_ = 'snappy')
            self.container.protocol("WM_DELETE_WINDOW", self.close)
            self.container.title(title)
            self.window = self.container

    def close(self, event = None):
        self.container.destroy()

def _get_root(window_type):
    if Tk_._default_root:
        return Tk_._default_root
    root = IPythonTkRoot(window_type = window_type)
    root.withdraw()
    return root
