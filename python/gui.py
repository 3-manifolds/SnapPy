# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import sys
import tkinter as Tk_
from tkinter import ttk as ttk
from tkinter.font import Font, families as font_families
if sys.version_info[0] < 3:
    from tkSimpleDialog import Dialog
    from SimpleDialog import SimpleDialog
else:
    from tkinter.simpledialog import Dialog, SimpleDialog

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
            self.windowBG = ttk_style.lookup('TFrame', 'background')
            self.groupBG = ttk_style.lookup('TLabelframe', 'background')
            try:
                # check if our Tk supports the new semantic colors
                test = Tk_._default_root.winfo_rgb('systemWindowBackgroundColor')
                self.subgroupBG = 'systemWindowBackgroundColor2'
            except:
                # This will be wrong in dark mode.
                self.subgroupBG = '#dbdbdb'
        else:
            self.windowBG = ttk_style.lookup('TLabelframe', 'background')
            self.groupBG = self.subgroupBG = self.windowBG
        self.font_info = fi = Font(font=ttk_style.lookup('TLabel', 'font')).actual()
        fi['size'] = abs(fi['size']) # Why would the size be negative???
