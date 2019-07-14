# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import sys
if sys.version_info[0] < 3: 
    import Tkinter as Tk_
    from tkinter import ttk as ttk
    from tkFont import Font
else:
    import tkinter as Tk_
    from tkinter import ttk as ttk
    from tkinter.font import Font

class SnapPyStyle:
    def __init__(self):
        self.ttk_style = ttk_style = ttk.Style()
        # The windowBG and groupBG colors can be used to match Tk objects to
        # Ttk containers.
        if sys.platform == 'darwin':
            try:
                # Make sure our Tk supports the new semantic colors 
                test = Tk_._default_root.winfo_rgb('systemWindowBackgroundColor')
                self.windowBG = 'systemWindowBackgroundColor'
                self.groupBG = 'systemWindowBackgroundColor1'
                self.subgroupBG = 'systemWindowBackgroundColor2'
            except:
                # These will be wrong in dark mode
                self.windowBG = '#ececec'
                self.groupBG = '#e3e3e3'
                self.subgroupBG = '#dbdbdb'
        else:
            self.windowBG = ttk_style.lookup('TLabelframe', 'background')
            self.groupBG = self.subgroupBG = self.windowBG
        self.font_info = fi = Font(font=ttk_style.lookup('TLabel', 'font')).actual()
        fi['size'] = abs(fi['size']) # Why would the size be negative???
