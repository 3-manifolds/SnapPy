# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import sys
if sys.version_info[0] < 3: 
    import Tkinter as Tk_
    import ttk
    from tkFont import Font
else:
    import tkinter as Tk_
    from tkinter import ttk as ttk
    from tkinter.font import Font

class _SnapPyStyle:
    def __init__(self):
        self.ttk_style = ttk_style = ttk.Style()
        self.WindowBG = self.GroupBG = ttk_style.lookup('TLabelframe', 'background')
        self.font_info = fi = Font(font=ttk_style.lookup('TLabel', 'font')).actual()
        fi['size'] = abs(fi['size']) # Why would the size be negative???

    def configure(self):
        ttk_style = self.ttk_style

def SnapPyStyle(root):
    if root is None:
        if Tk_._default_root is None:
            root = Tk_.Tk(className='snappy')
            root.withdraw()
        else:
            root = Tk_._default_root
    try:
        return root.style
    except AttributeError:
        root.style = style = _SnapPyStyle()
        style.configure()
        return style 
