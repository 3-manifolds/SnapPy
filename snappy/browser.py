# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import sys, os
try:
    import Tkinter as Tk_
    import ttk
except ImportError:
    import tkinter as Tk_
    from tkinter import ttk

oplus = '\x2a01'
ZZ = '\x2124'

# The ttk.LabelFrame is designed to go in a standard window.
# If placed in a ttk.Notebook it will have the wrong background
# color, since the notebook has a darker background than a
# standard window.  This hack fixes that, by overlaying a label with
# the correct background.
  
class NBLabelframe(ttk.Labelframe):
    def __init__(self, master, text=''):
        ttk.Labelframe.__init__(self, master, text=' ')
        self.overlay = Tk_.Label(self, text=text,
                bg=GroupBG,
                padx=12,
                anchor=Tk_.W,
                relief=Tk_.FLAT,
                borderwidth=1,
                highlightbackground=GroupBG,
                highlightcolor=GroupBG)
        self.overlay.place(relwidth=1, x=0, y=-2, bordermode="outside")

if sys.platform == 'darwin':
    WindowBG = 'SystemDialogBackgroundActive'
    GroupBG = '#e0e0e0'
    BrowserBG = '#a8a8a8'
    ST_args = {
        'selectborderwidth' : 0,
        'highlightbackground' : GroupBG,
        'highlightcolor' : GroupBG,
        'readonlybackground' : GroupBG,
        'relief' : Tk_.FLAT,
        'state' : 'readonly'}
else:
    NBLabelframe = ttk.Labelframe
 
class SelectableText(NBLabelframe):
    def __init__(self, master, labeltext=''):
        NBLabelframe.__init__(self, master, text=labeltext)
        self.var = Tk_.StringVar(master)
        self.value = Tk_.Entry(self, textvariable=self.var, **ST_args)
        self.value.pack()

    def set(self, value):
        self.var.set(value)

    def get(self):
        return self.var.get()

class Browser(Tk_.Toplevel):
    def __init__(self, master, manifold):
        self.manifold = manifold
        Tk_.Toplevel.__init__(self, master)
        self.protocol("WM_DELETE_WINDOW", self.close)
        if sys.platform == 'darwin':
            this_dir =  os.path.dirname(__file__)
            Tk_path = os.path.join(this_dir,
                "darwin-tk" + str(Tk_.TkVersion))
            master.tk.call('lappend', 'auto_path', Tk_path)
            master.tk.call('package', 'require', 'mactoolbar')
            self.tk.call('mactoolbar::createbutton',
                'garbage string', 'Hi', 'Hello',
                os.path.join(this_dir, 'info_icon.gif'),
                lambda : None)
            self.tk.call('mactoolbar::create', self._w)
            self.tk.call('set', 'tk::mac::useCompatibilityMetrics', '0')
            self.tk.call('tk::unsupported::MacWindowStyle',
                'style', self._w, 'document',
                ('standardDocument', 'unifiedTitleAndToolbar')
                )
            #print self.tk.call( 'tk::unsupported::MacWindowStyle',
            #    'style', ._w)
        self.config(bg=BrowserBG)
        self.style = ttk.Style(self)
        self.notebook = nb = ttk.Notebook(self)
        self.title(manifold.name())
        self.build_invariants()
        self.dirichlet_frame = Tk_.Frame(self)
        self.horoball_frame = Tk_.Frame(self)
        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(0, weight=1)
        nb.add(self.dirichlet_frame, text='Dirichlet domain')
        nb.add(self.horoball_frame, text='Cusps')
        nb.grid(row=0, column=0, sticky=Tk_.NSEW)
        self.statusbar = Tk_.Frame(self, height=18, bg=BrowserBG)
        self.statusbar.grid(row=1, column=0, sticky=Tk_.NSEW)
        self.update_info()
        # temporary
        self.geometry('600x400')
        
    def build_invariants(self):
        self.invariant_frame = frame = Tk_.Frame(self, bg=GroupBG)
        self.volume = SelectableText(frame, labeltext='Volume')
        self.volume.grid(padx=10, pady=10)
        self.cs = SelectableText(frame, labeltext='Chern-Simons Invariant')
        self.cs.grid(padx=10, pady=10)
        self.homology = SelectableText(frame, labeltext='First Homology')
        self.homology.grid(padx=10, pady=10)
        self.notebook.add(self.invariant_frame,
                          text='Invariants', padding=[0])

    def update_info(self):
        self.volume.set(repr(self.manifold.volume()))
        self.cs.set(repr(self.manifold.chern_simons()))
        group = repr(self.manifold.homology())
        # For some reason these fancy fonts slow things way down!
        #group = group.replace('+', oplus).replace('Z', ZZ)             
        self.homology.set(group)
        
    def close(self):
        self.destroy()
        
if __name__ == '__main__':
    from snappy import *
    M = Manifold('m125')
    root = Tk_.Tk()
    root.withdraw()
    browser = Browser(root, M)
    root.wait_window(browser)
    
    
