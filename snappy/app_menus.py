# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import sys
try:
    import Tkinter as Tk_
    import ttk
except ImportError:
    import tkinter as Tk_
    from tkinter import ttk

from snappy import filedialog
import tempfile, png, os
    
OSX_shortcuts = {'Open...'    : 'Command-o',
                 'Save'       : 'Command-s',
                 'Save as...' : 'Command-Shift-s',
                 'Cut'        : 'Command-x',
                 'Copy'       : 'Command-c',
                 'Paste'      : 'Command-v',
                 'Close'     : 'Command-w', 
                 'Left'       : '←',
                 'Up'         : '↑',
                 'Right'      : '→',
                 'Down'       : '↓'}

OSX_shortcut_events = {'Open...' : '<Command-o>',
                 'Save'          : '<Command-s>',
                 'Save as...'    : '<Command-S>',
                 'Cut'           : '<Command-x>',
                 'Copy'          : '<Command-c>',
                 'Paste'         : '<Command-v>',
                 'Close'     : '<Command-w>'
                 }
                 
Linux_shortcuts = {'Open...'    : 'Cntl+O',
                   'Save'       : 'Cntl+S',
                   'Save as...' : 'Cntl+Shift+S',
                   'Cut'        : 'Cntl+X',
                   'Copy'       : 'Cntl+W',
                   'Paste'      : 'Cntl+V',
                   'Left'       : '←',
                   'Up'         : '↑',
                   'Right'      : '→',
                   'Down'       : '↓'}

Linux_shortcut_events = {'Open...'   : '<Control-o>',
                         'Save'      : '<Control-s>',
                         'Save as...' : '<Control-S>',
                         'Cut'       : '<Control-x>',
                         'Copy'      : '<Control-w>',
                         'Paste'     : '<Control-v>',
                          }

if sys.platform == 'darwin' :
    scut = OSX_shortcuts
    scut_events = OSX_shortcut_events
elif sys.platform == 'linux2' :
    scut = Linux_shortcuts
    scut_events = Linux_shortcut_events
else: # fall back choice
    scut = Linux_shortcuts
    scut_events = Linux_shortcut_events

def add_menu(root, menu, label, command, state='active'):
    accelerator = scut.get(label, '')
    menu.add_command(label=label, accelerator=accelerator,
            command=command, state=state)
    if scut_events.get(label, None) and state != 'disabled':
        root.bind(scut_events[label], command)

def really_disable_menu_items(menubar):
    """
    On OS X, menu items aren't greying out as they should.
    """
    for label, menu in menubar.children.items():
        try:
            for i in range(menu.index('end') + 1):
                if menu.type(i) == 'command':
                    if menu.entrycget(i, 'state') == 'disabled':
                        menu.entryconfig(i, state='disabled')
        except:
            pass

def add_edit_menu_with_disabled_items(menubar, window):
    edit_menu = Tk_.Menu(menubar, name='edit')
    add_menu(window, edit_menu, 'Cut', None, state='disabled')
    add_menu(window, edit_menu, 'Copy', None, state='disabled')
    add_menu(window, edit_menu, 'Paste', None, state='disabled')
    add_menu(window, edit_menu, 'Delete', None, state='disabled')
    menubar.add_cascade(label='Edit', menu=edit_menu)

def add_window_and_snappy_help(self):
    if self.window_master is not None:
        Window_menu = self.window_master.menubar.children['window']
        self.menubar.add_cascade(label='Window', menu=Window_menu)
        Help_menu = Tk_.Menu(self.menubar, name="help")
        Help_menu.add_command(label='Help on SnapPy ...',
                              command=self.window_master.howto)
        self.menubar.add_cascade(label='Help', menu=Help_menu)

def togl_save_image(self):
    savefile = filedialog.asksaveasfile(
        mode='wb',
        title='Save Image As PNG Image File',
        defaultextension = '.png',
        filetypes = [
            ("PNG image files", "*.png *.PNG", ""),
            ("All files", "")])
    self.window.update()
    self.widget.redraw()
    self.window.update()
    self.widget.redraw()
    if savefile:
        ppm_file = tempfile.mktemp() + ".ppm"
        PI = Tk_.PhotoImage()
        self.widget.tk.call(self.widget._w, 'takephoto', PI.name)
        PI.write(ppm_file, format='ppm')
        infile = open(ppm_file, 'rb')
        format, width, height, depth, maxval = \
                png.read_pnm_header(infile, ('P5','P6','P7'))
        greyscale = depth <= 2
        pamalpha = depth in (2,4)
        supported = [2**x-1 for x in range(1,17)]
        mi = supported.index(maxval)
        bitdepth = mi+1
        writer = png.Writer(width, height,
                        greyscale=greyscale,
                        bitdepth=bitdepth,
                        alpha=pamalpha)
        writer.convert_pnm(infile, savefile)
        savefile.close()
        infile.close()
        os.remove(ppm_file)

def browser_menus(self):
    self.menubar = menubar = Tk_.Menu(self.window)
    Python_menu = Tk_.Menu(menubar, name="apple")
    Python_menu.add_command(label='About SnapPy ...')
    Python_menu.add_separator()
    Python_menu.add_command(label='Preferences ...', state='disabled')
    Python_menu.add_separator()
    if sys.platform == 'linux2' and self.window_master is not None:
        Python_menu.add_command(label='Quit SnapPy', command=
                                self.window_master.close)
    menubar.add_cascade(label='SnapPy', menu=Python_menu)
    File_menu = Tk_.Menu(menubar, name='file')
    add_menu(self.window, File_menu, 'Open...', None, 'disabled')
    add_menu(self.window, File_menu, 'Save as...', self.save)
    File_menu.add_separator()
    add_menu(self.window, File_menu, 'Close', self.close)
    menubar.add_cascade(label='File', menu=File_menu)
    add_edit_menu_with_disabled_items(menubar, self.window)
    add_window_and_snappy_help(self)
        

def dirichlet_menus(self):
    self.menubar = menubar = Tk_.Menu(self.window)
    Python_menu = Tk_.Menu(menubar, name="apple")
    Python_menu.add_command(label='About SnapPy ...')
    Python_menu.add_separator()
    Python_menu.add_command(label='Preferences ...', state='disabled')
    Python_menu.add_separator()
    if sys.platform == 'linux2' and self.window_master is not None:
        Python_menu.add_command(label='Quit SnapPy', command=
                                self.window_master.close)
    menubar.add_cascade(label='SnapPy', menu=Python_menu)
    File_menu = Tk_.Menu(menubar, name='file')
    add_menu(self.window, File_menu, 'Open...', None, 'disabled')
    add_menu(self.window, File_menu, 'Save as...', None, 'disabled')
    File_menu.add_command(label='Save Image...', command=self.save_image)
    File_menu.add_separator()
    add_menu(self.window, File_menu, 'Close', command=self.close)
    menubar.add_cascade(label='File', menu=File_menu)
    add_edit_menu_with_disabled_items(menubar, self.window)
    if self.window_master:
        Window_menu = self.window_master.menubar.children['window']
        menubar.add_cascade(label='Window', menu=Window_menu)
    Help_menu = Tk_.Menu(menubar, name="help")
    Help_menu.add_command(label='Help on PolyhedronViewer ...',
                          command=self.widget.help)
    menubar.add_cascade(label='Help', menu=Help_menu)

def horoball_menus(self):
    self.menubar = menubar = Tk_.Menu(self.window)
    Python_menu = Tk_.Menu(menubar, name="apple")
    Python_menu.add_command(label='About SnapPy ...')
    Python_menu.add_separator()
    Python_menu.add_command(label='Preferences ...',  state='disabled')
    Python_menu.add_separator()
    if sys.platform == 'linux2' and self.window_master is not None:
        Python_menu.add_command(label='Quit SnapPy',
                                command=self.window_master.close)
    menubar.add_cascade(label='SnapPy', menu=Python_menu)
    File_menu = Tk_.Menu(menubar, name='file')
    File_menu.add_command(
        label='Open...', accelerator=scut['Open...'], state='disabled')
    File_menu.add_command(
        label='Save as...', accelerator=scut['Save as...'], state='disabled')
    Print_menu = Tk_.Menu(menubar, name='print')
    File_menu.add_command(label='Save Image...', command=self.save_image)
    File_menu.add_separator()
    File_menu.add_command(label='Close', command=self.close)
    menubar.add_cascade(label='File', menu=File_menu)
    add_edit_menu_with_disabled_items(menubar, self.window)
    View_menu = Tk_.Menu(menubar, name='view')
    View_menu.add_checkbutton(label='parallelogram',
                              command=self.view_check,
                              variable=self.pgram_var)
    View_menu.add_checkbutton(label='Ford edges',
                              command=self.view_check,
                              variable=self.Ford_var)
    View_menu.add_checkbutton(label='triangulation',
                              command=self.view_check,
                              variable=self.tri_var)
    View_menu.add_checkbutton(label='horoballs',
                              command=self.view_check,
                              variable=self.horo_var)
    View_menu.add_checkbutton(label='labels',
                              command=self.view_check,
                              variable=self.label_var)
    menubar.add_cascade(label='View', menu=View_menu)
    if self.window_master is not None:
        Window_menu = self.window_master.menubar.children['window']
        menubar.add_cascade(label='Window', menu=Window_menu)
    Help_menu = Tk_.Menu(menubar, name="help")
    Help_menu.add_command(label='Help on HoroballViewer ...',
                          command=self.widget.help)
    menubar.add_cascade(label='Help', menu=Help_menu)


def link_menus(self):
    self.menubar = menubar = Tk_.Menu(self.window)
    Python_menu = Tk_.Menu(menubar, name="apple")
    Python_menu.add_command(label='About SnapPy ...')
    Python_menu.add_separator()
    Python_menu.add_command(label='Preferences ...', state='disabled')
    Python_menu.add_separator()
    if sys.platform == 'linux2' and self.window_master is not None:
        Python_menu.add_command(label='Quit SnapPy', command=
                                self.window_master.close)
    menubar.add_cascade(label='SnapPy', menu=Python_menu)
    File_menu = Tk_.Menu(menubar, name='file')
    add_menu(self.window, File_menu, 'Open...', None, 'disabled')
    add_menu(self.window, File_menu, 'Save as...', None, 'disabled')
    self.build_save_image_menu(menubar, File_menu)
    File_menu.add_separator()
    add_menu(self.window, File_menu, 'Close', command=self.close)
    menubar.add_cascade(label='File', menu=File_menu)
    add_edit_menu_with_disabled_items(menubar, self.window)
    add_window_and_snappy_help(self)
