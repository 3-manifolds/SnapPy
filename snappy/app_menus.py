# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import sys
try:
    import Tkinter as Tk_
    import ttk
except ImportError:
    import tkinter as Tk_
    from tkinter import ttk

OSX_shortcuts = {'Open'   : '\t\t⌘O',
                 'Save'   : '\t\t⌘S',
                 'SaveAs' : '\t⌘⇧S',
                 'Cut'    : '\t⌘X',
                 'Copy'   : '\t⌘C',
                 'Paste'  : '\t⌘V'}

Linux_shortcuts = {'Open'   : '',
                   'Save'   : '',
                   'SaveAs' : '',
                   'Cut'    : '     Cntl+X',
                   'Copy'   : '',
                   'Paste'  : '  Cntl+V'}

if sys.platform == 'darwin' :
    scut = OSX_shortcuts
elif sys.platform == 'linux2' :
    scut = Linux_shortcuts
else: # fall back choice
    scut = Linux_shortcuts

def togl_save_image(self):
    savefile = filedialog.asksaveasfile(
        mode='wb',
        title='Save Image As PNG Image File',
        defaultextension = '.png',
        filetypes = [
            ("PNG image files", "*.png *.PNG", ""),
            ("All files", "")])
    self.to_front()
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
    File_menu.add_command(
        label='Open...' + scut['Open'], state='disabled')
    File_menu.add_command(
        label='Save as...' + scut['SaveAs'], command=self.save)
    File_menu.add_separator()
    File_menu.add_command(label='Close', command=self.close)
    menubar.add_cascade(label='File', menu=File_menu)
    Edit_menu = Tk_.Menu(menubar, name='edit')
    Edit_menu.add_command(
        label='Cut' + scut['Cut'], state='disabled')
    Edit_menu.add_command(
        label='Copy' + scut['Copy'], state='disabled')
    Edit_menu.add_command(
        label='Paste' + scut['Paste'], state='disabled')
    Edit_menu.add_command(
        label='Delete', state='disabled')
    menubar.add_cascade(label='Edit', menu=Edit_menu)
    if self.window_master is not None:
        Window_menu = self.window_master.menubar.children['window']
        menubar.add_cascade(label='Window', menu=Window_menu)
        Help_menu = Tk_.Menu(menubar, name="help")
        Help_menu.add_command(label='Help on SnapPy ...',
                              command=self.window_master.howto)
        menubar.add_cascade(label='Help', menu=Help_menu)

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
    File_menu.add_command(
        label='Open...' + scut['Open'], state='disabled')
    File_menu.add_command(
        label='Save as...' + scut['SaveAs'], state='disabled')
    File_menu.add_command(label='Save Image...', command=self.save_image)
    File_menu.add_separator()
    File_menu.add_command(label='Close', command=self.close)
    menubar.add_cascade(label='File', menu=File_menu)
    Edit_menu = Tk_.Menu(menubar, name='edit')
    Edit_menu.add_command(
        label='Cut' + scut['Cut'], state='disabled')
    Edit_menu.add_command(
        label='Copy' + scut['Copy'], state='disabled')
    Edit_menu.add_command(
        label='Paste' + scut['Paste'], state='disabled')
    Edit_menu.add_command(
        label='Delete', state='disabled')
    menubar.add_cascade(label='Edit', menu=Edit_menu)
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
        label='Open...' + scut['Open'], state='disabled')
    File_menu.add_command(
        label='Save as...' + scut['SaveAs'], state='disabled')
    Print_menu = Tk_.Menu(menubar, name='print')
    File_menu.add_command(label='Save Image...', command=self.save_image)
    File_menu.add_separator()
    File_menu.add_command(label='Close', command=self.close)
    menubar.add_cascade(label='File', menu=File_menu)
    Edit_menu = Tk_.Menu(menubar, name='edit')
    Edit_menu.add_command(
        label='Cut' + scut['Cut'], state='disabled')
    Edit_menu.add_command(
        label='Copy' + scut['Copy'], state='disabled')
    Edit_menu.add_command(
        label='Paste' + scut['Paste'], state='disabled')
    Edit_menu.add_command(
        label='Delete', state='disabled')
    menubar.add_cascade(label='Edit', menu=Edit_menu)
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
