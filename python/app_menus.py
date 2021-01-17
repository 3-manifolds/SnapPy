# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from builtins import range
import sys, os, webbrowser
from urllib.request import pathname2url
from .gui import *
from . import __file__ as snappy_dir
from .infodialog import about_snappy, InfoDialog

OSX_shortcuts = {'Open...'    : 'Command-o',
                 'Save'       : 'Command-s',
                 'Save as...' : 'Command-Shift-s',
                 'Cut'        : 'Command-x',
                 'Copy'       : 'Command-c',
                 'Paste'      : 'Command-v',
                 'Close'      : 'Command-w',
                 'Left'       : '←',
                 'Up'         : '↑',
                 'Right'      : '→',
                 'Down'       : '↓'}

OSX_shortcut_events = {'Open...'     : '<Command-o>',
                       'Save'        : '<Command-s>',
                       'Save as...'  : '<Command-S>',
                       'Cut'         : '<Command-x>',
                       'Copy'        : '<Command-c>',
                       'Paste'       : '<Command-v>',
                       'Close'       : '<Command-w>'
                   }

Linux_shortcuts = {'Open...'    : 'Cntl+O',
                   'Save'       : 'Cntl+S',
                   'Save as...' : 'Cntl+Shift+S',
                   'Cut'        : 'Cntl+X',
                   'Copy'       : 'Cntl+Shift+C',
                   'Paste'      : 'Cntl+Shift+V',
                   'Left'       : '←',
                   'Up'         : '↑',
                   'Right'      : '→',
                   'Down'       : '↓'}

Linux_shortcut_events = {'Open...'    : '<Control-o>',
                         'Save'       : '<Control-s>',
                         'Save as...' : '<Control-S>',
                         'Cut'        : '<Control-x>',
                         'Copy'       : '<Control-C>',
                         'Paste'      : '<Control-V>',
                          }

help_report_bugs_label = 'Report bugs ...'
help_polyhedron_viewer_label = 'Polyhedron Viewer Help ...'
help_horoball_viewer_label = 'Horoball Viewer Help ...'

if sys.platform == 'darwin' :
    scut = OSX_shortcuts
    scut_events = OSX_shortcut_events
elif sys.platform == 'linux2' or sys.platform == 'linux' :
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

class EditMenu(Tk_.Menu):
    """Edit Menu cascade containing Cut, Copy, Paste and Delete. To use,
    provide a callback function which returns a dict specifying which
    editing functions should be enabled.  The keys should be chosen
    from the list ['Cut', 'Copy', 'Paste, 'Delete'] and the values
    should be functions to be call for the corresponding actions.  If
    a key is missing, the corresponding menu entry will be disabled.

    """
    entries = ['Cut', 'Copy', 'Paste', 'Delete']

    def __init__(self, menubar, callback):
        Tk_.Menu.__init__(self, menubar, name='snappyedit', postcommand=self.configure)
        self.get_actions = callback
        self.add_entry('Cut', lambda event=None: self.actions['Cut']())
        self.add_entry('Copy', lambda event=None: self.actions['Copy']())
        self.add_entry('Paste', lambda event=None: self.actions['Paste']())
        self.add_entry('Delete', lambda event=None: self.actions['Delete']())
        self.actions = {}

    def add_entry(self, label, command):
        accelerator = scut.get(label, '')
        self.add_command(label=label, accelerator=accelerator,
                         command=command, state='disabled')

    def configure(self):
        """Called before the menu is opened."""
        self.actions = self.get_actions()
        for entry in self.entries:
            if self.actions.get(entry, None):
                self.entryconfig(entry, state='normal')
            else:
                self.entryconfig(entry, state='disabled')

class HelpMenu(Tk_.Menu):
    """Help Menu cascade.  Always contains the main SnapPy help entry.
    Additional help entries for specific tools, such as a Dirichlet
    viewer, may be added or removed.

    """
    def __init__(self, menubar):
        # on OS X setting name='help' makes this a system help menu.
        Tk_.Menu.__init__(self, menubar, name='help')
        if sys.platform != 'darwin':
            self.add_command(label = 'SnapPy Help ...', command=self.show_SnapPy_help)
        self.add_command(label = help_report_bugs_label,
                         command=self.show_bugs_page)
        self.extra_commands = {}

    def show_SnapPy_help(self):
        self.show_page('index.html')

    def show_bugs_page(self):
        self.show_page('bugs.html')

    def show_page(self, page):
        path = os.path.join(os.path.dirname(snappy_dir), 'doc', page)
        if os.path.exists(path):
            url = 'file:' + pathname2url(path)
            try:
                webbrowser.open_new_tab(url)
            except webbrowser.Error:
                from tkinter import messagebox
                messagebox.showwarning('Error', 'Failed to open the documentation file.')
        else:
            from tkinter import messagebox
            messagebox.showwarning('Not found!',
                                   'The file %s does not exist.'%path)

    def extra_command(self, label, command):
        self.extra_commands[label] = command

    def activate(self, labels):
        """Manage extra help entries.
        Pass the labels of the extra commands to be activated.
        """
        end = self.index(Tk_.END)
        if sys.platform == 'darwin':
            self.delete(0, self.index(Tk_.END))
        elif end > 0:
            self.delete(1, self.index(Tk_.END))
        for label in labels:
            if label in self.extra_commands:
                self.add_command(label=label, command=self.extra_commands[label])

class ListedWindow():
    """
    Mixin class that allows emulation of the Apple Window menu. Windows
    register when they open by calling the class method register_window.  They
    call unregister_window when they close.  The class maintains a list of all
    openwindows.Participating windows should be subclasses of ListedWindow, as
    should objects which need access to its list of all windows in the app, such
    as the Preferences object.
    """
    window_list = []
    prefs = {}

    @classmethod
    def register_window(cls, window):
        assert isinstance(window, ListedWindow)
        cls.window_list.append(window)
        window.apply_prefs()

    @classmethod
    def unregister_window(cls, window):
        try:
            cls.window_list.remove(window)
        except ValueError:
            pass

    def bring_to_front(self):
        self.deiconify()
        self.lift()
        self.focus_force()

    def apply_prefs(self):
        # Subclasses should override this if they use preferences.
        pass

class WindowMenu(Tk_.Menu, ListedWindow):
    """
    Menu with a postcommand which shows all listed windows.
    """
    def __init__(self, menubar):
        Tk_.Menu.__init__(self, menubar, name='window', postcommand=self.build_entries)

    def build_entries(self):
        if sys.platform == 'darwin':
            return
        self.delete(0, self.index(Tk_.END))
        for object in self.window_list:
            self.add_command(label=object.menu_title, command=object.bring_to_front)

def browser_menus(self):
    """
    Menus for the browser window.  Used as Browser.build_menus.
    Creates a menubar attribute for the browser.
    """
    self.menubar = menubar = Tk_.Menu(self)
    Python_menu = Tk_.Menu(menubar, name="apple")
    Python_menu.add_command(label='About SnapPy ...',
                            command=lambda : about_snappy(self))
    Python_menu.add_separator()
    if sys.platform in ('linux2', 'linux') and self.main_window is not None:
        Python_menu.add_command(label='Quit SnapPy', command=self.main_window.close)
    menubar.add_cascade(label='SnapPy', menu=Python_menu)
    File_menu = Tk_.Menu(menubar, name='file')
    add_menu(self, File_menu, 'Open...', None, 'disabled')
    add_menu(self, File_menu, 'Save as...', self.save)
    File_menu.add_separator()
    add_menu(self, File_menu, 'Close', self.close)
    menubar.add_cascade(label='File', menu=File_menu)
    menubar.add_cascade(label='Edit ', menu=EditMenu(menubar, self.edit_actions))
    if sys.platform == 'darwin':
        menubar.add_cascade(label='View', menu=Tk_.Menu(menubar, name='view'))
    menubar.add_cascade(label='Window', menu=WindowMenu(menubar))
    help_menu = HelpMenu(menubar)

    def dirichlet_help():
        InfoDialog(self, 'Viewer Help', self.dirichlet_viewer.widget.help_text)
    help_menu.extra_command(label=help_polyhedron_viewer_label, command=dirichlet_help)

    def horoball_help():
        InfoDialog(self, 'Viewer Help', self.horoball_viewer.widget.help_text)
    help_menu.extra_command(label=help_horoball_viewer_label, command=horoball_help)
    menubar.add_cascade(label='Help', menu=help_menu)

def plink_menus(self):
    """Menus for the SnapPyLinkEditor."""
    self.menubar = menubar = Tk_.Menu(self.window)
    Python_menu = Tk_.Menu(menubar, name="apple")
    Python_menu.add_command(label='About PLink...', command=self.about)
    Python_menu.add_separator()
    if sys.platform in ('linux2', 'linux') and self.main_window is not None:
        Python_menu.add_command(label='Quit SnapPy', command=self.main_window.close)
    menubar.add_cascade(label='SnapPy', menu=Python_menu)
    File_menu = Tk_.Menu(menubar, name='file')
    add_menu(self.window, File_menu, 'Open...', self.load)
    add_menu(self.window, File_menu, 'Save as...', self.save)
    self.build_save_image_menu(menubar, File_menu) # Add image save menu
    File_menu.add_separator()
    if self.callback:
        add_menu(self.window, File_menu, 'Close', self.done)
    else:
        add_menu(self.window, File_menu, 'Exit', self.done)
    menubar.add_cascade(label='File', menu=File_menu)
    Edit_menu = Tk_.Menu(menubar, name='snappyedit')
    add_menu(self.window, Edit_menu, 'Cut', None, state='disabled')
    add_menu(self.window, Edit_menu, 'Copy', None, state='disabled')
    add_menu(self.window, Edit_menu, 'Paste', None, state='disabled')
    add_menu(self.window, Edit_menu, 'Delete', None, state='disabled')
    menubar.add_cascade(label='Edit ', menu=Edit_menu)
    self._add_info_menu()
    self._add_tools_menu()
    self._add_style_menu()
    menubar.add_cascade(label='Window', menu=WindowMenu(menubar))
    Help_menu = Tk_.Menu(menubar, name="help")
    menubar.add_cascade(label='Help', menu=HelpMenu(menubar))
    Help_menu.add_command(label='PLink Help ...', command=self.howto)
    self.window.config(menu=menubar)

def dirichlet_menus(self):
    """
    Menus for the standalone Dirichlet viewer.  Called by the view Frame, not the
    master Toplevel.
    """
    self.menubar = menubar = Tk_.Menu(self.master)
    Python_menu = Tk_.Menu(menubar, name="apple")
    Python_menu.add_command(label='About SnapPy ...',
                            command=lambda : about_snappy(self.master))
    Python_menu.add_separator()
    if sys.platform in ('linux2', 'linux') and self.main_window is not None:
        Python_menu.add_command(label='Quit SnapPy', command=
                                self.main_window.close)
    menubar.add_cascade(label='SnapPy', menu=Python_menu)
    File_menu = Tk_.Menu(menubar, name='file')
    add_menu(self.master, File_menu, 'Open...', None, 'disabled')
    add_menu(self.master, File_menu, 'Save as...', None, 'disabled')
    File_menu.add_command(label='Save Image...', command=self.master.save_image)
    Export_menu = Tk_.Menu(File_menu, name='export')
    File_menu.add_cascade(label='Export as STL...', menu=Export_menu)
    Export_menu.add_command(label='Export STL', command=self.export_stl)
    Export_menu.add_command(label='Export Cutout STL', command=self.export_cutout_stl)
    File_menu.add_separator()
    add_menu(self.master, File_menu, 'Close', command=self.master.close)
    menubar.add_cascade(label='File', menu=File_menu)
    menubar.add_cascade(label='Edit ', menu=EditMenu(menubar, self.master.edit_actions))
    menubar.add_cascade(label='Window', menu=WindowMenu(menubar))
    help_menu = HelpMenu(menubar)
    help_menu.extra_command(label=help_polyhedron_viewer_label, command=self.widget.help)
    help_menu.activate([help_polyhedron_viewer_label, help_report_bugs_label])
    self.menubar.add_cascade(label='Help', menu=help_menu)

def horoball_menus(self):
    """
    Menus for the standalone Horoball viewer.  Called by the view Frame, not the
    master Toplevel.
    """
    self.menubar = menubar = Tk_.Menu(self.master)
    Python_menu = Tk_.Menu(menubar, name="apple")
    Python_menu.add_command(label='About SnapPy ...',
                            command=lambda : about_snappy(self.master))
    Python_menu.add_separator()
    if sys.platform in ('linux2', 'linux') and self.main_window is not None:
        Python_menu.add_command(label='Quit SnapPy',
                                command=self.main_window.close)
    menubar.add_cascade(label='SnapPy', menu=Python_menu)
    File_menu = Tk_.Menu(menubar, name='file')
    add_menu(self.master, File_menu, 'Open...', None, 'disabled')
    add_menu(self.master, File_menu, 'Save as...', None, 'disabled')
    File_menu.add_command(label='Save Image...', command=self.master.save_image)
    File_menu.add_separator()
    add_menu(self.master, File_menu, 'Close', command=self.master.close)
    menubar.add_cascade(label='File', menu=File_menu)
    menubar.add_cascade(label='Edit ', menu=EditMenu(menubar, self.master.edit_actions))
    menubar.add_cascade(label='Window', menu=WindowMenu(menubar))
    help_menu = HelpMenu(menubar)
    help_menu.extra_command(label=help_horoball_viewer_label, command=self.widget.help)
    help_menu.activate([help_horoball_viewer_label, help_report_bugs_label])
    self.menubar.add_cascade(label='Help', menu=help_menu)

def inside_view_menus(self):
    """Menus for the standalone Inside viewer.  Called by the view Frame, not the
    master Toplevel."""
    self.menubar = menubar = Tk_.Menu(self.master)
    Python_menu = Tk_.Menu(menubar, name="apple")
    Python_menu.add_command(label='About SnapPy ...',
                            command=lambda : about_snappy(self.master))
    Python_menu.add_separator()
    if sys.platform in ('linux2', 'linux') and self.main_window is not None:
        Python_menu.add_command(label='Quit SnapPy', command=
                                self.main_window.close)
    menubar.add_cascade(label='SnapPy', menu=Python_menu)
    File_menu = Tk_.Menu(menubar, name='file')
    add_menu(self.master, File_menu, 'Open...', None, 'disabled')
    add_menu(self.master, File_menu, 'Save as...', None, 'disabled')
    File_menu.add_command(label='Save Image...', command=self.master.save_image)
    File_menu.add_separator()
    add_menu(self.master, File_menu, 'Close', command=self.master.close)
    menubar.add_cascade(label='File', menu=File_menu)
    menubar.add_cascade(label='Edit ', menu=EditMenu(menubar, self.master.edit_actions))
    menubar.add_cascade(label='Window', menu=WindowMenu(menubar))
    help_menu = HelpMenu(menubar)
    #help_menu.extra_command(label=help_polyhedron_viewer_label, command=self.widget.help)
    #help_menu.activate([help_polyhedron_viewer_label, help_report_bugs_label])
    help_menu.activate([help_report_bugs_label])
    self.menubar.add_cascade(label='Help', menu=help_menu)
