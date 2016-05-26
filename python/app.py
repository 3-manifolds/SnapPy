# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import os, sys, re, tempfile, time, png, webbrowser

# An IPython shell to use in the terminal.
try:
    from IPython.terminal.embed import InteractiveShellEmbed
except ImportError:
    from IPython.frontend.terminal.embed import InteractiveShellEmbed
InteractiveShellEmbed.readline_use = False
InteractiveShellEmbed.autoindent = False
# Temporary hack to enable colors without readline.
# IPython says this will be removed in the future.
InteractiveShellEmbed.colors_force = True

import snappy
from snappy.tkterminal import TkTerm
from snappy.app_menus import HelpMenu, EditMenu, WindowMenu
from snappy.app_menus import dirichlet_menus, horoball_menus, plink_menus
from snappy.app_menus import togl_save_image, add_menu, scut
from snappy import filedialog
from snappy import SnapPeaFatalError
from snappy.polyviewer import PolyhedronViewer
from snappy.horoviewer import HoroballViewer
from snappy.browser import Browser
from snappy.SnapPy import SnapPea_interrupt, msg_stream
from snappy.preferences import Preferences, PreferenceDialog
from snappy.infodialog import about_snappy
from snappy.phone_home import update_needed

try:
    import Tkinter as Tk_
    import tkMessageBox
    from tkMessageBox import askyesno
    from tkFont import Font

except ImportError: # Python 3
    import tkinter as Tk_
    import tkinter.messagebox as tkMessageBox
    from tkinter.messagebox import askyesno 
    from tkinter.font import Font

from plink import LinkEditor
from plink.smooth import Smoother

class SnapPyTerm(TkTerm, WindowMenu):

    def __init__(self, the_shell):
        self.main_window = self
        self.menu_title = 'SnapPy Shell'
        WindowMenu.register(self)
        TkTerm.__init__(self, the_shell, name='SnapPy Command Shell')
        self.prefs = SnapPyPreferences(self)
        self.start_interaction()
        if sys.platform == 'darwin':
            assert str(self.window) == "."
            # Under OS X, the main window shouldn't be closable.
            self.window.protocol('WM_DELETE_WINDOW', lambda : self.window.iconify())
            self.window.createcommand("::tk::mac::OpenDocument", self.OSX_open_filelist)
        else:
            self.window.tk.call('namespace', 'import', '::tk::dialog::file::')
            self.window.tk.call('set', '::tk::dialog::file::showHiddenBtn',  '1')
            self.window.tk.call('set', '::tk::dialog::file::showHiddenVar',  '0')

    def add_bindings(self):
        self.window.bind('<<Paste>>', self.edit_paste)

    def build_menus(self):
        window = self.window
        self.menubar = menubar = Tk_.Menu(window)
        Python_menu = Tk_.Menu(menubar, name="apple")
        Python_menu.add_command(label='About SnapPy...',
                                command=lambda : about_snappy(window))
        if sys.platform == 'darwin':
            window.createcommand('::tk::mac::ShowPreferences', self.edit_prefs)
        else:
            Python_menu.add_separator()
            Python_menu.add_command(label='Preferences...',
                                    command=self.edit_prefs)
            Python_menu.add_separator()
            Python_menu.add_command(label='Quit SnapPy', command=self.close)
        menubar.add_cascade(label='SnapPy', menu=Python_menu)
        File_menu = Tk_.Menu(menubar, name='file')
        add_menu(window, File_menu, 'Open...', self.open_file)
        add_menu(window, File_menu, 'Open link...', self.open_link_file)
        add_menu(window, File_menu, 'Save', self.save_file, state='disabled')
        add_menu(window, File_menu, 'Save as...', self.save_file_as)
        menubar.add_cascade(label='File', menu=File_menu)
        menubar.add_cascade(label='Edit ', menu=EditMenu(menubar, self.edit_actions))
        menubar.add_cascade(label='Window', menu=WindowMenu(menubar))
        help_menu = HelpMenu(menubar)
        if sys.platform == 'darwin':
            window.createcommand('::tk::mac::ShowHelp', help_menu.show_SnapPy_help)
        menubar.add_cascade(label='Help', menu=help_menu)
        
    def edit_prefs(self):
        apple_menu = self.menubar.children['apple']
        apple_menu.entryconfig(2, state='disabled')
        dialog = PreferenceDialog(self.window, self.prefs)
        if dialog.okay:
            answer = askyesno('Save?',
                              'Do you want to save these settings?')
            if answer:
                self.prefs.write_prefs()
        apple_menu.entryconfig(2, state='active')

    def OSX_open_filelist(self, *args):
        for arg in args:
            sys.stderr.write(repr(arg)+'\n')

    def open_file(self, event=None):
        openfile = filedialog.askopenfile(
            parent=self.window,
            title='Run Saved Transcript In Current Namespace',
            defaultextension='.py',
            filetypes = [
                ("Python and text files", "*.py *.ipy *.txt", "TEXT"),
                ("All text files", "", "TEXT"),
                ("All files", "")])
        if openfile:
            lines = openfile.readlines()
            openfile.close()
            if re.search("%\s*([vV]irtual)*\s*[lL]ink\s*[Pp]rojection", lines[0]):
                tkMessageBox.showwarning('Bad file',
                                         'This is a SnapPea link projection file, not a session transcript.')
            elif re.search("%\s*[tT]riangulation", lines[0]):
                tkMessageBox.showwarning('Bad file',
                                         'This is a SnapPea triangulation file, not a session transcript.')
            elif re.search("%\s*Generators", lines[0]):
                tkMessageBox.showwarning('Bad file',
                                         'This is a SnapPea generator file, not a session transcript.')
            else:
                for line in lines:
                    if line.startswith('#') or len(line) == 1:
                        continue
                    self.write(line)
                    self.interact_handle_input(line)
                    self.interact_prompt()

    def open_link_file(self):
        openfile = filedialog.askopenfile(
            parent=self.window,
            title='Load Link Projection File',
            defaultextension='.lnk',
            filetypes = [
                ("Link and text files", "*.lnk *.txt", "TEXT"),
                ("All text files", "", "TEXT"),
                ("All files", "")])
        if openfile:
            if not re.search("%\s*([vV]irtual)*\s*[lL]ink\s*[Pp]rojection", openfile.readline()):
                tkMessageBox.showwarning('Bad file',
                                         'This is not a SnapPea link projection file')
                openfile.close()
            else:
                name = openfile.name
                openfile.close()
                line = "Manifold()\n"
                self.write(line)
                self.interact_handle_input(line)
                self.interact_prompt()
                M = self.IP.user_ns['_']
                M.LE.load(file_name=name)

    def save_file_as(self, event=None):
        savefile = filedialog.asksaveasfile(
            parent=self.window,
            mode='w',
            title='Save Transcript as a Python script',
            defaultextension='.py',
            filetypes = [
                ("Python and text files", "*.py *.ipy *.txt", "TEXT"),
                ("All text files", "", "TEXT"),
                ("All files", "")])
        if savefile:
            savefile.write("""\
#!/usr/bin/env/python
# This script was saved by SnapPy on %s.
"""%time.asctime())
            inputs = self.IP.history_manager.input_hist_raw
            results = self.IP.history_manager.output_hist
            for n in range(1,len(inputs)):
                savefile.write('\n'+re.sub('\n+','\n',inputs[n]) +'\n')
                try:
                    output = repr(results[n]).split('\n')
                except KeyError:
                    continue
                for line in output:
                    savefile.write('#' + line + '\n')
            savefile.close()

    def save_file(self, event=None):
        self.window.bell()
        self.write2('Save As\n')

# These classes assume that the global variable "terminal" exists

class SnapPyBrowser(Browser, WindowMenu):
    def __init__(self, manifold):
        Browser.__init__(self, manifold, terminal.window)
        self.prefs = terminal.prefs
        self.menu_title = self.window.title()
        WindowMenu.register(self)
        self.main_window = terminal

    def close(self, event=None):
        WindowMenu.unregister(self)
        self.window.destroy()

class SnapPyLinkEditor(LinkEditor, WindowMenu):
    def __init__(self, root=None, no_arcs=False, callback=None, cb_menu='',
                 manifold=None, file_name=None):
        self.manifold = manifold
        self.main_window = terminal
        LinkEditor.__init__(self, root=terminal.window, no_arcs=no_arcs,
                            callback=callback, cb_menu=cb_menu,
                            manifold=manifold, file_name=file_name)
        self.set_title()
        WindowMenu.register(self)
        self.window.focus_set()
        self.window.update_idletasks()
        self.window.after_idle(self.set_title)

    def set_title(self):
        # Try to determine the variable associated to the manifold:
        title = 'Plink Editor'
        if self.IP:
            ns = self.IP.user_ns
            names = [name for name in ns
                     if ns[name] is self.manifold]
            if names:
                names.sort(key=lambda x : '}'+x if x.startswith('_') else x)
                title += ' - %s' % names[0]
            else:
                count = self.IP.execution_count
                if ns['_'] is self.manifold:
                    title += ' - Out[%d]'%count
        self.window.title(title)
        self.menu_title = title

    build_menus = plink_menus
 
    def copy_info(self):
        if not self.infotext.selection_present():
           self.infotext.selection_range(0, Tk_.END)
        self.infotext.focus()
        self.infotext.event_generate('<<Copy>>')

    def load(self, event=None, file_name=None):
        LinkEditor.load(self, file_name)

    def save(self, event=None):
        LinkEditor.save(self)

class SnapPyPolyhedronViewer(PolyhedronViewer, WindowMenu):
    def __init__(self, facedicts, root=None, title='Polyhedron Viewer'):
        self.main_window = terminal
        PolyhedronViewer.__init__(self, facedicts, root=terminal.window,
                                  title=title)
        self.menu_title = self.window.title()
        WindowMenu.register(self)

    def add_help(self):
        pass

    build_menus = dirichlet_menus

    def edit_actions(self):
        return {}
    
    def close(self):
        self.polyhedron.destroy()
        WindowMenu.unregister(self)
        self.window.destroy()

    def save_image(self):
        togl_save_image(self)

class SnapPyHoroballViewer(HoroballViewer, WindowMenu):
    def __init__(self, nbhd, which_cusp=0, cutoff=None,
                 root=None, title='Horoball Viewer'):
        self.main_window = terminal
        HoroballViewer.__init__(self, nbhd, which_cusp=which_cusp,
                                cutoff=cutoff, root=terminal.window,
                                title=title, prefs = terminal.prefs)
        self.menu_title = self.window.title()
        WindowMenu.register(self)
        self.view_check()

    build_menus = horoball_menus

    def add_help(self):
        pass

    def edit_actions(self):
        return {}

    def close(self):
        self.widget.activate()
        self.scene.destroy()
        WindowMenu.unregister(self)
        self.window.destroy()

    def save_image(self):
        togl_save_image(self)

class SnapPyPreferences(Preferences):
    def __init__(self, terminal):
        self.terminal = terminal
        Preferences.__init__(self, terminal.text)
        self.apply_prefs()

    def apply_prefs(self):
        self.terminal.set_font(self['font'])
        self.terminal.window.update_idletasks()
        changed = self.changed()
        IP = self.terminal.IP
        self.terminal.quiet = True
        if 'autocall' in changed:
            if self.prefs_dict['autocall']:
                IP.magics_manager.magics['line']['autocall'](2)
            else:
                IP.magics_manager.magics['line']['autocall'](0)
        if 'automagic' in changed:
            if self.prefs_dict['automagic']:
                IP.magics_manager.magics['line']['automagic']('on')
            else:
                IP.magics_manager.magics['line']['automagic']('off')
        self.terminal.quiet = False

app_banner = """
 Hi.  It's SnapPy.  
 SnapPy is based on the SnapPea kernel, written by Jeff Weeks.
 Type "Manifold?" to get started.
"""

help_banner = """Type X? for help with X.
Use the Help menu or type help() to view the SnapPy documentation."""

class SnapPyExit:
    """
    Replacement for the IPython ExitAutocall class
    """
    def __repr__(self):
        return 'Please use the SnapPy menu to quit.'
    __str__ = __repr__

    def __call__(self):
        return self

# This hack avoids an unnecessary warning from IPython saying that
# _Helper is not included in the app2py site.py file.
class _Helper(object):
    pass
import site
site._Helper = _Helper

# This will be used for paging by IPython help.
def IPython_pager(self, text, start=0, screen_lines=0):
    terminal.page(text)

# This will be used for paging by pydoc help.
import pydoc

def pydoc_pager(text):
    terminal.page(pydoc.plain(text))

pydoc.getpager() # this call creates the global variable pydoc.pager
pydoc.pager = pydoc_pager

# This sets the "system menu" icon in the title bar to be the SnapPy
# icon (in Windows and ??KDE??)

def set_icon(window):
    if sys.platform == 'win32':
        try:
            ico = os.path.join(os.path.dirname(snappy.__file__), 'SnapPy.ico')
            window.iconbitmap(default=ico)
        except:
            pass

def main():
    global terminal
    the_shell = InteractiveShellEmbed.instance(
        banner1=app_banner + update_needed())
    terminal = SnapPyTerm(the_shell)
    the_shell.tkterm = terminal
    set_icon(terminal.window)
    the_shell.set_hook('show_in_pager', IPython_pager)
    SnapPy_ns = dict([(x, getattr(snappy,x)) for x in snappy.__all__])
    SnapPy_ns['exit'] = SnapPy_ns['quit'] = SnapPyExit()
    SnapPy_ns['pager'] = None
    helper = pydoc.Helper(input=terminal, output=terminal)
    helper.__call__ = lambda x=None : helper.help(x) if x else SnapPy_help()
    helper.__repr__ = lambda : help_banner
    SnapPy_ns['help'] = helper
    the_shell.user_ns.update(SnapPy_ns)
    snappy.browser.main_window = terminal
    LP, HP = snappy.SnapPy, snappy.SnapPyHP
    LP.LinkEditor = HP.LinkEditor = SnapPyLinkEditor
    SnapPyLinkEditor.IP = the_shell
    LP.PolyhedronViewer = HP.PolyhedronViewer = SnapPyPolyhedronViewer
    LP.HoroballViewer = HP.HoroballViewer = SnapPyHoroballViewer
    LP.Browser = HP.Browser = SnapPyBrowser
    LP.msg_stream.write = HP.msg_stream.write = terminal.write2
    LP.UI_callback = HP.UI_callback = terminal.SnapPea_callback
    if not snappy.SnapPy._within_sage:
        snappy.pari.UI_callback = terminal.PARI_callback
    terminal.window.lift()
    terminal.window.mainloop()

if __name__ == "__main__":
    main()
