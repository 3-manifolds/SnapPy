# -*- coding: utf-8 -*-
import os
import sys
import re
import time
from collections.abc import Mapping  # Python 3.5 or newer
from IPython.core.displayhook import DisplayHook
from tkinter.messagebox import askyesno

from .gui import *
from . import filedialog
from .exceptions import SnapPeaFatalError
from .app_menus import HelpMenu, EditMenu, WindowMenu, ListedWindow
from .app_menus import dirichlet_menus, horoball_menus, inside_view_menus, plink_menus
from .app_menus import add_menu, scut, open_html_docs
from .browser import Browser
from .horoviewer import HoroballViewer
from .infowindow import about_snappy, InfoWindow
from .polyviewer import PolyhedronViewer
from .raytracing.inside_viewer import InsideViewer
from .settings import Settings, SettingsDialog
from .phone_home import update_needed
from .SnapPy import SnapPea_interrupt, msg_stream
from .shell import SnapPyInteractiveShellEmbed
from .tkterminal import TkTerm, snappy_path

from plink import LinkEditor
from plink.smooth import Smoother

if 'SNAPPYHOME' in os.environ:
    if sys.platform == 'win32':
        os.environ['USERPROFILE'] = os.environ['SNAPPYHOME']
    else:
        os.environ['HOME'] = os.environ['SNAPPYHOME']


class SnapPyTerm(TkTerm, ListedWindow):
    """
    The main window of the SnapPy app, which runs an embedded IPython shell.
    """

    def __init__(self):
        self.ipython_shell = shell = SnapPyInteractiveShellEmbed.instance(
            banner1=app_banner + update_needed())
        shell.output = self
        shell.set_hook('show_in_pager', IPython_pager)
        self.main_window = self
        self.menu_title = 'SnapPy Shell'
        self.register_window(self)
        TkTerm.__init__(self, shell, name='SnapPy Command Shell')
        self.settings = SnapPySettings(self)
        self.start_interaction()
        if sys.platform == 'darwin':
            assert str(self.window) == "."
            # Under OS X, the main window shouldn't be closable.
            self.window.protocol('WM_DELETE_WINDOW', lambda : self.window.iconify())
            self.window.createcommand("::tk::mac::OpenDocument", self.OSX_open_filelist)
        else:
            self.window.tk.call('namespace', 'import', '::tk::dialog::file::')
            self.window.tk.call('set', '::tk::dialog::file::showHiddenBtn', '1')
            self.window.tk.call('set', '::tk::dialog::file::showHiddenVar', '0')
        self.encoding = None

    def add_bindings(self):
        self.window.bind('<<Paste>>', self.edit_paste)

    def about_window(self):
        window = self.window
        if not hasattr(window, 'about_snappy'):
            window.about_snappy = about_snappy(window)
        else:
            window.about_snappy.deiconify()
            window.about_snappy.lift()
            window.about_snappy.focus_force()

    def build_menus(self):
        window = self.window
        self.menubar = menubar = Tk_.Menu(window)
        Python_menu = Tk_.Menu(menubar, name="apple")
        Python_menu.add_command(label='About SnapPy...',
                                command=self.about_window)
        if sys.platform == 'darwin':
            window.createcommand('::tk::mac::ShowPreferences', self.edit_settings)
            # By default, the Quit menu command terminates the Python process.
            # That is not so friendly if cleanup is needed.
            window.createcommand('::tk::mac::Quit', self.close)
        else:
            Python_menu.add_separator()
            Python_menu.add_command(label='Settings...',
                                    command=self.edit_settings)
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
        if sys.platform == 'darwin':
            menubar.add_cascade(label='View', menu=Tk_.Menu(menubar, name='view'))
        menubar.add_cascade(label='Window', menu=WindowMenu(menubar))
        help_menu = HelpMenu(menubar)
        if sys.platform == 'darwin':
            window.createcommand('::tk::mac::ShowHelp', help_menu.show_SnapPy_help)
        menubar.add_cascade(label='Help', menu=help_menu)

    def edit_settings(self):
        terminal.can_quit = False
        if sys.platform == 'darwin':
            self.window.deletecommand('::tk::mac::ShowPreferences')
        else:
            apple_menu = self.menubar.children['apple']
            apple_menu.entryconfig(2, state='disabled')
        dialog = SettingsDialog(self.window, self.settings)
        terminal.add_blocker(dialog,
            'Changes to your settings will be lost if you quit SnapPy now.')
        dialog.run()
        terminal.remove_blocker(dialog)
        if dialog.okay:
            answer = askyesno('Save?',
                              'Do you want to save these settings?')
            if answer:
                self.settings.write_settings()
        if sys.platform == 'darwin':
            self.window.createcommand('::tk::mac::ShowPreferences', self.edit_settings)
        else:
            apple_menu.entryconfig(2, state='active')
        self.can_quit = True

    def OSX_open_filelist(self, *args):
        for arg in args:
            sys.stderr.write(repr(arg)+'\n')

    def open_file(self, event=None):
        openfile = filedialog.askopenfile(
            parent=self.window,
            title='Run Saved Transcript In Current Namespace',
            defaultextension='.py',
            filetypes=[
                ("Python and text files", "*.py *.ipy *.txt", "TEXT"),
                ("All text files", "", "TEXT"),
                ("All files", "")])
        if openfile:
            lines = openfile.readlines()
            openfile.close()
            if re.search(r"%\s*([vV]irtual)*\s*[lL]ink\s*[Pp]rojection", lines[0]):
                tkMessageBox.showwarning('Bad file',
                    'This is a SnapPea link projection file, '
                    'not a session transcript.')
            elif re.search(r"%\s*[tT]riangulation", lines[0]):
                tkMessageBox.showwarning('Bad file',
                    'This is a SnapPea triangulation file, '
                    'not a session transcript.')
            elif re.search(r"%\s*Generators", lines[0]):
                tkMessageBox.showwarning('Bad file',
                    'This is a SnapPea generator file, '
                    'not a session transcript.')
            else:
                while lines[0][0] in ('#', '\n'):
                    lines.pop(0)
                for line in lines:
                    # Skip comments.
                    if line.startswith('#'):
                        continue
                    # Simulate the user typing this line of code; strip off
                    # the indentation since that has already been printed by
                    # interact_prompt.
                    self.write(line[:-1].lstrip(), mark=Tk_.INSERT, advance=False)
                    # Then simulate the user pressing the Enter key.
                    self.handle_return(event=None)

    def open_link_file(self, event=None):
        openfile = filedialog.askopenfile(
            title='Load Link Projection File',
            defaultextension='.lnk',
            filetypes=[
                ("Link and text files", "*.lnk *.txt", "TEXT"),
                ("All text files", "", "TEXT"),
                ("All files", "")])
        if openfile:
            if not re.search(r"%\s*([vV]irtual)*\s*[lL]ink\s*[Pp]rojection", openfile.readline()):
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
            filetypes=[
                ("Python and text files", "*.py *.ipy *.txt", "TEXT"),
                ("All text files", "", "TEXT"),
                ("All files", "")])
        if savefile:
            savefile.write("""\
#!/usr/bin/env/python
# This script was saved by SnapPy on %s.
""" % time.asctime())
            inputs = self.IP.history_manager.input_hist_raw
            results = self.IP.history_manager.output_hist
            for n in range(1,len(inputs)):
                savefile.write('\n'+re.sub('\n+','\n',inputs[n]) + '\n')
                try:
                    output = repr(results[n]).split('\n')
                except:
                    continue
                for line in output:
                    savefile.write('#' + line + '\n')
            savefile.close()

    def save_file(self, event=None):
        self.window.bell()
        self.write2('Save As\n')

# These classes assume that the global variable "terminal" exists


class SnapPyBrowser(Browser, ListedWindow):
    def __init__(self, manifold, root=None, main_window=None):
        Browser.__init__(self, manifold, root=root, main_window=terminal)
        self.settings = terminal.settings
        self.menu_title = self.title()
        self.register_window(self)
        self.dirichlet_viewer.help_button.configure(command=self.dirichlet_help)

    def close(self, event=None):
        self.unregister_window(self)
        self.destroy()

    def apply_settings(self):
        if self.inside_view:
            self.inside_view.apply_settings(self.main_window.settings)

    def dirichlet_help(self):
        if not hasattr(self, 'polyhedron_help'):
            self.polyhedron_help = InfoWindow(self, 'Polyhedron Viewer Help',
                                              self.dirichlet_viewer.widget.help_text,
                                              'polyhedron_help')
        else:
            self.polyhedron_help.deiconify()
            self.polyhedron_help.lift()
            self.polyhedron_help.focus_force()

    def horoball_help(self):
        if not hasattr(self, 'horoviewer_help'):
            self.horoviewer_help = InfoWindow(self, 'Horoball Viewer Help',
                                              self.horoball_viewer.widget.help_text,
                                              'horoviewer_help')
        else:
            self.horoviewer_help.deiconify()
            self.horoviewer_help.lift()
            self.horoviewer_help.focus_force()


class SnapPyLinkEditor(LinkEditor, ListedWindow):
    def __init__(self, root=None, no_arcs=False, callback=None, cb_menu='',
                 manifold=None, file_name=None):
        self.manifold = manifold
        self.main_window = terminal
        LinkEditor.__init__(self, root=terminal.window, no_arcs=no_arcs,
                            callback=callback, cb_menu=cb_menu,
                            manifold=manifold, file_name=file_name)
        self.set_title()
        self.register_window(self)
        self.window.focus_set()
        self.window.after_idle(self.set_title)

    def done(self, event=None):
        self.unregister_window(self)
        self.window.withdraw()

    def reopen(self):
        self.register_window(self)
        self.window.deiconify()

    def deiconify(self):
        self.window.deiconify()

    def lift(self):
        self.window.lift()

    def focus_force(self):
        self.window.focus_force()

    def set_title(self):
        # Try to determine the variable associated to the manifold:
        title = 'Plink Editor'
        if self.IP:
            ns = self.IP.user_ns
            names = [name for name in ns
                     if ns[name] is self.manifold]
            if names:
                names.sort(key=lambda x : '}'+x if x.startswith('_') else x)
                if names[0] == '_':
                    count = self.IP.execution_count
                    title += ' - Out[%d]' % count
                else:
                    title += ' - %s' % names[0]
            else:
                count = self.IP.execution_count
                if ns['_'] is self.manifold:
                    title += ' - Out[%d]' % count
        self.window.title(title)
        self.menu_title = title

    _build_menus = plink_menus

    def load(self, event=None, file_name=None):
        LinkEditor.load(self, file_name)

    def save(self, event=None):
        LinkEditor.save(self)

    def howto(self):
        open_html_docs('plink.html')

    __repr__ = object.__repr__


class SnapPyViewerWindow(ViewerWindow, ListedWindow):
    def __init__(self, *args, **kwargs):
        ViewerWindow.__init__(self, *args, **kwargs)
        self.main_window = terminal
        self.menu_title = self.title()
        self.register_window(self)

    def apply_settings(self):
        # The view's apply_settings method has a different signature.  It
        # expects to be passed a Settings object.
        self.view.apply_settings(self.main_window.settings)

    def close(self, event=None):
        self.unregister_window(self)
        self.view = None
        self.destroy()


class SnapPyPolyhedronViewer(PolyhedronViewer):

    build_menus = dirichlet_menus

    def __init__(self, *args, **kwargs):
        PolyhedronViewer.__init__(self, *args, **kwargs, main_window=terminal)
        self.help_button.configure(command=self.help_window)

    def help_window(self):
        window = self.parent
        if not hasattr(window, 'polyhedron_help'):
            window.polyhedron_help = InfoWindow(window,  'Polyhedron Viewer Help',
                                                self.widget.help_text, 'polyhedron_help')
        else:
            window.polyhedron_help.deiconify()
            window.polyhedron_help.lift()
            window.polyhedron_help.focus_force()


class SnapPyHoroballViewer(HoroballViewer):

    build_menus = horoball_menus

    def __init__(self, *args, **kwargs):
        HoroballViewer.__init__(self, *args, **kwargs, main_window=terminal)
        self.main_window = terminal

    def help_window(self):
        window = self.parent
        if not hasattr(window, 'horoball_help'):
            window.horoball_help = InfoWindow(window,  'Horoball Viewer Help',
                                              self.widget.help_text, 'horoball_help')
        else:
            window.horoball_help.deiconify()
            window.horoball_help.lift()
            window.horoball_help.focus_force()


class SnapPyInsideViewer(InsideViewer):

    build_menus = inside_view_menus


class SnapPySettings(Settings, ListedWindow):
    def __init__(self, terminal):
        self.terminal = terminal
        Settings.__init__(self, terminal.text)
        self.apply_settings()

    def apply_settings(self):
        self.terminal.set_font(self['font'])
        changed = self.changed()
        IP = self.terminal.IP
        self.terminal.quiet = True
        if 'autocall' in changed:
            if self.settings_dict['autocall']:
                IP.magics_manager.magics['line']['autocall'](2)
            else:
                IP.magics_manager.magics['line']['autocall'](0)
        if 'automagic' in changed:
            if self.settings_dict['automagic']:
                IP.magics_manager.magics['line']['automagic']('on')
            else:
                IP.magics_manager.magics['line']['automagic']('off')
        self.terminal.quiet = False
        for window in self.window_list:
            window.apply_settings()


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


class _Helper():
    pass


import site
site._Helper = _Helper

# This will be used for paging by IPython help.


def IPython_pager(self, text, start=0, screen_lines=0):
    if isinstance(text, Mapping):
        text = text['text/plain']
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
            import snappy
            ico = os.path.join(os.path.dirname(snappy.__file__), 'SnapPy.ico')
            window.iconbitmap(default=ico)
        except:
            pass
    if sys.platform == 'darwin':
        if not sys.executable.endswith('SnapPy.app/Contents/MacOS/python'):
            image_file = os.path.join(snappy_path, 'SnapPy.png')
            if os.path.exists(image_file):
                dock_icon = Tk_.PhotoImage(file=image_file)
                window.eval('wm iconphoto . -default %s' % dock_icon)


# from multiprocessing import Process
class SnapPyKernelServer():
    """
    Placeholder for a real SnapPyKernelServer which the app can use to
    asynchronously compute data about manifolds.
    """
    def __init__(self):
        self._process = Process(target=self.task)
        self._process.start()

    def task(self):
        while 1:
            # uncomment these two lines to watch the fake process running
            # print('.', end='')
            # sys.stdout.flush()
            time.sleep(2)

    def stop(self):
        self._process.terminate()

    def __del__(self):
        if self._process.is_alive():
            self._process.terminate()


def main():
    global terminal
    import snappy
    # kernel_server = SnapPyKernelServer()
    terminal = SnapPyTerm()
    sys.stdout = terminal
    set_icon(terminal.window)
    SnapPy_ns = dict([(x, getattr(snappy, x)) for x in snappy.__all__])
    # SnapPy_ns['kernel_server'] = kernel_server
    SnapPy_ns['exit'] = SnapPy_ns['quit'] = SnapPyExit()
    SnapPy_ns['pager'] = None
    helper = pydoc.Helper(input=terminal, output=terminal)
    helper.__call__ = lambda x=None : helper.help(x) if x else SnapPy_help()
    helper.__repr__ = lambda : help_banner
    SnapPy_ns['help'] = helper
    terminal.ipython_shell.user_ns.update(SnapPy_ns)
    snappy.browser.main_window = terminal
    LP, HP = snappy.SnapPy, snappy.SnapPyHP
    LP.LinkEditor = HP.LinkEditor = SnapPyLinkEditor
    SnapPyLinkEditor.IP = terminal.ipython_shell
    LP.ViewerWindow = HP.ViewerWindow = SnapPyViewerWindow
    LP.PolyhedronViewer = HP.PolyhedronViewer = SnapPyPolyhedronViewer
    LP.HoroballViewer = HP.HoroballViewer = SnapPyHoroballViewer
    snappy.ViewerWindow = SnapPyViewerWindow
    snappy.ViewerWindow.main_window = terminal
    snappy.InsideViewer = SnapPyInsideViewer
    snappy.InsideViewer.main_window = terminal
    LP.Browser = HP.Browser = SnapPyBrowser
    LP.msg_stream.write = HP.msg_stream.write = terminal.write2
    LP.UI_callback = HP.UI_callback = terminal.SnapPea_callback
    # if not snappy.SnapPy._within_sage:
    #    snappy.pari.UI_callback = terminal.PARI_callback
    terminal.window.lift()
    terminal.window.mainloop()
    # kernel_server.stop()


if __name__ == "__main__":
    main()
