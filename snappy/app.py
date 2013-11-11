# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import os, sys, IPython, datetime
try:
    from IPython.terminal.embed import InteractiveShellEmbed
except ImportError:
    from IPython.frontend.terminal.embed import InteractiveShellEmbed
InteractiveShellEmbed.readline_use = False
InteractiveShellEmbed.autoindent = False
# Temporary hack to enable colors without readline.
# IPython says this will be removed in the future.
InteractiveShellEmbed.colors_force = True
from IPython.utils import io
from IPython.core.autocall import IPyAutocall
import snappy
from snappy.app_menus import dirichlet_menus, horoball_menus, really_disable_menu_items
from snappy.app_menus import togl_save_image, add_menu, scut
from snappy import filedialog
from snappy import SnapPeaFatalError
from snappy.polyviewer import PolyhedronViewer
from snappy.horoviewer import HoroballViewer
from snappy.browser import Browser
from snappy.SnapPy import SnapPea_interrupt, msg_stream
from snappy.preferences import Preferences, PreferenceDialog
from snappy.infodialog import InfoDialog
from snappy.version import version as SnapPy_version
from snappy.phone_home import update_needed
snappy_path = os.path.dirname(snappy.__file__)
icon_file = os.path.join(snappy_path, 'info_icon.gif')
if sys.platform == 'win32':
    mousewheel_factor = -120
else:
    mousewheel_factor = 1

try:
    import Tkinter as Tk_
    import tkMessageBox
    from tkFont import Font
    from tkMessageBox import askyesno
    from urllib import pathname2url

except ImportError: # Python 3
    import tkinter as Tk_
    import tkinter.messagebox as tkMessageBox
    from tkinter.font import Font
    from tkinter.messagebox import askyesno 
    from urllib.request import pathname2url

import os, sys, re, webbrowser, signal, tempfile, time, png
from plink import LinkEditor
from plink.smooth import Smoother

debug_Tk = False
    
ansi_seqs = re.compile('(?:\x01*\x1b\[((?:[0-9]*;)*[0-9]*.)\x02*)*([^\x01\x1b]*)',
                       re.MULTILINE)

ansi_colors =  {'0;30m': 'Black',
                '0;31m': 'Red',
                '0;32m': 'Green',
                '0;33m': 'Brown',
                '0;34m': 'Blue',
                '0;35m': 'Purple',
                '0;36m': 'Cyan',
                '0;37m': 'LightGray',
                '1;30m': 'Black', #'DarkGray',
                '1;31m': 'DarkRed',
                '1;32m': 'SeaGreen',
                '1;33m': 'Yellow',
                '1;34m': 'LightBlue',
                '1;35m': 'MediumPurple',
                '1;36m': 'LightCyan',
                '1;37m': 'White'}

delims = re.compile(r'[\s\[\]\{\}\(\)\+\-\=\'`~!@#\$\^\&\*]+')

about_snappy = """
For information on how to use SnapPy, please see the
Help menu.

SnapPy is a user interface for the SnapPea kernel,
which was written by Jeff Weeks.  SnapPy was
written by Marc Culler and Nathan Dunfield and
is distributed under the GNU Public License,
version 2 or later.  Its home page is:
     http://snappy.computop.org/

The release number of this SnapPy is %s.

SnapPy is written in the Python language, using
Cython to incorporate the SnapPea kernel code.
The graphical interface uses Tcl/Tk, via Python's
Tkinter module.

Information, downloads, and source code for the
SnapPea kernel and for the user interfaces written
by Jeff Weeks are available at:
     http://www.geometrygames.org/SnapPea-old/
     http://www.geometrygames.org/SnapPea/

Copyright © 2009-%d, Marc Culler, Nathan
Dunfield, and others.
"""% (SnapPy_version, datetime.datetime.now().year)

class Tk(Tk_.Tk):
    def __init__(self, error_handler=None):
        Tk_.Tk.__init__(self, className='snappy')
        # Tkinter ignores exceptions raised by callbacks, but
        # calls this function to report their occurence.
        if error_handler:
            self.report_callback_exception = error_handler

class TkTerm:
    """
    A Tkinter terminal window that runs an IPython shell.
    Supports the IOStream interface, and can function as
    a replacement for sys.stdout / sys.stderr.
    Some ideas borrowed from code written by Eitan Isaacson, IBM Corp.
    """
    def __init__(self, the_shell, name='TkTerm'):
        if debug_Tk:
            self.window = window = Tk()
        else:
            self.window = window = Tk(self.report_callback_exception)
            io.stdout = self
            io.stderr = self
        # Prevent runt window
        try:
            window.tk.call('console', 'hide')
        except Tk_.TclError:
            pass
        # Global Tk options
        window.option_add('*Menu.tearOff', 0)
        # Construct the window
        window.title(name)
        window.protocol("WM_DELETE_WINDOW", self.close)
        self.icon = Tk_.PhotoImage(file=icon_file)
        self.frame = frame = Tk_.Frame(window)
        self.text = text = Tk_.Text(frame,
                                    foreground='Black',
                                    borderwidth=3,
                                    background='#ec0fffec0',
                                    highlightthickness=0,
                                    relief=Tk_.FLAT
                                )
        self.scroller = scroller = Tk_.Scrollbar(frame, command=text.yview)
        text.config(yscrollcommand=scroller.set)
        scroller.pack(side=Tk_.RIGHT, fill=Tk_.Y, pady=10)
        text.pack(fill=Tk_.BOTH, expand=Tk_.YES)
        frame.pack(fill=Tk_.BOTH, expand=Tk_.YES)
        text.focus_set()
        text.bind('<KeyPress>', self.handle_keypress)
        text.bind('<Return>', self.handle_return)
        text.bind('<Shift-Return>', lambda event : None)
        text.bind('<BackSpace>', self.handle_backspace)
        text.bind('<Delete>', self.handle_backspace)
        text.bind('<Tab>', self.handle_tab)
        text.bind('<Up>', self.handle_up)
        text.bind('<Shift-Up>', lambda event : None)
        text.bind('<Down>', self.handle_down)
        text.bind('<Shift-Down>', lambda event : None)
        text.bind('<<Cut>>', self.protect_text)
        text.bind('<<Paste>>', self.paste)
        text.bind('<<Clear>>', self.protect_text)
        if sys.platform != 'darwin':
            text.bind_all('<ButtonPress-2>', self.middle_mouse_down)
            text.bind_all('<ButtonRelease-2>', self.middle_mouse_up)
        text.bind('<Button-3>', lambda event:'break')
        text.bind('<Button-4>', lambda event:text.yview_scroll(-1, Tk_.UNITS))
        text.bind('<Button-5>', lambda event:text.yview_scroll(1, Tk_.UNITS))
        text.bind('<MouseWheel>', self.handle_mousewheel)
        if sys.platform == 'darwin':
            self.window.bind_all('<Command-Key-q>', self.close)
        elif sys.platform == 'linux2':
            self.window.bind_all('<Alt-Key-q>', self.close)
        self.add_bindings()
        # 'output_end' marks the end of the text written by us.
        # Everything above this position should be
        # immutable, and tagged with the "output" style.
        self.text.mark_set('output_end', Tk_.INSERT)
        self.text.mark_gravity('output_end', Tk_.LEFT)
        text.tag_config('output')
        # Make sure we don't override the cut-paste background.
        text.tag_lower('output') 
        # Remember where we were when tab was pressed.
        self.tab_index = None
        self.tab_count = 0
        # Manage history
        self.hist_pointer = 0
        self.hist_stem = ''
        self.editing_hist = False
        self.filtered_hist = []
        # Remember illegal pastes
        self.nasty = None
        self.nasty_text = None
        # Build style tags for colored text, 
        for code in ansi_colors:
            text.tag_config(code, foreground=ansi_colors[code])
        # and a style tag for messages.
        text.tag_config('msg', foreground='Red')
        self.build_menus()
        self.window.config(menu=self.menubar)
        self.edit_config(None)
        self.output_count = 0
        # pager support
        self.prompt_index = None
        self.scroll_back = False
        # Setup the IPython embedded shell
        self.IP = the_shell
        if not debug_Tk:
            # write and write_err will be removed soon
            the_shell.write = self.write
            the_shell.write_err = self.write
            the_shell._showtraceback = self.showtraceback
        the_shell.system = self.system
        sys.displayhook = the_shell.displayhook
        the_shell.more = False
        the_shell.magics_manager.magics['line']['colors']('LightBG')
        if the_shell.banner1:
            self.banner = the_shell.banner1
        else:
            cprt = 'Type "copyright", "credits" or "license" for more information.'
            self.banner = "Python %s on %s\n%s\n(%s)\n" %(
                sys.version, sys.platform, cprt,
                self.__class__.__name__)
        self.quiet = False
        # We will set this flag if the user types ^C
        self.interrupted = False
        # This flag will be set when IPython is running code
        self.running_code = False
        # This flag will be set when a SnapPea computation is aborted
        self.aborted_SnapPea = False
        # Let the UI update itself (and check for ^C) every second.
        if sys.platform != 'win32':
            signal.signal(signal.SIGALRM, lambda sig, frame : None )
            try:
                signal.setitimer(signal.ITIMER_REAL, 1.0, 1.0)
            except AttributeError: # itimer is not supported in python 2.5
                pass
        self.closed = False
        #self.start_interaction()

    # For subclasses to override:
    def build_menus(self):
        pass

    # For subclasses to override:
    def add_bindings(self):
        pass

    def system(self, cmd):
        output = self.IP.getoutput(cmd, split=False)
        self.write(output)
        
    def showtraceback(self, etype, evalue, stb):
        self.write(self.IP.InteractiveTB.stb2text(stb))

    def SnapPea_callback(self, interrupted=False):
        """
        Callback for SnapPea to keep the UI alive during long computations.
        """
        self.window.update()
        if interrupted:
            self.interrupted = False
            raise KeyboardInterrupt('SnapPea computation aborted')

    def PARI_callback(self):
        """
        Callback for CyPari to keep the UI alive during long computations.
        """
        self.window.update()
        if self.running_code:
            if self.interrupted:
                snappy.pari.abort()
                self.interrupted = False
                raise KeyboardInterrupt('PARI computation aborted')

    def interrupt(self):
        # Tell the ticker to raise a KeyboardInterrupt after the update
        if not self.running_code:
            raise KeyboardInterrupt('Halted')
        else:
            self.interrupted = True
            # Inform the SnapPea kernel about the interrupt.
            SnapPea_interrupt()
            
    def report_callback_exception(self, exc, value, traceback):
        # This is called when exceptions are caught by Tk.
#        self.write2('(Tk) ' + exc.__name__ +': ' + str(value) +'\n')
        sys.last_type = exc
        sys.last_value = value
        sys.last_traceback = traceback
        self.IP.showtraceback()
    
    def set_font(self, fontdesc):
        self.text.config(font=fontdesc)
        self.char_size = Font(self.text, fontdesc).measure('M')
        self.text.tag_add('all', '1.0', Tk_.END)
        self.text.tag_config('all', font=fontdesc)

    def close(self, event=None):
        self.window.update_idletasks()
        self.window.quit()
        self.closed = True

    def handle_mousewheel(self, event):
        # OS X scroll gestures are smoother if handled by the
        # scrollbar.
        self.window.event_generate(
            '<MouseWheel>',
            delta=-event.delta,
            state=event.state,
            rootx=event.x_root,
            rooty=event.y_root,
            x=event.x, y=event.y,
            time=event.time,
            root=self.scroller)
        # try to synchronize
        self.text.yview_scroll(0, Tk_.UNITS)
        return
    
    def handle_keypress(self, event):
        self.clear_completions()
        # OS X Tk > 8.4 sends weird strings for some keys 
        if len(event.char) != 1:
            return
        if event.char == '\001': # ^A
            self.text.mark_set(Tk_.INSERT, 'output_end')
            return 'break'
        if event.char == '\025': # ^U
            self.text.delete('output_end', Tk_.END)
            return 'break'
        if event.char == '\040': # space
            if self.text.compare(Tk_.INSERT, '<', 'output_end'):
                self.page_down()
                return 'break'
        if event.char == '\003': # ^C
            self.interrupt()
        if self.text.compare(Tk_.INSERT, '<', 'output_end'):
            self.text.mark_set(Tk_.INSERT, 'output_end')

    def handle_return(self, event):
        self.clear_completions()
        line=self.text.get('output_end', Tk_.END)
        self.text.tag_add('output', 'output_end', Tk_.END)
        self.text.mark_set('output_end', Tk_.END)
        if not self.running_code:
            self.send_line(line)
        return 'break'

    def handle_backspace(self, event):
        self.clear_completions()
        if self.text.compare(Tk_.INSERT, '<=', 'output_end'):
            self.window.bell()
            return 'break'
        if self.IP.indent_current_nsp >= 4:
            if self.text.get(Tk_.INSERT+'-4c', Tk_.INSERT) == '    ':
                self.text.delete(Tk_.INSERT+'-4c', Tk_.INSERT)
                return 'break'

    def handle_tab(self, event):
        if self.text.compare(Tk_.INSERT, '<', 'output_end'):
            return 'break'
        self.tab_index = self.text.index(Tk_.INSERT)
        self.tab_count += 1
        if self.tab_count > 2:
            self.clear_completions()
            return 'break'
        line = self.text.get('output_end', self.tab_index).strip('\n')
        word = delims.split(line)[-1]
        try:
            completions = self.IP.complete(word)[1]
        except TypeError:
            completions = []
        if word.find('_') == -1:
            completions = [x for x in completions 
                           if x.find('__') == -1 and x.find('._') == -1]
        if len(completions) == 0:
            self.window.bell()
            self.tab_count = 0
            return 'break'
        stem = self.stem(completions)
        if len(stem) > len(word):
            self.do_completion(word, stem)
        elif len(completions) > 60 and self.tab_count == 1:
            self.show_completions(['%s possibilities -- hit tab again to view them all'%len(completions)])
        else:
            self.show_completions(completions)
            if len(completions) <= 60:
                self.tab_count += 1
        return 'break'

    def do_completion(self, word, completion):
        tail = completion[len(word):]
        self.text.insert(self.tab_index, tail)
        self.text.delete(Tk_.INSERT, Tk_.END)
        self.tab_index = None
        self.tab_count = 0

    def show_completions(self, comps):
        self.text.delete(self.tab_index, Tk_.END)
        width = self.text.winfo_width()
        font = Font(self.text, self.text.cget('font'))
        charwidth = width//self.char_size
        biggest = 2 + max([len(x) for x in comps])
        num_cols = charwidth//biggest
        num_rows = (len(comps) + num_cols -1)//num_cols
        rows = []
        format = '%%-%ds'%biggest
        for n in range(0, num_rows):
            rows.append(''.join([format%x for x in comps[n:len(comps):num_rows]]))
        view = '\n'.join(rows)
        self.text.insert(self.tab_index, '\n'+view)
        self.text.mark_set(Tk_.INSERT, self.tab_index)
        self.window.update_idletasks()
        self.text.see(Tk_.END)

    def clear_completions(self):
        if self.tab_index:
            self.text.delete(self.tab_index, Tk_.END)
            self.tab_index = None
            self.tab_count = 0
 
    def stem(self, wordlist):
        if len(wordlist) == 1:
            return wordlist[0]
        for n in range(1,100):
            heads = set([x[:n] for x in wordlist])
            if len(heads) == 0:
                return ''
            elif len(heads) == 1:
                result = heads.pop()
            else:
                return result

    def write_history(self):
        self.text.see('output_end')
        self.window.update_idletasks()
        input = self.filtered_hist[-self.hist_pointer]
        input = input.replace('\n\n', '\n').strip('\n')
        if input.find('\n') > -1:
            margin = self.text.bbox('output_end')[0]
            margin -= Font(self.text).measure(' ')
            self.editing_hist = True
            self.text.tag_config('history',
                                 lmargin1=margin,
                                 lmargin2=margin,
                                 background='White')
            self.text.tag_lower('history')
            self.write(input +'\n', style=('history',), mutable=True)
            self.text.mark_set('history_end', Tk_.INSERT)
            self.text.mark_set(Tk_.INSERT, 'output_end')
        else:
            self.write(input, style=(), mutable=True)
        self.text.see(Tk_.INSERT)
            
    def handle_up(self, event):
        if self.text.compare(Tk_.INSERT, '<', 'output_end'):
            return
        insert_line = str(self.text.index(Tk_.INSERT)).split('.')[0] 
        prompt_line = str(self.text.index('output_end')).split('.')[0]
        if insert_line != prompt_line:
            return
        if self.hist_pointer == 0:
            input_history = self.IP.history_manager.input_hist_raw
            self.hist_stem = self.text.get('output_end', Tk_.END).strip()
            self.filtered_hist = [x for x in input_history
                                  if x.startswith(self.hist_stem)]
        if self.hist_pointer >= len(self.filtered_hist):
            self.window.bell()
            return 'break'
        self.text.delete('output_end', Tk_.END)
        self.hist_pointer += 1
        self.write_history()
        return 'break'

    def handle_down(self, event):
        if self.text.compare(Tk_.INSERT, '<', 'output_end'):
            return
        if self.editing_hist:
            insert_line = int(str(self.text.index(Tk_.INSERT)).split('.')[0]) 
            bottom_line = int(str(self.text.index('history_end')).split('.')[0])
            if insert_line < bottom_line - 1:
                return
        if self.hist_pointer == 0:
            self.window.bell()
            return 'break'
        self.text.delete('output_end', Tk_.END)
        self.hist_pointer -= 1
        if self.hist_pointer == 0:
            self.write(self.hist_stem.strip('\n'),
                       style=(), mutable=True)
        else:
            self.write_history()
        return 'break'

    def paste(self, event):
        """
        Prevent messing around with immutable text.
        """
        clip = primary = ''
        try:
            clip = event.widget.selection_get(selection="CLIPBOARD")
        except:
            pass
        try: 
            primary = event.widget.selection_get(selection="PRIMARY")
        except:
            pass
        paste = primary if primary else clip
        if self.text.compare(Tk_.INSERT, '<', 'output_end'):
            self.text.mark_set(Tk_.INSERT, 'output_end')
        self.text.insert(Tk_.INSERT, paste)
        self.text.see(Tk_.INSERT)
        try:
            self.text.tag_remove(Tk_.SEL, Tk_.SEL_FIRST, Tk_.SEL_LAST)
        except Tk_.TclError:
            pass
        return 'break'

    def protect_text(self, event):
        try:
            if self.text.compare(Tk_.SEL_FIRST, '<', 'output_end'):
                self.window.bell()
                return 'break'
        except:
            pass

    def middle_mouse_down(self, event):
        # Part 1 of a nasty hack to prevent pasting into the immutable text.
        # Needed because returning 'break' does not prevent the paste.
        if self.text.compare(Tk_.CURRENT, '<', 'output_end'):
            self.window.bell()
            try:
                self.nasty = str(self.text.index(Tk_.CURRENT))
            except AttributeError:
                return 'break'
            paste = event.widget.selection_get(selection="PRIMARY")
            self.nasty_text = paste.split()[0]
            return 'break'

    def middle_mouse_up(self, event):
        # Part 2 of the nasty hack.
        if self.nasty:
            # The CURRENT mark may be off by 1 from the actual paste index
            # This will probably fail sometimes.
            start = self.text.search(self.nasty_text, index=self.nasty+'-2c')
            if start:
                self.text.delete(start, Tk_.INSERT)
            self.nasty = None
            self.nasty_text = None
        try:
            self.text.tag_remove(Tk_.SEL, Tk_.SEL_FIRST, Tk_.SEL_LAST)
        except Tk_.TclError:
            pass
        return 'break'

    def start_interaction(self):
        """
        Print the banner and issue the first prompt.
        """
        self.text.image_create(Tk_.END, image=self.icon)
        banner_label = Tk_.Label(self.text, text=self.banner,
                                 background='#ec0fffec0',
                                 foreground='DarkGreen',
                                 anchor=Tk_.W,
                                 justify=Tk_.LEFT,
                                 font=self.prefs['font'])
        self.text.window_create(Tk_.END, window=banner_label)
        self.text.mark_set('output_end', '2.0')
         # Set a reasonable default directory for files to be saved to.
        try:
            home = os.environ['HOME']
        except KeyError:
            home = os.path.expanduser("~")
        desktop = os.path.join(home, "Desktop")
        default_save_dir = desktop if os.path.exists(desktop) else home
        self.IP.magics_manager.magics['line']['cd']("-q " + default_save_dir)
        # Create the prompt and go!
        self.interact_prompt()
        self.text.mark_set('output_end',Tk_.INSERT)
 
    def interact_prompt(self):
        """
        Print an input prompt, or a continuation prompt.
        """
        if self.IP.more:
            try:
                prompt = self.IP.prompt_manager.render('in2')
            except:
                self.IP.showtraceback()
        else:
            try:
                prompt = self.IP.separate_in + self.IP.prompt_manager.render('in')
            except:
                self.IP.showtraceback()
        self.write(prompt)
        if self.scroll_back:
            self.window.after_idle(self.do_scroll_back, self.prompt_index)
            self.scroll_back = False
        self.prompt_index = self.text.index(Tk_.END)
            
    def interact_handle_input(self, line):
        self.IP.input_splitter.push(line)
        self.IP.more = self.IP.input_splitter.push_accepts_more()
        if not self.IP.more:
            source_raw = self.IP.input_splitter.source_raw_reset()[1]
            self.IP.run_cell(source_raw, store_history=True)
        
    def send_line(self, line):
        """
        Send one line of input to the interpreter, which will write
        the result on our Text widget.  Then issue a new prompt.
        """
        self.write('\n')
        #line = line.decode(self.IP.stdin_encoding)
        if line[0] == '\n':
            line = line[1:]
        line = line.replace('\n\n','\n')
        self.running_code = True
        handle_input = self.interact_handle_input
        interact_prompt = self.interact_prompt
        try:
            handle_input(line)
        except KeyboardInterrupt:
            self.write('(IP) Keyboard Interrupt: ')
            self.IP.resetbuffer()
        self.running_code = False
        interact_prompt()
        if self.editing_hist and not self.IP.more:
            self.text.tag_delete('history')
            self.editing_hist = False
        self.text.see(Tk_.INSERT)
        self.text.mark_set('output_end',Tk_.INSERT)
        if self.IP.more:
            indent_current_str = self.IP._indent_current_str
            self.text.insert(Tk_.INSERT, indent_current_str(), ())
        self.text.delete(Tk_.INSERT, Tk_.END)
        self.hist_pointer = 0
        self.hist_stem = ''
                   
    def write(self, string, style=('output',), mutable=False):
        """
        Writes a string containing ansi color escape sequences to our
        Text widget, starting at the output_end mark.
        """
        if self.quiet:
            return
        if self.interrupted:
            self.interrupted = False
            raise KeyboardInterrupt('Writing')
        self.text.mark_set(Tk_.INSERT, 'output_end')
        pairs = ansi_seqs.findall(string)
        for pair in pairs:
            code, text = pair
            tags = (code,) + style if code else style
            if text:
                self.text.insert(Tk_.INSERT, text, tags)
                self.output_count += len(text)
        if mutable is False:
            self.text.mark_set('output_end', Tk_.INSERT)
        self.text.see(Tk_.INSERT)
        # Give the Text widget a chance to update itself every
        # so often (but let's not overdo it!)
        if self.output_count > 2000:
            self.output_count = 0
            self.text.update()
            if self.interrupted:
                self.interrupted = False
                self.interrupt()

    def write2(self, string):
        """
        Write method for messages.  These go in the "immutable"
        part, so as not to confuse the prompt.
        """
        self.window.tkraise()
        self.text.mark_set('save_insert', Tk_.INSERT)
        self.text.mark_set('save_end', 'output_end')
        self.text.mark_set(Tk_.INSERT, str(self.text.index('output_end'))+'linestart')
        self.text.insert(Tk_.INSERT, string, ('output', 'msg',))
        self.text.mark_set(Tk_.INSERT, 'save_insert')
        self.text.see('output_end')
        self.text.update_idletasks()

    def writelines(self, lines):
        lines = iter(lines)
        for line in lines:
            self.write(line)

    def page(self, text):
        """
        Our pager.  Just writes the text to the screen then jumps back to
        the beginning.  You can scroll down to read it.
        """
        self.scroll_back = True
        self.write(str(text))

    def do_scroll_back(self, index):
        self.window.after_idle(self.text.see, index)
        self.text.mark_set(Tk_.INSERT, index)

    def page_down(self):
        insert_line = int(str(self.text.index(Tk_.INSERT)).split('.')[0]) 
        prompt_line = int(str(self.text.index('output_end')).split('.')[0])
        height = prompt_line - insert_line
        if height < 1:
            return
        scroll_amount = min(int(self.text.cget('height')) - 1, height) 
        self.text.yview_scroll(scroll_amount, Tk_.UNITS)
        if scroll_amount == height:
            self.text.mark_set(Tk_.INSERT, 'output_end')
        else:
            self.text.mark_set(Tk_.INSERT, '%s.0'%(insert_line + scroll_amount))             
    def flush(self):
        """
        Required for a stdout / stderr proxy.
        """
        self.text.update_idletasks()

class ListedInstance(object):
    def __init__(self):
        self.focus_var = Tk_.IntVar()

    def to_front(self):
        self.window.deiconify()
        self.window.lift()
        self.window_master.update_window_list()

    def focus(self, event):
        self.focus_var.set(1)
        return 'break'
        pass

    def unfocus(self, event):
        self.focus_var.set(0)
        return 'break'

class SnapPyTerm(TkTerm, ListedInstance):

    def __init__(self, the_shell):
        self.window_master = self
        self.window_list=[]
        self.menu_title = 'SnapPy Shell'
        TkTerm.__init__(self, the_shell, name='SnapPy Command Shell')
        self.prefs = SnapPyPreferences(self)
        self.start_interaction()
        self.edit_config(None)
        # Under OS X, the window shouldn't be closable:
        if sys.platform == 'darwin':
            assert str(self.window) == "."
            self.window.createcommand("::tk::mac::OpenDocument",
                                  self.OSX_open_filelist)
            self.window.eval("::tk::unsupported::MacWindowStyle style .  document {verticalZoom horizontalZoom collapseBox resizable}")
            really_disable_menu_items(self.menubar)
        else:
            self.window.tk.call('namespace', 'import', '::tk::dialog::file::')
            self.window.tk.call('set', '::tk::dialog::file::showHiddenBtn',  '1')
            self.window.tk.call('set', '::tk::dialog::file::showHiddenVar',  '0')


    def add_bindings(self):
        self.text.bind_all('<ButtonRelease-1>', self.edit_config)
        self.window.bind('<FocusIn>', self.focus)
        self.window.bind('<FocusOut>', self.unfocus)
        self.focus_var = Tk_.IntVar(value=1)

    def build_menus(self):
        self.menubar = menubar = Tk_.Menu(self.window)
        Python_menu = Tk_.Menu(menubar, name="apple")
        Python_menu.add_command(label='About SnapPy...', command=self.about)
        Python_menu.add_separator()
        Python_menu.add_command(label='SnapPy Preferences...',
                                command=self.edit_prefs)
        if sys.platform in ('linux2', 'win32'):
            Python_menu.add_separator()
            Python_menu.add_command(label='Quit SnapPy', command=self.close)
        menubar.add_cascade(label='SnapPy', menu=Python_menu)
        File_menu = Tk_.Menu(menubar, name='file')
        add_menu(self.window, File_menu, 'Open...', self.open_file)
        add_menu(self.window, File_menu, 'Open link...', self.open_link_file)
        add_menu(self.window, File_menu, 'Save', self.save_file, state='disabled')
        add_menu(self.window, File_menu, 'Save as...', self.save_file_as)
        menubar.add_cascade(label='File', menu=File_menu)
        Edit_menu = Tk_.Menu(menubar, name='edit')
        add_menu(self.window, Edit_menu, 'Cut', 
                 lambda event=None: self.text.event_generate('<<Cut>>'))
        add_menu(self.window, Edit_menu, 'Copy',
                 lambda event=None: self.text.event_generate('<<Copy>>'))
        add_menu(self.window, Edit_menu, 'Paste',
                 lambda event=None: self.text.event_generate('<<Paste>>'))
        add_menu(self.window, Edit_menu, 'Delete',
                 lambda event=None: self.text.event_generate('<<Clear>>')) 
        menubar.add_cascade(label='Edit', menu=Edit_menu)
        self.window_menu = Window_menu = Tk_.Menu(menubar, name='window')
        self.update_window_list()
        menubar.add_cascade(label='Window', menu=Window_menu)
        Help_menu = Tk_.Menu(menubar, name="help")
        Help_menu.add_command(label='Help on SnapPy...', command=self.howto)
        menubar.add_cascade(label='Help', menu=Help_menu)

    def update_window_list(self):
        self.window_menu.delete(0,'end')
        for instance in [self] + self.window_list:
            self.window_menu.add_command(
                label=instance.menu_title,
                command=instance.to_front)

    def add_listed_instance(self, instance):
        self.window_list.append(instance)

    def delete_listed_instance(self, instance):
        self.window_list.remove(instance)

    def edit_prefs(self):
        apple_menu = self.menubar.children['apple']
        #entry = 2 if sys.platform == 'darwin' else 3
        apple_menu.entryconfig(2, state='disabled')
        dialog = PreferenceDialog(self.window, self.prefs)
        if dialog.okay:
            answer = askyesno('Save?',
                              'Do you want to save these settings?')
            if answer:
                self.prefs.write_prefs()
        apple_menu.entryconfig(2, state='active')

    def edit_config(self, event):
        edit_menu = self.menubar.children['edit']
        if sys.platform == 'darwin':
            activate, deactivate, rest = [2], [], [0, 1, 3]
        else:
            activate, deactivate, rest = [], [0, 3], [1, 2]
        try:
            self.text.get(Tk_.SEL_FIRST, Tk_.SEL_LAST)
            activate += rest
        except Tk_.TclError:
            deactivate += rest

        for i in activate:
            edit_menu.entryconfig(i, state='active')
        for i in deactivate:
            edit_menu.entryconfig(i, state='disabled')

    def OSX_open_filelist(self, *args):
        for arg in args:
            sys.stderr.write(repr(arg)+'\n')

    def open_file(self, event=None):
        openfile = filedialog.askopenfile(
            title='Run Saved Transcript In Current Namespace',
            defaultextension='.py',
            filetypes = [
                ("Python and text files", "*.py *.ipy *.txt", "TEXT"),
                ("All text files", "", "TEXT"),
                ("All files", "")])
        if openfile:
            lines = openfile.readlines()
            openfile.close()
            if re.search("%\s*[lL]ink\s*[Pp]rojection", lines[0]):
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
            title='Load Link Projection File',
            defaultextension='.lnk',
            filetypes = [
                ("Link and text files", "*.lnk *.txt", "TEXT"),
                ("All text files", "", "TEXT"),
                ("All files", "")])
        if openfile:
            if not re.search("%\s*[lL]ink\s*[Pp]rojection", openfile.readline()):
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

    def about(self):
        InfoDialog(self.window, 'About SnapPy', about_snappy)

    def howto(self):
        doc_file = os.path.join(os.path.dirname(snappy.__file__),
                                'doc', 'index.html')
        doc_path = os.path.abspath(doc_file)
        url = 'file:' + pathname2url(doc_path)
        try:
            webbrowser.open(url) 
        except:
            tkMessageBox.showwarning('Not found!', 'Could not open URL\n(%s)'%url)

# These classes assume that the global variable "terminal" exists

class SnapPyBrowser(Browser, ListedInstance):
    def __init__(self, manifold):
        Browser.__init__(self, manifold, terminal.window)
        self.menu_title = self.window.title()
        self.focus_var = Tk_.IntVar(self.window)
        self.window_master = terminal
        self.window_master.add_listed_instance(self)
        self.window_master.update_window_list()
        self.window.bind('<FocusIn>', self.focus)
        self.window.bind('<FocusOut>', self.unfocus)
        if sys.platform=='darwin':
            really_disable_menu_items(self.menubar)

    def close(self, event=None):
        window_list = self.window_master.window_list
        if self in window_list:
                window_list.remove(self)
        self.window_master.update_window_list()
        self.window.destroy()

class SnapPyLinkEditor(LinkEditor, ListedInstance):
    def __init__(self, root=None, no_arcs=False, callback=None, cb_menu='',
                 title='PLink Editor'):
        self.focus_var = Tk_.IntVar()
        self.window_master = terminal
        LinkEditor.__init__(self, root=terminal.window, no_arcs=no_arcs,
                            callback=callback, cb_menu=cb_menu,
                            title=title)
        self.menu_title = self.window.title()
        self.window_master.add_listed_instance(self)
        self.window_master.update_window_list()
        self.window.bind('<FocusIn>', self.focus)
        self.window.bind('<FocusOut>', self.unfocus)
        self.window.focus_set()
        self.window.update_idletasks()
        if sys.platform == 'darwin':
            really_disable_menu_items(self.menubar)



    def focus(self, event):
        self.focus_in(event)
        ListedInstance.focus(self, event)

    def unfocus(self, event):
        self.focus_out(event)
        ListedInstance.unfocus(self, event)

    def build_menus(self):
        self.menubar = menubar = Tk_.Menu(self.window)
        Python_menu = Tk_.Menu(menubar, name="apple")
        Python_menu.add_command(label='About PLink...', command=self.about)
        Python_menu.add_separator()
        Python_menu.add_command(label='Preferences...', state='disabled')
        Python_menu.add_separator()
        if sys.platform == 'linux2':
            Python_menu.add_command(label='Quit SnapPy', command=terminal.close)
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

        Edit_menu = Tk_.Menu(menubar, name='edit')

        add_menu(self.window, Edit_menu, 'Cut', None, state='disabled')
        add_menu(self.window, Edit_menu, 'Copy', None, state='disabled')
        add_menu(self.window, Edit_menu, 'Paste', None, state='disabled')
        add_menu(self.window, Edit_menu, 'Delete', None, state='disabled')
        menubar.add_cascade(label='Edit', menu=Edit_menu)
        self.build_plink_menus() # Application Specific Menus
        Tools_menu = self.tools_menu
        if self.callback:
            Tools_menu.add_separator()
            Tools_menu.add_command(label=self.cb_menu, command=self.do_callback)
        Window_menu = self.window_master.menubar.children['window']
        menubar.add_cascade(label='Window', menu=Window_menu)
        Help_menu = Tk_.Menu(menubar, name="help")
        Help_menu.add_command(label='Help on PLink ...', command=self.howto)
        menubar.add_cascade(label='Help', menu=Help_menu)
        self.window.config(menu=menubar)

    def to_front(self):
        self.reopen()
        self.window.tkraise()
        self.focus_var.set(1)
        self.window_master.update_window_list()

    def copy_info(self):
        if not self.infotext.selection_present():
           self.infotext.selection_range(0, Tk_.END)
        self.infotext.focus()
        self.infotext.event_generate('<<Copy>>')

    def load(self, event=None, file_name=None):
        LinkEditor.load(self, file_name)

    def save(self, event=None):
        LinkEditor.save(self)

class SnapPyPolyhedronViewer(PolyhedronViewer, ListedInstance):
    def __init__(self, facedicts, root=None, title='Polyhedron Viewer'):
        self.focus_var = Tk_.IntVar()
        self.window_master = terminal
        PolyhedronViewer.__init__(self, facedicts, root=terminal.window,
                                  title=title)
        self.menu_title = self.window.title()
        self.window_master.add_listed_instance(self)
        self.window_master.update_window_list()
        self.window.bind('<FocusIn>', self.focus)
        self.window.bind('<FocusOut>', self.unfocus)
        if sys.platform=='darwin':
            really_disable_menu_items(self.menubar)

    def add_help(self):
        pass

    build_menus = dirichlet_menus

    def close(self):
        self.polyhedron.destroy()
        self.window_master.window_list.remove(self)
        self.window_master.update_window_list()
        self.window.destroy()

    def save_image(self):
        togl_save_image(self)

class SnapPyHoroballViewer(HoroballViewer, ListedInstance):
    def __init__(self, nbhd, which_cusp=0, cutoff=None,
                 root=None, title='Horoball Viewer'):
        self.focus_var = Tk_.IntVar()
        self.window_master = terminal
        HoroballViewer.__init__(self, nbhd, which_cusp=which_cusp,
                                cutoff=cutoff, root=terminal.window,
                                title=title, prefs = terminal.prefs)
        self.menu_title = self.window.title()
        self.window_master.add_listed_instance(self)
        self.window_master.update_window_list()
        self.window.bind('<FocusIn>', self.focus)
        self.window.bind('<FocusOut>', self.unfocus)
        self.view_check()
        if sys.platform=='darwin':
            really_disable_menu_items(self.menubar)

    build_menus = horoball_menus

    def view_check(self):
        if self.horo_var.get():
            self.widget.set_background(0.3, 0.3, 0.4)
        else:
            self.widget.set_background(1.0, 1.0, 1.0)
        self.widget.tkRedraw()
        
    def close(self):
        self.widget.activate()
        self.scene.destroy()
        self.window_master.window_list.remove(self)
        self.window_master.update_window_list()
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
def IPython_pager(self, text):
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
    SnapPy_ns['Browser'] = SnapPyBrowser
    helper = pydoc.Helper(input=terminal, output=terminal)
    helper.__call__ = lambda x=None : helper.help(x) if x else terminal.howto()
    helper.__repr__ = lambda : help_banner
    SnapPy_ns['help'] = helper
    io.stdout = io.stderr = sys.stdout = sys.stderr = terminal
    SnapPy_ns['io'] = io
    the_shell.user_ns.update(SnapPy_ns)
    snappy.browser.window_master = terminal
    snappy.SnapPy.LinkEditor = SnapPyLinkEditor
    snappy.SnapPy.PolyhedronViewer = SnapPyPolyhedronViewer
    snappy.SnapPy.HoroballViewer = SnapPyHoroballViewer
    snappy.SnapPy.Browser = SnapPyBrowser
    snappy.SnapPy.msg_stream.write = terminal.write2
    snappy.SnapPy.UI_callback = terminal.SnapPea_callback
    if not snappy.SnapPy._within_sage:
        snappy.pari.UI_callback = terminal.PARI_callback
    terminal.window.lift()
    terminal.window.mainloop()

if __name__ == "__main__":
    main()
