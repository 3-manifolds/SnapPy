import __builtin__
import IPython
import Tkinter as Tk_
import tkFileDialog
import tkMessageBox
from tkFont import Font
from tkMessageBox import askyesno
import os, sys, re, webbrowser, signal
from plink import LinkEditor
from urllib import pathname2url
from pydoc import help
import time
import snappy
import snappy.version
from snappy import SnapPeaFatalError
from snappy.polyviewer import PolyhedronViewer
from snappy.horoviewer import HoroballViewer
from snappy.SnapPy import SnapPea_interrupt
from snappy.SnapPy import msg_stream
from snappy.shell import SnapPy_showtraceback
from snappy.preferences import Preferences, PreferenceDialog

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

OSX_shortcuts = {'Open'   : u'\t\t\u2318O',
                 'Save'   : u'\t\t\u2318S',
                 'SaveAs' : u'\t\u2318\u21e7S',
                 'Cut'    : u'\t\t\u2318X',
                 'Copy'   : u'\t\u2318C',
                 'Paste'  : u'\t\u2318V'}

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
    
class Tk(Tk_.Tk):
    def __init__(self, error_handler):
        Tk_.Tk.__init__(self)
        # Tkinter ignores exceptions raised by callbacks, but
        # calls this function to report their occurence.
        self.report_callback_exception = error_handler

class TkTerm:
    """
    A Tkinter terminal window that runs an IPython shell.
    Some ideas borrowed from code written by Eitan Isaacson, IBM Corp.
    """
    def __init__(self, the_shell, name='TkTerm'):
        self.shell = the_shell
        if debug_Tk:
            self.window = window = Tk_.Tk()
        else:
            self.window = window = Tk(self.report_callback_exception)
        try:
            window.tk.call('console', 'hide')
        except Tk_.TclError:
            pass
        # Supposedly the following is ignored by Windows and Mac
        window.option_add('*Dialog.msg.wrapLength', '20i')
        window.title(name)
        window.protocol("WM_DELETE_WINDOW", self.close)
        self.frame = frame = Tk_.Frame(window)
        self.text = text = Tk_.Text(frame,
                                    foreground='Black',
                                    borderwidth=3,
                                    background='#ec0fffec0',
                                    highlightthickness=0,
                                    relief=Tk_.FLAT
                                )
        self.scroller = scroller = Tk_.Scrollbar(frame, command=text.yview)
        text.config(yscrollcommand = scroller.set)
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
        text.bind_all('<ButtonPress-2>', self.middle_mouse_down)
        text.bind_all('<ButtonRelease-2>', self.middle_mouse_up)
        text.bind('<Button-3>', lambda event : 'break')
        text.bind('<Button-4>', lambda event : 'break')
        text.bind('<MouseWheel>', lambda event : text.yview_scroll(-1, Tk_.UNITS))
        if sys.platform == 'darwin':
            self.window.bind_all('<Command-Key-q>', self.close_event)
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
        self.output_count = 0
        self.IP = the_shell.IP
        self.IP.magic_colors('LightBG')
        if not debug_Tk:
            self.IP.write = self.write                 # used for the prompt
            IPython.Shell.Term.cout.write = self.write # used for output
            IPython.Shell.Term.cerr.write = self.write # used for tracebacks
            sys.stdout = self # also used for tracebacks (why???)
        sys.displayhook = self.IP.outputcache
        if the_shell.banner:
            self.banner = the_shell.banner
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
            signal.signal(signal.SIGALRM, self.UI_ticker)
            try:
                signal.setitimer(signal.ITIMER_REAL, 1.0, 1.0)
            except AttributeError: # itimer is not supported in python 2.5
                pass
        self.start_interaction()

    # For subclasses to override:
    def build_menus(self):
        pass

    # For subclasses to override:
    def add_bindings(self):
        pass

    def UI_ticker(self, signal, stackframe):
        self.interrupted = False
        self.window.update()
        if self.interrupted:
            if self.running_code:
                self.interrupted = False
                if self.aborted_SnapPea:
                    self.aborted_SnapPea = False
                    raise KeyboardInterrupt, 'SnapPea computation aborted'
                else:
                    raise KeyboardInterrupt, 'Running'

    def interrupt(self):
        # Tell the ticker to raise a KeyboardInterrupt after the update
        if not self.running_code:
            raise KeyboardInterrupt, 'Halted'
        else:
            self.interrupted = True
            # Inform the SnapPea kernel about the interrupt.
            self.aborted_SnapPea = SnapPea_interrupt()
            
    def report_callback_exception(self, exc, value, traceback):
        # This is called when exceptions are caught by Tk.
        self.write2('(Tk) ' + exc.__name__ +': ' + str(value) +'\n')
        sys.last_type = exc
        sys.last_value = value
        sys.last_traceback = traceback
#        self.IP.IPtraceback((exc, value, traceback))
    
    def set_font(self, fontdesc):
        self.text.config(font=fontdesc)
        self.char_size = Font(self.text, fontdesc).measure('M')
        self.text.tag_add('all', '1.0', Tk_.END)
        self.text.tag_config('all', font=fontdesc)

    def close(self):
        self.window.update_idletasks()
        self.window.quit()

    def close_event(self, event):
        self.close()

    def handle_keypress(self, event):
        self.clear_completions()
        if event.char == '\001':
            self.text.mark_set(Tk_.INSERT, 'output_end')
            return 'break'
        if event.char == '\025':
            self.text.delete('output_end', Tk_.END)
            return 'break'
        if event.char == '\003':
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
        self.tab_index = self.text.index(Tk_.INSERT)
        self.tab_count += 1
        if self.tab_count > 2:
            self.clear_completions()
            return 'break'
        line = self.text.get('output_end', self.tab_index).strip('\n')
        word = delims.split(line)[-1]
        if word[-1] != '_':
            completions = [x for x in self.IP.complete(word) 
                           if x.find('__') == -1 and x.find('._') == -1]
        else:
            completions = self.IP.complete(word)
        if len(completions) == 0:
            self.window.bell()
            return 'break'
        stem = self.stem(completions)
        if len(stem) > len(word):
            self.do_completion(word, stem)
        elif len(completions) > 60 and self.tab_count == 1:
            self.show_completions(['%s possibilities -- hit tab again to view them all'%len(completions)])
        else:
            self.show_completions(completions)
        return 'break'

    def do_completion(self, word, completion):
        tail = completion[len(word):]
        self.text.insert(self.tab_index, tail)
        self.text.delete(Tk_.INSERT, Tk_.END)
        self.tab_index = None

    def show_completions(self, comps):
        width = self.text.winfo_width()
        font = Font(self.text, self.text.cget('font'))
        charwidth = width/self.char_size
        biggest = 2 + max([len(x) for x in comps])
        num_cols = charwidth/biggest
        num_rows = (len(comps) + num_cols -1)/num_cols
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

    def history_check(self):
        if self.text.compare(Tk_.INSERT, '<', 'output_end'):
            return False
        insert_line = str(self.text.index(Tk_.INSERT)).split('.')[0] 
        prompt_line = str(self.text.index('output_end')).split('.')[0]
        if insert_line > prompt_line:
            return False
        self.history_index = self.text.index(Tk_.INSERT)
        return True

    def write_history(self):
        self.text.see('output_end')
        self.window.update_idletasks()
        input = self.filtered_hist[-self.hist_pointer]
        input = input.replace('\n\n', '\n').strip('\n')
        if input.find('\n') > -1:
            input = '\n'+input
            margin = self.text.bbox('output_end')[0]
            margin -= Font(self.text).measure(' ')
            self.editing_hist = True
            self.text.tag_config('history',
                                 lmargin1=margin,
                                 lmargin2=margin,
                                 background='White')
            self.text.tag_lower('history')
            self.write(input + '\n', style=('history',), mutable=True)
        else:
            self.write(input, style=(), mutable=True)
        self.text.mark_set(Tk_.INSERT, self.history_index)
        self.text.see(Tk_.INSERT)
            
    def handle_up(self, event):
        if self.history_check() is False:
            return
        insert_line = str(self.text.index(Tk_.INSERT)).split('.')[0] 
        prompt_line = str(self.text.index('output_end')).split('.')[0]
        if insert_line > prompt_line:
            return
        if self.hist_pointer == 0:
            self.hist_stem = self.text.get('output_end', Tk_.END).strip('\n')
            self.filtered_hist = [x for x in self.IP.input_hist_raw
                                  if x.startswith(self.hist_stem)]
        if self.hist_pointer >= len(self.filtered_hist):
            self.window.bell()
            return 'break'
        self.text.delete('output_end', Tk_.END)
        self.hist_pointer += 1
        self.write_history()
        return 'break'

    def handle_down(self, event):
        if self.history_check() is False:
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
        if self.text.compare(Tk_.SEL_FIRST, '<', 'output_end'):
            self.window.bell()
            return 'break'

    def middle_mouse_down(self, event):
        # Part 1 of a nasty hack to prevent pasting into the immutable text.
        # Needed because returning 'break' does not prevent the paste.
        if self.text.compare(Tk_.CURRENT, '<', 'output_end'):
            self.window.bell()
            self.nasty = str(self.text.index(Tk_.CURRENT))
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
        self.text.tag_config('banner', foreground='DarkGreen')
        self.write(self.banner, style=('output', 'banner'))
        # Set a reasonable default directory for files to be saved to.
        home = os.path.expanduser("~")
        desktop = home + os.sep + "Desktop"
        default_save_dir = desktop if os.path.exists(desktop) else home
        self.IP.magic_cd("-q " + default_save_dir)
        # Create the prompt and go!
        self.IP.interact_prompt()
        self.text.mark_set('output_end',Tk_.INSERT)
 
    def send_line(self, line):
        """
        Send one line of input to the interpreter, who will write
        the result on our Text widget.  Then issue a new prompt.
        """
        self.write('\n')
        #line = line.decode(self.IP.stdin_encoding)
        if line[0] == '\n':
            line = line[1:]
        line = line.replace('\n\n','\n')
        self.running_code = True
        try:
            self.IP.interact_handle_input(line)
        except KeyboardInterrupt:
            self.write('(IP) Keyboard Interrupt: ')
            self.IP.resetbuffer()
        self.running_code = False
        self.IP.interact_prompt()
        if self.editing_hist and not self.IP.more:
            self.text.tag_delete('history')
            self.editing_hist = False
        self.text.see(Tk_.INSERT)
        self.text.mark_set('output_end',Tk_.INSERT)
        if self.IP.more:
            self.text.insert(Tk_.INSERT, self.IP.indent_current_str(), ())
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
            raise KeyboardInterrupt, 'Writing'
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

    def flush(self):
        """
        Since we are pretending to be an IOTerm.
        """
        self.text.update_idletasks()

class ListedInstance(object):
    def __init__(self):
        self.focus_var = Tk_.IntVar()

    def to_front(self):
        self.window.tkraise()
        self.focus_var.set(1)
        self.window_master.update_window_list()

    def focus(self, event):
        self.focus_var.set(1)
        return 'break'

    def unfocus(self, event):
        self.focus_var.set(0)
        return 'break'

class SnapPyTerm(TkTerm, ListedInstance):

    def __init__(self, the_shell):
        self.window_master = self
        self.window_list=[]
        self.title='SnapPy Shell'
        TkTerm.__init__(self, the_shell, name='SnapPy Command Shell')
        self.prefs = SnapPyPreferences(self)
        self.edit_config(None)
        # Under OS X, the window shouldn't be closeable:
        if sys.platform == 'darwin':
            assert str(self.window) == "."
            self.window.createcommand("::tk::mac::OpenDocument",
                                  self.OSX_open_filelist)
            self.window.eval("::tk::unsupported::MacWindowStyle style .  document {verticalZoom horizontalZoom collapseBox resizable}")
            
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
        File_menu.add_command(
            label='Open...' + scut['Open'],
            command=self.open_file)
        File_menu.add_command(
            label='Open link...',
            command=self.open_link_file)
        File_menu.add_command(
            label='Save' + scut['Save'],
            command=self.save_file, state='disabled')
        File_menu.add_command(
            label='Save as...' + scut['SaveAs'],
            command=self.save_file_as)
        menubar.add_cascade(label='File', menu=File_menu)
        Edit_menu = Tk_.Menu(menubar, name='edit')
        Edit_menu.add_command(
            label='Cut' + scut['Cut'],
            command=lambda : self.text.event_generate('<<Cut>>')) 
        Edit_menu.add_command(
            label='Copy' + scut['Copy'],
            command=lambda : self.text.event_generate('<<Copy>>'))  
        Edit_menu.add_command(
            label='Paste' + scut['Paste'],
            command=lambda : self.text.event_generate('<<Paste>>'))
        Edit_menu.add_command(
            label='Delete',
            command=lambda : self.text.event_generate('<<Clear>>')) 
        menubar.add_cascade(label='Edit', menu=Edit_menu)
        self.window_menu = Window_menu = Tk_.Menu(menubar, name='window')
        self.update_window_list()
        menubar.add_cascade(label='Window', menu=Window_menu)
        Help_menu = Tk_.Menu(menubar, name="help")
        Help_menu.add_command(label='Help on SnapPy...', command=self.howto)
        menubar.add_cascade(label='Help', menu=Help_menu)
        self.window.config(menu=menubar)

    def update_window_list(self):
        self.window_menu.delete(0,'end')
        for instance in [self] + self.window_list:
            self.window_menu.add_checkbutton(
                label=instance.title,
                variable=instance.focus_var,
                command=instance.to_front)

    def add_listed_instance(self, instance):
        self.window_list.append(instance)

    def delete_listed_instance(self, instance):
        self.window_list.remove(instance)

    def edit_prefs(self):
        apple_menu = self.menubar.children['apple']
        entry = 2 if sys.platform == 'darwin' else 3
        apple_menu.entryconfig(entry, state='disabled')
        dialog = PreferenceDialog(self.window, self.prefs)
        if dialog.okay:
            answer = askyesno('Save?',
                              'Do you want to save these settings?')
            if answer:
                self.prefs.write_prefs()
        apple_menu.entryconfig(entry, state='active')

    def edit_config(self, event):
        edit_menu = self.menubar.children['edit']
        try:
            self.text.get(Tk_.SEL_FIRST, Tk_.SEL_LAST)
            for n in (0,1,2,3):
                edit_menu.entryconfig(n, state='active')
        except Tk_.TclError:
            for n in (0,1,3):
                edit_menu.entryconfig(n, state='disabled')

    def OSX_open_filelist(self, *args):
        for arg in args:
            print >> sys.stderr, arg

    def open_file(self):
        openfile = tkFileDialog.askopenfile(
            title='Run Saved Transcript In Current Namespace',
            defaultextension='.py',
            filetypes = [
                ("Python and text files", "*.py *.ipy *.txt", "TEXT"),
                ("All text files", "", "TEXT"),
                ("All files", "")])
        if openfile:
            lines = openfile.readlines()
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
                    self.IP.interact_handle_input(line)
                    self.IP.interact_prompt()
                    
        openfile.close()

    def open_link_file(self):
        openfile = tkFileDialog.askopenfile(
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
                self.IP.interact_handle_input(line)
                self.IP.interact_prompt()
                M = self.IP.output_hist[self.IP.outputcache.prompt_count - 1]
                M.LE.load(openfile.name)

    def save_file_as(self):
        savefile = tkFileDialog.asksaveasfile(
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
            inputs = self.IP.input_hist
            results = self.IP.output_hist
            for n in range(1,len(inputs)):
                savefile.write('\n'+re.sub('\n+','\n',inputs[n]))
                try:
                    output = repr(results[n]).split('\n')
                except KeyError:
                    continue
                for line in output:
                    savefile.write('# ' + line + '\n')
            savefile.close()

    def save_file(self):
        self.window.bell()
        self.write2('Save As\n')

    def about(self):
        tkMessageBox.showinfo('About SnapPy', u"""
SnapPy is a user interface for the SnapPea kernel,
which was written by Jeff Weeks.  SnapPy was
written by Marc Culler and Nathan Dunfield and
is distributed under the GNU Public License,
version 2.  Its home page is:
     http://www.math.uic.edu/~t3m

The release number of this version is %s.

SnapPy is written in the Python language, using
Cython to incorporate the SnapPea kernel code.
The graphical interface uses Tcl/Tk, via Python's
Tkinter module.

For information on how to use SnapPy, please see the
Help menu.

Information, downloads, and source code for the
SnapPea kernel and for the user interfaces written
by Jeff Weeks are available at:
     http://www.geometrygames.org/SnapPea-old/
     http://www.geometrygames.org/SnapPea/

Copyright \u00a9 2009, Marc Culler, Nathan Dunfield
and others.
"""%snappy.version.version)

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

class SnapPyLinkEditor(LinkEditor, ListedInstance):
    def __init__(self, root=None, no_arcs=False, callback=None, cb_menu='',
                 title='PLink Editor'):
        self.focus_var = Tk_.IntVar()
        self.window_master = terminal
        LinkEditor.__init__(self, terminal.window, no_arcs, callback,
                            cb_menu, title)
        self.window.bind('<FocusIn>', self.focus)
        self.window.bind('<FocusOut>', self.unfocus)

    def build_menus(self):
        self.menubar = menubar = Tk_.Menu(self.window)
        Python_menu = Tk_.Menu(menubar, name="apple")
        Python_menu.add_command(label='About PLink ...', command=self.about)
        Python_menu.add_separator()
        Python_menu.add_command(label='Preferences ...', state='disabled')
        Python_menu.add_separator()
        if sys.platform == 'linux2':
            Python_menu.add_command(label='Quit SnapPy', command=terminal.close)
        menubar.add_cascade(label='SnapPy', menu=Python_menu)
        File_menu = Tk_.Menu(menubar, name='file')
        File_menu.add_command(
            label='Open...' + scut['Open'],
            command=self.load)
        File_menu.add_command(
            label='Save as...' + scut['SaveAs'],
            command=self.save)
        Print_menu = Tk_.Menu(menubar, name='print')
        Print_menu.add_command(label='monochrome',
                               command=lambda : self.save_image(color_mode='mono'))
        Print_menu.add_command(label='color', command=self.save_image)
        File_menu.add_cascade(label='Save Image', menu=Print_menu)
        File_menu.add_separator()
        if self.callback:
            File_menu.add_command(label='Close', command=self.done)
        else:
            File_menu.add_command(label='Exit', command=self.done)
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
        # Application Specific Menus
        PLink_menu = Tk_.Menu(menubar)
        PLink_menu.add_command(label='Make alternating',
                       command=self.make_alternating)
        PLink_menu.add_command(label='Reflect', command=self.reflect)
        PLink_menu.add_command(label='Clear', command=self.clear)
        Info_menu = Tk_.Menu(PLink_menu)
        Info_menu.add_command(label='DT code', command=self.dt_normal)
        Info_menu.add_command(label='DT for Snap', command=self.dt_snap)
        Info_menu.add_command(label='Gauss code', command=self.not_done)
        Info_menu.add_command(label='PD code', command=self.not_done)
        PLink_menu.add_cascade(label='Info', menu=Info_menu)
        PLink_menu.add_separator()
        if self.callback:
            PLink_menu.add_command(label=self.cb_menu, command=self.do_callback)
        menubar.add_cascade(label='PLink', menu=PLink_menu)
        #
        Window_menu = self.window_master.menubar.children['window']
        self.window_master.add_listed_instance(self)
        self.window_master.update_window_list()
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

class SnapPyPolyhedronViewer(PolyhedronViewer, ListedInstance):
    def __init__(self, facedicts, root=None, title=u'Polyhedron Viewer'):
        self.focus_var = Tk_.IntVar()
        self.window_master = terminal
        PolyhedronViewer.__init__(self, facedicts, root=terminal.window,
                                  title=title)
        self.window.bind('<FocusIn>', self.focus)
        self.window.bind('<FocusOut>', self.unfocus)

    def add_help(self):
        pass

    def build_menus(self):
        self.menubar = menubar = Tk_.Menu(self.window)
        Python_menu = Tk_.Menu(menubar, name="apple")
        Python_menu.add_command(label='About SnapPy ...')
        Python_menu.add_separator()
        Python_menu.add_command(label='Preferences ...', state='disabled')
        Python_menu.add_separator()
        if sys.platform == 'linux2':
            Python_menu.add_command(label='Quit SnapPy', command=terminal.close)
        menubar.add_cascade(label='SnapPy', menu=Python_menu)
        File_menu = Tk_.Menu(menubar, name='file')
        File_menu.add_command(
            label='Open...' + scut['Open'], state='disabled')
        File_menu.add_command(
            label='Save as...' + scut['SaveAs'], state='disabled')
        Print_menu = Tk_.Menu(menubar, name='print')
        Print_menu.add_command(label='monochrome',
                               command=lambda : self.save_image(color_mode='mono'),
                               state='disabled')
        Print_menu.add_command(label='color',
                               command=lambda : self.save_image(color_mode='color'),
                               state='disabled')
        File_menu.add_cascade(label='Save Image', menu=Print_menu)
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
        Window_menu = self.window_master.menubar.children['window']
        self.window_master.add_listed_instance(self)
        self.window_master.update_window_list()
        menubar.add_cascade(label='Window', menu=Window_menu)
        Help_menu = Tk_.Menu(menubar, name="help")
        Help_menu.add_command(label='Help on PolyhedronViewer ...',
                              command=self.widget.help)
        menubar.add_cascade(label='Help', menu=Help_menu)
        self.window.config(menu=menubar)

    def close(self):
        self.window_master.window_list.remove(self)
        self.window_master.update_window_list()
        self.window.destroy()

class SnapPyHoroballViewer(HoroballViewer, ListedInstance):
    def __init__(self, nbhd, cutoff=0.1, which_cusp=0,
                 root=None, title='Horoball Viewer'):
        self.focus_var = Tk_.IntVar()
        self.window_master = terminal
        HoroballViewer.__init__(self, nbhd, cutoff, which_cusp,
                                root=terminal.window, title=title)
        self.window.bind('<FocusIn>', self.focus)
        self.window.bind('<FocusOut>', self.unfocus)

    def build_menus(self):
        self.menubar = menubar = Tk_.Menu(self.window)
        Python_menu = Tk_.Menu(menubar, name="apple")
        Python_menu.add_command(label='About SnapPy ...')
        Python_menu.add_separator()
        Python_menu.add_command(label='Preferences ...',  state='disabled')
        Python_menu.add_separator()
        if sys.platform == 'linux2':
            Python_menu.add_command(label='Quit SnapPy', command=terminal.close)
        menubar.add_cascade(label='SnapPy', menu=Python_menu)
        File_menu = Tk_.Menu(menubar, name='file')
        File_menu.add_command(
            label='Open...' + scut['Open'], state='disabled')
        File_menu.add_command(
            label='Save as...' + scut['SaveAs'], state='disabled')
        Print_menu = Tk_.Menu(menubar, name='print')
        Print_menu.add_command(label='monochrome',
                               command=lambda : self.save_image(color_mode='mono'),
                               state='disabled')
        Print_menu.add_command(label='color',
                               command=lambda : self.save_image(color_mode='color'),
                               state='disabled')
        File_menu.add_cascade(label='Save Image', menu=Print_menu)
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
                                  command=self.widget.tkRedraw,
                                  variable=self.pgram_var)
        View_menu.add_checkbutton(label='Ford edges',
                                  command=self.widget.tkRedraw,
                                  variable=self.Ford_var)
        View_menu.add_checkbutton(label='triangulation',
                                  command=self.widget.tkRedraw,
                                  variable=self.tri_var)
        menubar.add_cascade(label='View', menu=View_menu)
        Window_menu = self.window_master.menubar.children['window']
        self.window_master.add_listed_instance(self)
        self.window_master.update_window_list()
        menubar.add_cascade(label='Window', menu=Window_menu)
        Help_menu = Tk_.Menu(menubar, name="help")
        Help_menu.add_command(label='Help on HoroballViewer ...',
                              command=self.widget.help)
        menubar.add_cascade(label='Help', menu=Help_menu)
        self.window.config(menu=menubar)

    def close(self):
        self.window_master.window_list.remove(self)
        self.window_master.update_window_list()
        self.window.destroy()

class SnapPyPreferences(Preferences):
    def __init__(self, terminal):
        self.terminal = terminal
        Preferences.__init__(self, terminal.text)
        self.apply_prefs()

    def apply_prefs(self):
        self.terminal.set_font(self['font'])
        self.terminal.window.update_idletasks()
        changed = self.changed()
        IP = self.terminal.shell.IP
        self.terminal.quiet = True
        if 'autocall' in changed:
            if self.prefs_dict['autocall']:
                IP.magic_autocall(2)
            else:
                IP.magic_autocall(0)
        if 'automagic' in changed:
            if self.prefs_dict['automagic']:
                IP.magic_automagic('on')
            else:
                IP.magic_automagic('off')
        if 'tracebacks' in changed:
            IP.tracebacks = self.prefs_dict['tracebacks']
            if IP.tracebacks:
                self.terminal.write('\nTracebacks are enabled\n')
            else:
                self.terminal.write('\nTracebacks are disabled\n')
        self.terminal.quiet = False

app_banner = """
    Hi.  It's SnapPy.  
    SnapPy is based on the SnapPea kernel, written by Jeff Weeks.
    Type "Manifold?" to get started.
    """

class Quitter(object):
    """
    Replacement for the IPython Quitter
    """

    def __init__(self,shell,name):
        self.shell = shell
        self.name = name
        
    def __repr__(self):
        return 'Please use the SnapPy menu to quit.'
    __str__ = __repr__

    def __call__(self):
        return self

def main():
    argv = None  if not sys.platform == 'win32' else ["-ipythondir", os.path.join(os.environ['USERPROFILE'], ".ipython")]
    the_shell = IPython.Shell.IPShellEmbed(argv=argv, banner=app_banner)
    the_shell.IP.tracebacks = False
    the_shell.IP.IP_showtraceback = the_shell.IP.showtraceback
    the_shell.IP.showtraceback = SnapPy_showtraceback
    # To restore tracebacks: __IP.tracebacks = True
    __builtin__.exit = Quitter(the_shell,'exit')
    __builtin__.quit = Quitter(the_shell,'quit')
    the_shell.IP.magic_exit = lambda x: __builtin__.exit
    the_shell.IP.magic_quit = lambda x: __builtin__.quit
    
    SnapPy_ns = dict([(x, getattr(snappy,x)) for x in snappy.__all__])
    SnapPy_ns['help'] = help
    the_shell.IP.user_ns.update(SnapPy_ns)
    os.environ['TERM'] = 'dumb'
    global terminal
    terminal = SnapPyTerm(the_shell)
    the_shell.IP.tkterm = terminal
    snappy.SnapPy.LinkEditor = SnapPyLinkEditor
    snappy.SnapPy.PolyhedronViewer = SnapPyPolyhedronViewer
    snappy.SnapPy.HoroballViewer = SnapPyHoroballViewer
    snappy.SnapPy.msg_stream.write = terminal.write2
    terminal.window.mainloop()

if __name__ == "__main__":
    main()
