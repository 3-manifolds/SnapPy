# -*- coding: utf-8 -*-
import os, sys, re, signal, IPython
from builtins import range
from IPython.utils import io
from IPython.core.autocall import IPyAutocall
import snappy
from urllib.request import pathname2url
from .gui import *

snappy_path = os.path.abspath(os.path.dirname(snappy.__file__))
icon_file = os.path.join(snappy_path, 'info_icon.gif')
debug_Tk = False
ansi_seqs = re.compile(r'(?:\x01*\x1b\[((?:[0-9]*;)*[0-9]*.)\x02*)*([^\x01\x1b]*)',
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

class Tk(Tk_.Tk):
    def __init__(self, error_handler=None):
        Tk_.Tk.__init__(self, className='snappy')
        # Tkinter ignores exceptions raised by callbacks, but
        # calls this function to report their occurrence.
        if error_handler:
            self.report_callback_exception = error_handler
        # In Python 2.7 the _default root does not get set correctly.
        if not Tk_._default_root:
            Tk_._default_root = self

#  Some ideas for the TkTerm class were borrowed from code written by
#  Eitan Isaacson, IBM Corp.

class TkTerm:
    """
    A Tkinter terminal window that runs an IPython shell.  This class
    supports the IOStream interface, and can function as a replacement
    for sys.stdout / sys.stderr.
    """
    def __init__(self, shell, name='TkTerm'):
        if debug_Tk:
            self.window = window = Tk()
            io.stdout = sys.stdout = self
        else:
            self.window = window = Tk(self.report_callback_exception)
            self.encoding = sys.stdout.encoding
            self.saved_io = (sys.stdout, sys.stderr)
            io.stdout = io.stderr = sys.stdout = sys.stderr = self
        self._input_buffer = ''
        self._current_indent = 0
        # Global Tk options
        window.option_add('*Menu.tearOff', 0)
        window.title(name)
        window.protocol("WM_DELETE_WINDOW", self.close)
        self.icon = Tk_.PhotoImage(file=icon_file)
        self.frame = frame = Tk_.Frame(window)
        self.text = text = Tk_.Text(frame,
                                    width=85,
                                    height=24,
                                    foreground='Black',
                                    borderwidth=3,
                                    background='#ec0fffec0',
                                    insertbackground='#ff0000',
                                    highlightthickness=0,
                                    relief=Tk_.FLAT
                                )
        self.set_font(Font(text, text.cget('font')))
        self.scroller = scroller = Tk_.Scrollbar(frame, command=text.yview)
        text.config(yscrollcommand=scroller.set)
        scroller.pack(side=Tk_.RIGHT, fill=Tk_.Y, pady=10)
        text.pack(fill=Tk_.BOTH, expand=Tk_.YES)
        frame.pack(fill=Tk_.BOTH, expand=Tk_.YES)
        text.focus_set()
        window.bind('<FocusIn>', lambda event=None: text.focus_set())
        text.bind('<KeyPress>', self.handle_keypress)
        text.bind('<KeyRelease>', self.handle_keyrelease)
        if not sys.platform == 'darwin':
            text.bind('<Control-c>', self.handle_keypress)
        text.bind('<Return>', self.handle_return)
        text.bind('<Shift-Return>', self.handle_shift_return)
        text.bind('<BackSpace>', self.handle_backspace)
        text.bind('<Delete>', self.handle_backspace)
        text.bind('<Tab>', self.handle_tab)
        text.bind('<Up>', self.handle_up)
        text.bind('<Shift-Up>', lambda event : None)
        text.bind('<Down>', self.handle_down)
        text.bind('<Shift-Down>', lambda event : None)
        text.bind('<Home>', self.go_to_beginning)
        text.bind('<<Cut>>', self.protect_text)
        text.bind('<<Copy>>', self.edit_copy)
        text.bind('<<Paste>>', self.edit_paste)
        text.bind('<<Clear>>', self.protect_text)
        text.bind('<ButtonRelease>', self.mouse_up)
        if sys.platform != 'darwin':
            text.bind_all('<ButtonPress-2>', self.middle_mouse_down)
            text.bind_all('<ButtonRelease-2>', self.middle_mouse_up)
        text.bind('<Button-3>', lambda event:'break')
        text.bind('<Button-4>', lambda event:text.yview_scroll(-1, Tk_.UNITS))
        text.bind('<Button-5>', lambda event:text.yview_scroll(1, Tk_.UNITS))
        if sys.platform == 'darwin':
            self.window.bind_all('<Command-Key-q>', self.close)
            self.window.bind_all('<Command-Key-Left>', self.go_to_beginning)
        elif sys.platform == 'linux2' or sys.platform == 'linux':
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
        self.multiline = False
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
        self.output_count = 0
        # Setup the IPython embedded shell
        self.IP = shell
        # IPython >= 5 uses this to write pygments tokenized strings (if it exists)
        shell.pt_cli = self
        shell.system = self.system
        sys.displayhook = shell.displayhook
        shell.more = False
        # Figure out the size of the first input prompt
        self._input_prompt()
        if shell.banner1:
            self.banner = shell.banner1
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
        # This flag will be set when printing a traceback
        self.showing_traceback = False
        self.closed = False
        self._saved_index = Tk_.END

    # Emulate a ListedWindow.  We are listed, even though we are unique.
    def bring_to_front(self):
        self.window.deiconify()
        self.window.lift()
        self.window.focus_force()

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
        traceback = '\n' + self.IP.InteractiveTB.stb2text(stb)
        if etype == KeyboardInterrupt:
            self.write('KeyboardInterrupt: ', style='msg')
            self.write('%s'%evalue)
        else:
            self.write(traceback)
        self.reset()
        self.showing_traceback = True
        self.send_code('\n')

    def SnapPea_callback(self, interrupted=False):
        """
        Callback for SnapPea to keep the UI alive during long computations.
        """
        if interrupted:
            self.interrupted = False
            raise KeyboardInterrupt('SnapPea computation aborted')

    def PARI_callback(self):
        """
        Callback for CyPari to keep the UI alive during long computations.
        """
        if self.running_code:
            if self.interrupted:
                snappy.pari.abort()
                self.interrupted = False
                raise KeyboardInterrupt('PARI computation aborted')

    def interrupt(self):
        if not self.running_code:
            source_raw = self.reset()
            self.IP.history_manager.store_inputs(self.IP.execution_count, source_raw)
            self.IP.execution_count += 1
            raise KeyboardInterrupt('Halted')
        else:
            self.interrupted = True
            # Inform the SnapPea kernel about the interrupt.
            snappy.SnapPy.SnapPea_interrupt()
            snappy.SnapPyHP.SnapPea_interrupt()

    def report_callback_exception(self, exc, value, traceback):
        """
        Called when exceptions are caught by Tk.
        """
        sys.last_type = exc
        sys.last_value = value
        sys.last_traceback = traceback
        self.IP.showtraceback()

    def set_font(self, fontdesc):
        self.text.config(font=fontdesc)
        normal_font = Font(self.text, self.text.cget('font'))
        self.bold_font = bold_font = Font(self.text, self.text.cget('font'))
        self.bold_font.config(weight='bold')
        self.char_size = normal_font.measure('M')
        text = self.text
        text.tag_config('output', font=normal_font)
        text.tag_config('Prompt', foreground='#0000cc', font=normal_font)
        text.tag_config('PromptNum', foreground='#0000bb', font=bold_font)
        text.tag_config('OutPrompt', foreground='#cc0000', font=normal_font)
        text.tag_config('OutPromptNum', foreground='#bb0000', font=bold_font)

    def close(self, event=None):
        self.window.quit()
        self.closed = True

    def handle_keypress(self, event):
        self.clear_completions()
        protected = self.text.compare(Tk_.INSERT, '<', 'output_end')
        # We only respond to ASCII keys: no accents, no emojis.
        if len(event.char) > 1:
            return
        keysym = event.keysym
        if keysym == 'Left':
            # Don't go past the prompt
            if self.text.compare(Tk_.INSERT, '<', 'output_end+1c'):
                return 'break'
            return
        # Check for a control character
        if event.state & 4:
            # Ctrl+C is an interrupt
            if keysym == 'c':
                self.interrupt()
            # Ctrl+Shift+C copies on all platforms, even macOS
            if keysym == 'C':
                self.edit_copy()
            # Ctrl+Shift+V pastes on all platforms, even macOS
            if keysym == 'V':
                self.edit_paste()
            # emacs shortcuts (Ctrl-k is built-in)
            if keysym == 'a':
                self.text.mark_set(Tk_.INSERT, 'output_end')
                return 'break'
            if event.keysym == 'e':
                self.text.mark_set(Tk_.INSERT, Tk_.END)
                return 'break'
            if event.keysym == 'u':
                self.text.delete('output_end', Tk_.END)
                return 'break'
        # space pages down when viewing protected output
        if keysym == 'space' and protected:
            self.page_down()
            return 'break'
        # Typing in the protected area should not do anything.
        if event.char and protected:
            self.text.tag_remove(Tk_.SEL, '1.0', Tk_.END)
            return 'break'
        # Apple Fn-left sends \0121 and means ^A
        if event.char == '\0121': # Apple Fn-Left
            self.text.mark_set(Tk_.INSERT, 'output_end')
            return 'break'

    def handle_keyrelease(self, event):
        if self.editing_hist and self.multiline:
            self.text.tag_add('history', 'output_end', Tk_.INSERT)

    def handle_return(self, event):
        self.clear_completions()
        line = self.text.get('output_end', Tk_.END)
        if self.editing_hist:
            insert = self.text.index(Tk_.INSERT)
            newline = self.text.index('output_end')
            tail = self.text.get('%s.0'%insert.split('.')[0], Tk_.END)
            if tail.strip() != '':
                    return
            prompt_size = self.text.index('output_end-1c').split('.')[1]
            first = self.text.index('output_end').split('.')[0]
            last = self.text.index(Tk_.END).split('.')[0]
            prompt = self._continuation_prompt(int(prompt_size)).pop()
            for line_number in range(int(first) + 1, int(last) - 1):
                self.write(prompt[1], style=prompt[0], advance=False,
                           mark='%s.0'%line_number)
            self.text.delete(newline+'-1c', newline)
        else:
            self.text.insert(Tk_.END, '\n')
        self.text.tag_add('output', 'output_end', Tk_.END)
        self.text.mark_set('output_end', Tk_.END)
        if not self.running_code:
            self.send_code(line)
        return 'break'

    def handle_shift_return(self, event):
        if self.editing_hist or self.hist_pointer == 0:
            return 'break'
        l1, c1 = map(int, self.text.index('output_end').split('.'))
        l2, c2 = map(int, self.text.index(Tk_.INSERT).split('.'))
        index = '%d.%d'%(l1 + 1, c2 - c1)
        self.text.delete('output_end', Tk_.END)
        self.write_history(force_multiline=True)
        self.text.mark_set(Tk_.INSERT, index)

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
            self.write(self.hist_stem.strip('\n'), style=(), advance=False)
            self.editing_hist = False
            self.text.tag_delete('history')
        else:
            self.write_history()
        return 'break'

    def handle_backspace(self, event):
        self.clear_completions()
        if self.text.compare(Tk_.INSERT, '<=', 'output_end'):
            self.window.bell()
            return 'break'
        if self._current_indent >= 4:
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

    def go_to_beginning(self, event):
        self.text.mark_set(Tk_.INSERT, 'output_end')
        return 'break'

    def do_completion(self, word, completion):
        tail = completion[len(word):]
        self.text.insert(self.tab_index, tail)
        self.tab_index = Tk_.END
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
        self.text.see(Tk_.END)

    def clear_completions(self):
        if self.tab_index:
            self.text.delete(self.tab_index, Tk_.END)
            self.tab_index = None
            self.tab_count = 0

    def stem(self, wordlist):
        if len(wordlist) == 1:
            return wordlist[0]
        result = ''
        for n in range(1,100):
            heads = set(w[:n] for w in wordlist)
            if len(heads) > 1:
                return result
            elif len(heads) == 1:
                result = heads.pop()
        return wordlist[0][:100]

    def write_history(self, force_multiline=False):
        self.text.see('output_end')
        input = self.filtered_hist[-self.hist_pointer]
        input = re.sub('\n+', '\n', input).rstrip()
        if input.find('\n') > 0 or force_multiline:
            self.editing_hist = True
            self.multiline = True
            # Add a newline to place the editing box below the prompt.
            self.write('\n')
            self.text.tag_config('history',
                                 background='White')
            self.text.tag_lower('history')
            self.write(input + '\n\n', style=('history',), advance=False)
            self.text.mark_set('history_end', Tk_.INSERT)
            self.text.mark_set(Tk_.INSERT, 'output_end')
        else:
            self.multiline = False
            self.write(input, style=(), advance=False)
        self.text.see(Tk_.INSERT)

    def protect_text(self, event):
        try:
            if self.text.compare(Tk_.SEL_FIRST, '<=', 'output_end'):
                self.window.bell()
                return 'break'
        except:
            pass

    def edit_actions(self):
        "Return a dictionary of allowable edit actions."
        result = {}
        try:
            selectable = self.text.compare(Tk_.SEL_FIRST, '<', Tk_.SEL_LAST)
        except:
            selectable = False
        try:
            protected = self.text.compare(Tk_.INSERT, '<', 'output_end')
        except:
            protected = True
        try:
            clip = self.text.clipboard_get()
        except:
            clip = None
        if selectable:
            result['Copy'] = self.edit_copy
        if clip:
            result['Paste'] = self.edit_paste
        if protected:
            return result
        else:
            if selectable:
                result['Delete'] = self.edit_delete
                result['Cut'] = self.edit_cut
        return result

    def edit_cut(self):
        try:
            self.text.clipboard_clear()
            self.text.clipboard_append(self.text.selection_get())
            self.text.delete(Tk_.SEL_FIRST, Tk_.SEL_LAST)
        except:
            pass

    def edit_copy(self, event=None):
        try:
            self.text.clipboard_clear()
            self.text.clipboard_append(self.text.selection_get())
            self.text.tag_remove(Tk_.SEL, '1.0', Tk_.END)
            if self.text.compare(Tk_.INSERT, '<', 'output_end'):
                self.text.mark_set(Tk_.INSERT, Tk_.END)
        except:
            pass

    def edit_paste(self, event=None):
        text = self.text
        try:
            clip = self.text.clipboard_get()
        except:
            return
        try:
            start = text.index(Tk_.SEL_FIRST)
            text.delete(Tk_.SEL_FIRST, Tk_.SEL_LAST)
            text.insert(Tk_.SEL_FIRST, clip)
        except:
            if self.text.compare(Tk_.INSERT, '<', 'output_end'):
                self.text.mark_set(Tk_.INSERT, 'output_end')
            text.insert(Tk_.INSERT, clip)
            self.text.see(Tk_.INSERT)
        try:
            self.text.tag_remove(Tk_.SEL, Tk_.SEL_FIRST, Tk_.SEL_LAST)
        except:
            pass
        # prevent multiple pastes with ^V
        return 'break'

    def edit_delete(self):
        try:
            self.text.delete(Tk_.SEL_FIRST, Tk_.SEL_LAST)
        except:
            pass

    def mouse_up(self, event):
        self.text.mark_set(Tk_.INSERT, Tk_.END)

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
        self.text.insert(Tk_.END, '\n')
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

    def _input_prompt(self):
        result = [('Prompt', 'In['),
                  ('PromptNum', '%d'%self.IP.execution_count),
                  ('Prompt', ']: ')]
        self._prompt_size = sum(len(token[1]) for token in result)
        return result

    def _continuation_prompt(self, size):
        prompt_text = ' '*(size - 5) + '...: '
        return [( 'Prompt', prompt_text)]

    def interact_prompt(self):
        """
        Print an input prompt or a continuation prompt.
        """
        if self.showing_traceback:
            self.showing_traceback = False
            return
        try:
            if self.IP.more:
                prompt_tokens = self._continuation_prompt(self._prompt_size)
            else:
                self.write(self.IP.separate_in)
                prompt_tokens = self._input_prompt()
        except:
            self.IP.showtraceback()
        for style, text in prompt_tokens:
            self.write(text, style)

    def interact_handle_input(self, code, script=False):
        transformer = self.IP.input_transformer_manager
        assert code.endswith('\n')
        code = re.sub('\n+', '\n', code).rstrip() + '\n'
        self._input_buffer += code
        status, indent = transformer.check_complete(self._input_buffer)
        if script:
            # We are running a script.  If adding a newline would make the
            # buffer complete and the indent 0 then we want to run the code.
            # So add the newline to force the buffer to be run.
            Xstatus, Xindent = transformer.check_complete(self._input_buffer + '\n')
            if Xstatus == 'complete' and Xindent is None:
                self._input_buffer += '\n'
                status, indent = Xstatus, Xindent
        if status == 'incomplete':
            self.IP.more = True
            self._current_indent = indent or 0
            return
        if status == 'invalid':
            # Force display of the SyntaxError
            self.IP.run_cell(self._input_buffer, store_history=True)
            self.reset()
            return
        # The code is complete, but we only run it if the indent level is 0 or
        # if there is an empty line at the end.
        if self._current_indent == 0 or self._input_buffer.endswith('\n\n'):
            self.editing_hist = False
            self.multiline = False
            self.text.tag_delete('history')
            self._input_buffer = self._input_buffer.rstrip() + '\n'
            if self._input_buffer.count('\n') > 1:
                self.write('\n')
            self.running_code = True
            result = self.IP.run_cell(self._input_buffer, store_history=True)
            self.write('\n')
            self.reset()
        else:
            self.IP.more = True

    def reset(self):
        result = self._input_buffer
        self._input_buffer = ''
        self._current_indent = 0
        self.text.delete('output_end',Tk_.INSERT)
        self.editing_hist = False
        self.multiline = False
        self.running_code = False
        self.IP.more = False

        self.text.tag_delete('history')
        return result

    def send_code(self, code):
        """
        Accumulate some lines of code and either issue a continuation
        prompt or execute the code, print the result and issue a new
        input prompt.
        """
        try:
            self.interact_handle_input(code)
        except KeyboardInterrupt:
            self.write('(IP) Keyboard Interrupt: ')
            self.reset()
        self.interact_prompt()
        self.text.see(Tk_.INSERT)
        self.text.mark_set('output_end',Tk_.INSERT)
        if self.IP.more:
            self.text.insert(Tk_.INSERT, ' '*self._current_indent, ())
        self.hist_pointer = 0
        self.hist_stem = ''

    def write(self, string, style=('output',), advance=True,
              mark='output_end', see=True):
        """
        Writes a string containing ansi color escape sequences to our
        Text widget, starting at the specified mark.  If the advance
        option is True, advance the mark to the end of the output.
        """
        if self.quiet:
            return
        #if self.interrupted:
        #    self.interrupted = False
        #    raise KeyboardInterrupt('Writing')
        self.text.mark_set(Tk_.INSERT, mark)
        pairs = ansi_seqs.findall(string)
        for pair in pairs:
            code, text = pair
            tags = (code,) + style if code else style
            if text:
                self.text.insert(Tk_.INSERT, text, tags)
                self.output_count += len(text)
        if advance:
            self.text.mark_set(mark, Tk_.INSERT)
        if see:
            self.text.see(Tk_.INSERT)
        # Give the Text widget a chance to update itself every
        # so often (but let's not overdo it!)
        if self.output_count > 2000:
            self.output_count = 0
            self.text.update_idletasks()
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

    def writelines(self, lines):
        lines = iter(lines)
        for line in lines:
            self.write(line)

    def page(self, text):
        """
        Our pager.  Just writes the text to the screen then jumps back to
        the beginning.  You can scroll down to read it.
        """
        index = self.text.index(Tk_.INSERT)
        self.write('\n'+str(text), see=False)
        self.window.after_idle(self.text.see, index)

    def page_down(self):
        insert_line = int(str(self.text.index(Tk_.INSERT)).split('.')[0])
        prompt_line = int(str(self.text.index('output_end')).split('.')[0])
        height = prompt_line - insert_line
        if height < 1:
            return
        scroll_amount = min(int(self.text.cget('height')) - 1, height)
        self.text.yview_scroll(scroll_amount, Tk_.UNITS)
        # if scroll_amount == height:
        #     self.text.mark_set(Tk_.INSERT, 'output_end')
        # else:
        #     self.text.mark_set(Tk_.INSERT, '%s.0'%(insert_line + scroll_amount))


    def flush(self):
        """
        Required for a stdout / stderr proxy.
        """
        self.text.update_idletasks()
