# -*- coding: utf-8 -*-
import os
import sys
import re
from urllib.request import pathname2url
import tkinter as Tk_
from tkinter import ttk
from tkinter.font import Font
from tkinter.font import families as font_families
from tkinter.messagebox import askyesno
from IPython.utils import io

debug_Tk = True
ansi_seqs = re.compile(r'(?:\x01*\x1b\[((?:[0-9]*;)*[0-9]*.)\x02*)*([^\x01\x1b]*)',
                       re.MULTILINE)
ansi_colors = {'0;30m': 'Black',
               '0;31m': 'Red',
               '0;32m': 'Green',
               '0;33m': 'Brown',
               '0;34m': 'Blue',
               '0;35m': 'Purple',
               '0;36m': 'Cyan',
               '0;37m': 'LightGray',
               '1;30m': 'Black',  # 'DarkGray',
               '1;31m': 'DarkRed',
               '1;32m': 'SeaGreen',
               '1;33m': 'Yellow',
               '1;34m': 'LightBlue',
               '1;35m': 'MediumPurple',
               '1;36m': 'LightCyan',
               '1;37m': 'White'}
delims = re.compile(r'[\s\[\]\{\}\(\)\+\-\=\'`~!@#\$\^\&\*]+')


class FontChoice:
    def __init__(self, family, size, weight, slant):
        self.family = family
        self.size = size
        self.weight = weight
        self.slant = slant
        self.rest = f'{self.weight} {self.slant}'

    def as_tuple(self):
        size = self.size
        if sys.platform == 'darwin' and Tk_.TkVersion >= 9.0:
            size = int(size/1.3)
        return (self.family, size, self.rest)

    def __repr__(self):
        return 'FontChoice' + repr((self.family, self.size, self.weight, self.slant))

    def bold(self):
        return FontChoice(self.family, self.size, 'bold', self.slant)


def default_terminal_font():
    size = 13 if sys.platform == 'darwin' else 11
    family = Font(font='TkFixedFont').actual()['family']
    if sys.platform == 'win32':
        # Default is Courier New which is ugly and appears blurry.
        available = font_families()
        for better in ['Consolas', 'Cascadia Mono SemiLight']:
            if better in available:
                family = better

    return FontChoice(family, size, 'normal', 'roman')


class Tk(Tk_.Tk):
    def __init__(self, error_handler=None):
        Tk_.Tk.__init__(self, className='snappy')
        # Tkinter ignores exceptions raised by callbacks, but
        # calls this function to report their occurrence.
        if error_handler:
            self.report_callback_exception = error_handler


class TkTerminalBase:
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
            io.stdout = io.stderr = sys.stdout = sys.stderr = self
        self._input_buffer = ''
        self._current_indent = 0
        # Global Tk options
        window.option_add('*Menu.tearOff', 0)
        window.title(name)
        window.protocol("WM_DELETE_WINDOW", self.close)
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
        self.set_font(default_terminal_font())
        self.scroller = scroller = Tk_.Scrollbar(frame, command=text.yview)
        text.config(yscrollcommand=scroller.set)
        scroller.pack(side=Tk_.RIGHT, fill=Tk_.Y, pady=10)
        text.pack(fill=Tk_.BOTH, expand=Tk_.YES)
        frame.pack(fill=Tk_.BOTH, expand=Tk_.YES)
        text.focus_set()
        window.bind('<FocusIn>', lambda event=None: text.focus_set())
        text.bind('<Control-c>', self.handle_control_c)
        text.bind('<KeyPress>', self.handle_keypress)
        text.bind('<KeyRelease>', self.handle_keyrelease)
        text.bind('<Return>', self.handle_return)
        text.bind('<Shift-Return>', self.handle_shift_return)
        text.bind('<BackSpace>', self.handle_backspace)
        text.bind('<Delete>', self.handle_backspace)
        text.bind('<Tab>', self.handle_tab)
        text.bind('<Up>', self.handle_up)
        text.bind('<Down>', self.handle_down)
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
            self.text.bind('<Command-Key-Up>', self.jump_up)
        elif sys.platform == 'linux2' or sys.platform == 'linux':
            self.window.bind_all('<Alt-Key-q>', self.close)
        self.add_bindings()
        # 'output_end' marks the end of the text written by us.
        # Everything above this position should be
        # immutable, and tagged with the "output" style.
        # Normally it is set at the end of the input prompt.
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
            self.banner = "Python %s on %s\n%s\n(%s)\n" % (
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
        # This flag is set to prevent quitting the app.
        self.blockers = {}
        self.can_quit = True
        self.close_callback = lambda :None

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
            self.write('%s' % evalue)
        else:
            self.write(traceback)
        self.reset()
        self.showing_traceback = True
        self.process_return('\n')

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
        else:
            self.interrupted = True
            # Inform the SnapPea kernel about the interrupt.
            snappy.SnapPy.SnapPea_interrupt()
            snappy.SnapPyHP.SnapPea_interrupt()
        # For some reason we need to raise a KeyboardInterrupt here.  Otherwise
        # the call to KeyboardInterrupt('Halted') will be ignored.
        raise KeyboardInterrupt('Interrupted')

    def report_callback_exception(self, exc, value, traceback):
        """
        Called when exceptions are caught by Tk.
        """
        sys.last_type = exc
        sys.last_value = value
        sys.last_traceback = traceback
        self.IP.showtraceback()

    def set_font(self, font_choice):
        normal_tuple = font_choice.as_tuple()
        bold_tuple = font_choice.bold().as_tuple()
        normal_font = Font(font=font_choice.as_tuple())
        self.char_size = normal_font.measure('M')
        self.text.config(font=normal_tuple)
        text = self.text
        text.tag_config('output', font=normal_tuple)
        text.tag_config('Prompt', foreground='#0000cc', font=normal_tuple)
        text.tag_config('PromptNum', foreground='#0000bb', font=bold_tuple)
        text.tag_config('OutPrompt', foreground='#cc0000', font=normal_tuple)
        text.tag_config('OutPromptNum', foreground='#bb0000', font=bold_tuple)

    def add_blocker(self, window, message):
        self.blockers[window] = message

    def remove_blocker(self, window):
        self.blockers.pop(window)

    def close(self, event=None):
        can_quit = True
        topmost = None
        for blocker in self.blockers:
            if blocker.attributes('-topmost'):
                topmost = blocker
                blocker.attributes('-topmost', False)
            message = self.blockers[blocker]
            answer  = askyesno('Quit',
                message + '\nDo you want to quit anyway?')
            if not answer:
                can_quit = False
                break;
        if topmost:
            topmost.attributes('-topmost', True)
        if can_quit:
            for blocker in self.blockers:
                blocker.destroy()
            self.window.quit()
            self.closed = True

    def handle_control_c(self, event):
        self.interrupt()
        self.interrupted = False

    def handle_keypress(self, event):
        self.clear_completions()
        protected = self.text.compare(Tk_.INSERT, '<', 'output_end')
        # We only respond to ASCII keys: no accents, no emojis.
        if len(event.char) > 1:
            return
        keysym = event.keysym
        if keysym == 'Left':
            # Don't go into the prompt or continuation prompt
            line, pos = map(int, self.text.index(Tk_.INSERT).split('.'))
            if pos <= self._prompt_size:
                prompt_line = int(self.text.index('output_end').split('.')[0])
                insert_line = int(self.text.index(Tk_.INSERT).split('.')[0])
                if line > prompt_line:
                    self.text.mark_set(Tk_.INSERT, '%s.0-1c' % insert_line)
                else:
                    self.window.bell()
                return 'break'
            return
        if keysym == 'Right':
            # Don't go into a continuation prompt or past the end
            insert_pos = int(self.text.index('%s+1c' % Tk_.INSERT).split('.')[1])
            if insert_pos < self._prompt_size:
                insert_line = int(self.text.index(Tk_.INSERT).split('.')[0])
                if self.text.compare(Tk_.INSERT, '>=', '%s-2c' % Tk_.END):
                    self.window.bell()
                else:
                    self.text.mark_set(Tk_.INSERT, '%d.%d' % (
                        insert_line + 1, self._prompt_size))
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

    def handle_return(self, event=None):
        # If the input consists of one complete line of code we run it,
        # regardless of where the insert cursor is located.  Otherwise we only
        # run the code if the cursor is at the end of the input
        self.clear_completions()
        if self.running_code:
            return 'break'
        cursor = int(self.text.index(Tk_.INSERT).split('.')[0])
        last = int(self.text.index(Tk_.END).split('.')[0])
        first = int(self.text.index('output_end').split('.')[0])
        if cursor == first == last - 1:  # single line input
            self.text.mark_set(Tk_.INSERT, Tk_.END)
        self.text.insert(Tk_.INSERT, '\n')
        cell = self.text.get('output_end', Tk_.INSERT)
        try:
            self.interact_handle_input(cell)
        except KeyboardInterrupt:
            self.write('(IP) Keyboard Interrupt: ')
            self.reset()
            self.hist_pointer = 0
            self.hist_stem = ''
        self.interact_prompt()
        self.text.see(Tk_.INSERT)
        if self.IP.more:
            self.text.insert(Tk_.INSERT, ' '*self._current_indent, ())
        self.hist_pointer = 0
        self.hist_stem = ''
        return 'break'

    def handle_shift_return(self, event):
        self.text.mark_set(Tk_.INSERT, Tk_.END)
        return self.handle_return()

    def jump_up(self, event):
        return self.handle_up(event, jump=True)

    def handle_up(self, event, jump=False):
        if self.text.compare(Tk_.INSERT, '<', 'output_end'):
            return
        insert_line_number = str(self.text.index(Tk_.INSERT)).split('.')[0]
        prompt_line_number = str(self.text.index('output_end')).split('.')[0]
        if insert_line_number != prompt_line_number:
            return
        if self.hist_pointer == 0:
            input_history = self.IP.history_manager.input_hist_raw
            self.hist_stem = self.text.get('output_end', Tk_.END).strip()
            self.filtered_hist = [x for x in input_history
                                  if x and x.startswith(self.hist_stem)]
        if self.hist_pointer >= len(self.filtered_hist):
            self.window.bell()
            return 'break'
        self.text.delete('output_end', Tk_.END)
        self.hist_pointer += 1
        self.write_history()
        if jump:
            self.text.mark_set(Tk_.INSERT, 'output_end')
        return 'break'

    def handle_down(self, event):
        if self.text.compare(Tk_.INSERT, '<', 'output_end'):
            return
        if self.editing_hist:
            insert_line = int(str(self.text.index(Tk_.INSERT)).split('.')[0])
            bottom_line = int(str(self.text.index('history_end')).split('.')[0])
            if insert_line < bottom_line:
                return
        if self.hist_pointer == 0:
            return
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
        line, pos = map(int, self.text.index(Tk_.INSERT).split('.'))
        first_line = int(self.text.index('output_end').split('.')[0])
        if pos <= self._prompt_size + 1 and line != first_line:
            start = '%d.end' % (line - 1)
            end = '%d.end' % (line)
            self.text.mark_set(Tk_.INSERT, start)
            self.text.delete(start, '%s.%s' % (line, self._prompt_size))
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
            stem, completions = self.IP.complete(word)
        except TypeError:
            completions = []
        if word.find('_') == -1:
            completions = [x for x in completions
                           if x.find('__') == -1 and x.find('._') == -1]
        # No meaningful completions.  Ring the bell.
        if len(completions) == 0 or len(completions) == 1 and completions[0] == stem:
            self.window.bell()
            self.tab_count = 0
            return 'break'
        # Only one completion.  Use it.
        if len(completions) == 1:
            self.do_completion(stem, completions[0])
            self.tab_count = 0
            return 'break'
        max_stem = self.max_stem(stem, completions)
        # Add the maximal stem of all completions if it extends the word,
        if len(max_stem) > len(word):
            self.do_completion(stem, max_stem)
            self.tab_count = 0
            return 'break'
        # Show the possible completions, with a warning if there are lots.
        if len(completions) > 60 and self.tab_count == 1:
            self.show_completions('', '',
                ['%s possibilities -- hit tab again to view them all' %
                len(completions)])
        else:
            self.show_completions(word, stem, completions)
            if len(completions) <= 60:
                self.tab_count += 1
        return 'break'

    def go_to_beginning(self, event):
        self.text.mark_set(Tk_.INSERT, 'output_end')
        return 'break'

    def do_completion(self, stem, completion):
        tail = completion[len(stem):]
        self.text.insert(self.tab_index, tail)
        self.tab_index = Tk_.END
        self.tab_count = 0

    def show_completions(self, word, stem, comps):
        n = len(stem)
        comps = [word + c[n:] for c in comps]
        self.text.delete(self.tab_index, Tk_.END)
        width = self.text.winfo_width()
        charwidth = width // self.char_size
        biggest = 2 + max([len(x) for x in comps])
        num_cols = max(charwidth // biggest, 1)
        num_rows = (len(comps) + num_cols - 1) // num_cols
        rows = []
        format = '%%-%ds' % biggest
        for n in range(num_rows):
            rows.append(''.join(format % x for x in comps[n:len(comps):num_rows]))
        view = '\n'.join(rows)
        self.text.insert(self.tab_index, '\n'+view)
        self.text.mark_set(Tk_.INSERT, self.tab_index)
        self.text.see(Tk_.END)

    def clear_completions(self):
        if self.tab_index:
            self.text.delete(self.tab_index, Tk_.END)
            self.tab_index = None
            self.tab_count = 0

    def max_stem(self, stem, completions):
        if len(completions) == 1:
            return completions[0]
        result = stem
        for n in range(len(stem) + 1, 100):
            heads = {w[:n] for w in completions}
            if len(heads) > 1:
                return result
            elif len(heads) == 1:
                result = heads.pop()
        return result

    def write_continuation_prompt(self):
        prompt_tokens = self._continuation_prompt(self._prompt_size)
        for style, text in prompt_tokens:
            self.write(text, style, mark=Tk_.INSERT, advance=True)

    def write_history(self, force_multiline=False):
        self.text.see('output_end')
        input = self.filtered_hist[-self.hist_pointer]
        input = re.sub('\n+', '\n', input).rstrip()
        lines = input.split('\n')
        if len(lines) > 1 or force_multiline:
            self.editing_hist = True
            self.multiline = True
            self.text.delete('output_end', Tk_.END)
            self.write(lines[0] + '\n', mark=Tk_.INSERT, advance=True)
            prompt_tokens = self._continuation_prompt(self._prompt_size)
            for line in lines[1:-1]:
                self.write_continuation_prompt()
                self.write(line + '\n', mark=Tk_.INSERT)
            self.write_continuation_prompt()
            self.write(lines[-1], mark=Tk_.INSERT)
            self.text.mark_set('history_end', Tk_.INSERT)
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
        Display a banner and prepare to begin interaction.
        """
        # Subclasses should override this method
        banner_label = Tk_.Label(self.text,
            text="Please override the start_interaction method.",
            anchor=Tk_.W,
            justify=Tk_.LEFT)
        self.text.window_create(Tk_.END, window=banner_label)
        self.text.insert(Tk_.END, '\n')
        self.text.mark_set('output_end', '2.0')

    def _input_prompt(self):
        result = [('Prompt', 'In['),
                  ('PromptNum', '%d' % self.IP.execution_count),
                  ('Prompt', ']: ')]
        self._prompt_size = sum(len(token[1]) for token in result)
        return result

    def _continuation_prompt(self, size):
        prompt_text = ' '*(size - 5) + '...: '
        return [( 'Prompt', prompt_text)]

    def interact_prompt(self):
        """
        Print an input prompt or a continuation prompt.  For an input
        prompt set the output_end mark at the end of the prompt.
        """
        if self.showing_traceback:
            self.showing_traceback = False
            return
        try:
            if self.IP.more or self.editing_hist:
                self.write_continuation_prompt()
            else:
                if int(self.text.index('output_end').split('.')[1]) != 0:
                    self.write('\n\n', mark='output_end')
                prompt_tokens = self._input_prompt()
                for style, text in prompt_tokens:
                    self.write(text, style, mark='output_end')
        except:
            self.IP.showtraceback()

    def clean_code(self, code):
        """
        Remove blank lines and continuation prompts.
        """
        if not code.strip():
            return '\n'
        lines = list(code.split('\n'))
        clean_lines = [lines[0].lstrip()]
        for line in lines[1:]:
            try:
                clean_lines.append(line.split(': ', 1)[1])
            except IndexError:
                pass
        return '\n'.join(clean_lines)

    def interact_handle_input(self, cell, script=False):
        """
        Validate the code in the cell.  If the code is valid but
        incomplete, set the indent that should follow a continuation
        prompt and set the 'more' flag.  If the code is valid and
        complete then run the code.
        """
        assert cell.endswith('\n')
        if not cell.strip():
            self._current_indent = 0
            return
        if script:
            self._input_buffer += cell
        else:
            self._input_buffer = self.clean_code(cell)
        transformed_cell = self.IP.transform_cell(self._input_buffer)
        status, indent = self.IP.check_complete(transformed_cell)
        self._current_indent = indent or 0
        if status == 'incomplete':
            self.IP.more = True
            return
        if status == 'invalid':
            # Force display of the SyntaxError
            self.text.mark_set(Tk_.INSERT, self.text.index('output_end-1line'))
            self.IP.run_cell(self._input_buffer, store_history=True)
            self.text.insert('output_end-1line', '\n')
            self.reset()
            self.text.delete('output_end', Tk_.END)
            return
        # The code is complete, but we only run it if the cursor is on the
        # last line of the input and that line is blank.
        insert_line = int(self.text.index(Tk_.INSERT).split('.')[0])
        prompt_line = int(self.text.index('output_end').split('.')[0])
        tail = self.text.get('%d.%d' % (insert_line, self._prompt_size), Tk_.END)
        if not tail.strip():
            self.text.tag_delete('history')
            self._input_buffer = self._input_buffer.rstrip() + '\n'
            self.text.delete(Tk_.INSERT, Tk_.END)
            self.text.insert(Tk_.INSERT, '\n')
            self.text.mark_set('output_end', Tk_.INSERT)
            self.multiline = False
            self.running_code = True
            self.editing_hist = False
            last_line = insert_line - 1
            if last_line > prompt_line:
                # Delete the last continuation prompt.
                self.text.delete('%d.0' % last_line, '%d.0 lineend' % last_line)
            self.IP.run_cell(self._input_buffer, store_history=True)
            # Add a newline after the output.
            self.write('\n')
            self.text.tag_add('output', 'output_end', Tk_.END)
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

    def write(self, string, style=('output',), advance=True,
              mark='output_end', see=True):
        """
        Writes a string containing ansi color escape sequences to our
        Text widget, starting at the specified mark.  If the advance
        option is True, advance the mark to the end of the output.
        """
        if self.quiet:
            return
        # if self.interrupted:
        #    self.interrupted = False
        #    raise KeyboardInterrupt('Writing')
        if mark != Tk_.INSERT:
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
        # so often.  In order to force processing of key events
        # we have to call the evil update method, unfortunately.
        if self.output_count > 2000:
            self.output_count = 0
            self.text.update()
            if self.interrupted:
                self.interrupted = False
                raise KeyboardInterrupt('Halted')

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
