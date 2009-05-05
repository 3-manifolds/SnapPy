import IPython
import Tkinter as Tk_
import os, sys
import re

DefaultFonts = {'darwin': ('Monaco', 18, 'normal'),
                'linux2': ('fixed', 18, 'normal')
                }
def default_font():
    try:
        return DefaultFonts[sys.platform]
    except:
        return 'Helvetica'

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
                '1;30m': 'DarkGray',
                '1;31m': 'DarkRed',
                '1;32m': 'SeaGreen',
                '1;33m': 'Yellow',
                '1;34m': 'LightBlue',
                '1;35m': 'MediumPurple',
                '1;36m': 'LightCyan',
                '1;37m': 'White'}

class TkTerm:
    """
    A Tkinter terminal window that runs an IPython shell.
    Some ideas borrowed from code written by Eitan Isaacson, IBM Corp.
    """
    def __init__(self, the_shell, name='TkTerm', root=None):
        self.banner = the_shell.banner
        self.window = window = Tk_.Tk(root)
        window.protocol("WM_DELETE_WINDOW", self.close)
        self.frame = frame = Tk_.Frame(window)
        self.text = text = Tk_.Text(frame,
                                font=default_font()
                                )
        self.scroller = scroller = Tk_.Scrollbar(frame, command=text.yview)
        text.config(yscrollcommand = scroller.set)
        scroller.pack(side=Tk_.RIGHT, fill=Tk_.Y, pady=10)
        text.pack(fill=Tk_.BOTH, expand=Tk_.YES)
        frame.pack(fill=Tk_.BOTH, expand=Tk_.YES)
        text.focus_set()
        text.bind('<KeyPress>', self.handle_keypress)
        text.bind('<Return>', self.handle_return)
        text.bind('<BackSpace>', self.handle_backspace)
        text.bind('<Delete>', self.handle_backspace)
        text.tag_config('output', background='White')
        for code in ansi_colors:
            text.tag_config(code, foreground=ansi_colors[code])
        self.end_index = self.text.index(Tk_.INSERT)
        self.banner = the_shell.banner
        self.IP = the_shell.IP
        # This stuff should not all be necessary, but it is.
        self.IP.write = self.write
        self.IP.write_err = self.write
        IPython.Shell.Term.cout = self
        IPython.Shell.Term.cerr = self
        sys.stdout = self
        sys.stderr = self
        sys.displayhook = self.IP.outputcache
        #
        self.start_interaction()

    def close(self):
        self.live = False
        self.window.update_idletasks()
        self.window.destroy()

    def close_event(self, event):
        self.close()

    def handle_keypress(self, event):
        if event.char == '\003':
            raise KeyboardInterrupt
        if event.char == '\004':
            self.close()
            return
        if self.text.compare(Tk_.INSERT, '<', self.end_index):
            self.text.mark_set(Tk_.INSERT, self.end_index)

    def handle_return(self, event):
        line=self.text.get(self.end_index, Tk_.INSERT)
        self.text.tag_add('output', self.end_index, Tk_.INSERT)
        self.end_index = Tk_.INSERT
        self.send_line(line)
        return 'break'

    def handle_backspace(self, event):
        if self.text.compare(Tk_.INSERT, '<=', self.end_index):
            self.window.bell()
            return 'break'

    def start_interaction(self):
        """
        Print the interpreter's banner and issue the first prompt.
        """
        self.write(self.banner)
        self.IP.interact_prompt()
        self.end_index = self.text.index(Tk_.INSERT)
 
    def send_line(self, line):
        """
        Send one line of input to the interpreter, which will write
        the result on our Text widget.  Then issue a new prompt.
        """
        self.write('\n')
        if self.IP.more:
            self.IP.rl_do_indent = True
        else:
            self.IP.rl_do_indent = False
        line = line.decode(self.IP.stdin_encoding)
        self.IP.interact_handle_input(line)
        self.IP.interact_prompt()
        self.text.see(Tk_.INSERT)
        self.end_index = self.text.index(Tk_.INSERT)

    def write(self, string):
        """
        Writes a string containing ansi color escape sequences to our
        Text widget.
        """
        pairs = ansi_seqs.findall(string)
        for pair in pairs:
            code, text = pair
            tags = (code, 'output') if code else ('output',)
            if text:
                self.text.insert(Tk_.INSERT, text, tags)
        self.end_index = self.text.index(Tk_.INSERT)

    def flush(self):
        """
        Since we are pretending to be an IOTerm.
        """
        pass

if __name__ == "__main__":
    import SnapPy
    from SnapPy.SnapPy_shell import the_shell
    SnapPy_ns = dict([(x, getattr(SnapPy,x)) for x in SnapPy.__all__])
    the_shell.IP.user_ns.update(SnapPy_ns)
    os.environ['TERM'] = 'dumb'
    terminal = TkTerm(the_shell)
    terminal.window.mainloop()
