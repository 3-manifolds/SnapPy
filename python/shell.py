# -*- coding: utf-8 -*-
# SnapPy's subclass of the IPython InteractiveShellEmbed
from IPython.core.displayhook import DisplayHook
try:
    from IPython.terminal.embed import InteractiveShellEmbed
except ImportError:
    from IPython.frontend.terminal.embed import InteractiveShellEmbed

class SnapPyPromptDisplayHook(DisplayHook):
    """
    A DisplayHook used when displaying SnapPy's output prompts.  This
    subclass overrides one method in order to write the output prompt
    into the SnapPy console instead of sys.stdout.
    """

    def write_output_prompt(self):
        output = self.shell.output
        output.write(self.shell.separate_out)
        # If we're not displaying a prompt, it effectively ends with a newline,
        # because the output will be left-aligned.
        self.prompt_end_newline = True
        if self.do_full_cache:
            tokens = self.shell.prompts.out_prompt_tokens()
            prompt_txt = ''.join(s for t, s in tokens)
            if prompt_txt and not prompt_txt.endswith('\n'):
                # Ask for a newline before multiline output
                self.prompt_end_newline = False
            for token, text in tokens:
                output.write(text, style=(token[0],))

class SnapPyInteractiveShellEmbed(InteractiveShellEmbed):
    """
    An embedded IPython shell which can use a TkTerm for all of its
    input and output, including the output prompt.
    """
    readline_use = False
    autoindent_use = False
    colors_force = True

    def __init__(self, *args, **kwargs):
        super(InteractiveShellEmbed, self).__init__(*args, **kwargs)
        # Currently, using jedi means no completions at all.
        self.Completer.use_jedi = False

    def _displayhook_class_default(self):
        return SnapPyPromptDisplayHook
