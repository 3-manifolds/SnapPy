import sys, platform
import tkinter.filedialog as tkFileDialog

def askopenfile(**options):
    if sys.platform == 'darwin' and platform.mac_ver()[0] < '10.15.2':
        options.pop('parent', None)
    return tkFileDialog.askopenfile(**options)

def asksaveasfile(mode='w',**options):
    """
    Ask for a filename to save as, and returned the opened file.
    Modified from tkFileDialog to more intelligently handle
    default file extensions.
    """
    if sys.platform == 'darwin':
        if platform.mac_ver()[0] < '10.15.2':
            options.pop('parent', None)
        if 'defaultextension' in options and not 'initialfile' in options:
            options['initialfile'] = 'untitled' + options['defaultextension']

    return tkFileDialog.asksaveasfile(mode=mode, **options)

if __name__ == '__main__':
    asksaveasfile(defaultextension='.txt')
