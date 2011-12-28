import tkFileDialog, sys

askopenfile = tkFileDialog.askopenfile

def asksaveasfile(mode='w',**options):
    """
    Ask for a filename to save as, and returned the opened file.
    Modified from tkFileDialog to more intelligently handle
    default file extensions. 
    """
    if sys.platform == 'darwin':
        if options.has_key('defaultextension') and not options.has_key('initialfile'):
            options['initialfile'] = 'untitled' + options['defaultextension']

    return tkFileDialog.asksaveasfile(mode=mode, **options)

if __name__ == '__main__':
    asksaveasfile(defaultextension='.txt')

    
