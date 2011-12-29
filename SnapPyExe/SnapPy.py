import snappy
import os, IPython

#----------------------------------------------------------------------------
class HomeDirError(Exception):
    pass

def get_home_dir():
    """Return the closest possible equivalent to a 'home' directory.

    We first try $HOME.  Absent that, on NT it's $HOMEDRIVE\$HOMEPATH.

    Currently only Posix and NT are implemented, a HomeDirError exception is
    raised for all other OSes. """

    isdir = os.path.isdir
    env = os.environ

    # first, check py2exe distribution root directory for _ipython.
    # This overrides all. Normally does not exist.

    print "HACK HACK HACK the HACK IS BACK!"
    # NMD 2010/2/3: HACK IS RIGHT HERE!  
    if False and hasattr(sys, "frozen"): #Is frozen by py2exe
        if '\\library.zip\\' in IPython.__file__.lower():#libraries compressed to zip-file
            root, rest = IPython.__file__.lower().split('library.zip')
        else: 
            root=os.path.join(os.path.split(IPython.__file__)[0],"../../")
        root=os.path.abspath(root).rstrip('\\')
        if isdir(os.path.join(root, '_ipython')):
            os.environ["IPYKITROOT"] = root
        return root
    try:
        homedir = env['HOME']
        if not isdir(homedir):
            # in case a user stuck some string which does NOT resolve to a
            # valid path, it's as good as if we hadn't foud it
            raise KeyError
        return homedir
    except KeyError:
        if os.name == 'posix':
            raise HomeDirError,'undefined $HOME, IPython can not proceed.'
        elif os.name == 'nt':
            # For some strange reason, win9x returns 'nt' for os.name.
            try:
                homedir = os.path.join(env['HOMEDRIVE'],env['HOMEPATH'])
                if not isdir(homedir):
                    homedir = os.path.join(env['USERPROFILE'])
                    if not isdir(homedir):
                        raise HomeDirError
                return homedir
            except KeyError:
                try:
                    # Use the registry to get the 'My Documents' folder.
                    import _winreg as wreg
                    key = wreg.OpenKey(wreg.HKEY_CURRENT_USER,
                                       "Software\Microsoft\Windows\CurrentVersion\Explorer\Shell Folders")
                    homedir = wreg.QueryValueEx(key,'Personal')[0]
                    key.Close()
                    if not isdir(homedir):
                        e = ('Invalid "Personal" folder registry key '
                             'typically "My Documents".\n'
                             'Value: %s\n'
                             'This is not a valid directory on your system.' %
                             homedir)
                        raise HomeDirError(e)
                    return homedir
                except HomeDirError:
                    raise
                except:
                    return 'C:\\'
        elif os.name == 'dos':
            # Desperate, may do absurd things in classic MacOS. May work under DOS.
            return 'C:\\'
        else:
            raise HomeDirError,'support for your operating system not implemented.'


IPython.utils.path.get_home_dir = get_home_dir

from snappy.app import main
main()
