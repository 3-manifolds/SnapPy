import os, subprocess

sage_root = os.environ['SAGE_ROOT']
sage_local = os.path.join(sage_root, 'local')
sage_libdynload = os.path.join(sage_local, 'lib', 'python3.7', 'lib-dynload')
old_libssl = '/usr/local/lib/libssl.1.1.dylib'
new_libssl = os.path.join(sage_local, 'lib', 'libssl.1.1.dylib')
old_libcrypto = '/usr/local/lib/libcrypto.1.1.dylib'
new_libcrypto = os.path.join(sage_local, 'lib', 'libcrypto.1.1.dylib')
ssl_lib = os.path.join(sage_libdynload, '_ssl.so')

from . import __path__
tk_path = '/Library/Frameworks/Tk.framework/Versions/8.6'

def tk_installed():
    if os.path.exists(tk_path):
        tkconfig = os.path.join(tk_path, 'tkConfig.sh')
        output = subprocess.run(['grep', 'TK_PATCH_LEVEL', tkconfig],
                                   capture_output=True).stdout
        tk_patch = int(output.split(b"'")[1][1:])
        if tk_patch != 10:
            return False
    else:
        return False
    return True

if __name__ == '__main__':
    here = __path__[0]
    print("Installing Tk 8.6.10 in /Library/Frameworks -- please choose the defaults.")
    subprocess.run(['open', '-W', os.path.join(here, 'TclTk-8.6.10.pkg')])
    print("Adding the _tkinter and _ssl modules to Sage")
    subprocess.run(['cp', '-r', os.path.join(here, 'local/'), sage_local]) 
    print('Rewriting id and dependency paths in _ssl, libssl and libcrypto.')
    subprocess.run(['install_name_tool', '-id', ssl_lib, ssl_lib])
    subprocess.run(['install_name_tool', '-change', old_libssl, new_libssl, ssl_lib])
    subprocess.run(['install_name_tool', '-change', old_libcrypto, new_libcrypto, ssl_lib])
    subprocess.run(['install_name_tool', '-id', new_libssl, new_libssl])
    subprocess.run(['install_name_tool', '-change', old_libcrypto, new_libcrypto, new_libssl])
    subprocess.run(['install_name_tool', '-id', new_libcrypto, new_libcrypto])
