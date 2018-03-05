Building releases for Linux
===================

One needs two VMs, one 32bit and one 64bit. Currently, these are
Ubuntu 12.04 LTS (Penguin).  

VMWare Fusion setup: See info in sysadmin folder. 

Now, in the VM do::

    sudo hostname name-of-machine
    sudo apt-get install python-dev python2.6-dev python3.2-dev
    python-setuptools python-tk python3.2-tk g++ emacs libglu1-mesa-dev mercurial openssh-server


On VM::

   cd ~
   mkdir .ssh bin

and add the following to ".bashrc"::

    PATH="$HOME/bin:$PATH"
    umask 2

On host::

   scp ~/.ssh/authorized_keys2 192.168.3.*:.ssh/
   scp ~/bin/hg_t3m_clone 192.168.3.*:bin/

On VM::

   cd ~
      hg_t3m_clone SnapPy
   cd bin; ln -s ../SnapPy/release_tools/update_SnapPy_linux.py update_SnapPy.py; cd ../
   wget -nd https://raw.github.com/pypa/virtualenv/master/virtualenv.py
   python2.6 virtualenv.py python26
   python2.7 virtualenv.py python27
   python3.2 virtualenv.py python32
   python26/bin/python -m easy_install cython
   python27/bin/python -m easy_install cython
   python32/bin/python -m easy_install cython
   python26/bin/python -m easy_install sphinx
   python27/bin/python -m easy_install sphinx
   python32/bin/python -m easy_install sphinx
   cd SnapPy
   sh build_pari.sh

Now install PyPNG manually for Python 3.2, using the included
"pypng-0.0.12-p1.tar.gz" with a patched "setup.py".  Then, in the VM do::
   
   update_SnapPy.py


   


