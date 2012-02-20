Building releases for Linux
===================

One needs two VMs, one 32bit and one 64bit. Currently, these are
Ubuntu 11.10 (Ocelot); when updates for these cease (early 2013),
switch to 12.04 a LTS release good for five years.

VMWare setup:

* VMWare tools didn't install properly; tried installing package
"open-vm-tools" which didn't really help; next time perhaps download a
prebuild VM with it installed.  

* In VM->Settings->Network select NAT and, under "Advanced options"
have it generate a MAC Address.  

* Edit:: 

  /Library/Application Support/VMware Fusion/vmnet8/dhcpd.conf

adding blocks like::

       host ubuntu64 {
       	    hardware ethernet 00:50:56:2e:54:ab;
	    fixed-address 192.168.3.10;
	}
	host ubuntu32 {
     	     hardware ethernet 00:50:56:26:82:43;
	     fixed-address 192.168.3.11;
	}

where the first part "192.168.3." should match that earlier in the
file, and the suffix should be something in the range 10-100.  Then
do::

	sudo "/Library/Application Support/VMware Fusion/boot.sh"  --restart

For further details, see::

    http://www.thirdbit.net/articles/2008/03/04/dhcp-on-vmware-fusion/

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
   scp ~/work/SnapPy/building_releases/update_SnapPy_linux.py  192.168.3.*:bin/update_SnapPy.py

On VM::

   cd ~
   hg_t3m_clone SnapPy
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


   


