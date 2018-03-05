SnapPea 2.0 File Formats

Strictly speaking, the SnapPea kernel is platform independent and
defines no file formats.  Nevertheless, while writing the Macintosh
SnapPea UI I tried to design file formats which would be suitable
on all platforms (Macintosh, Windows & Unix).  In particular, all
SnapPea files are standard text files.  I strongly encourage the use
of these file formats on all platforms, and can provide unix-style code
(fprintf/fscanf) for reading them and passing them to the SnapPea kernel.

SnapPea uses three types of files:

     triangulation files
     matrix generator files
     link projection files

Each is described in a separate document in this directory.

The only format the user needs to be aware of is that for matrix
generators, which is essentially just a list of matrices.
