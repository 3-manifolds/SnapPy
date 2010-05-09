/*
This hack avoids errors when importing Cython modules which were compiled on
a system with a different size for unicode characters.  This will have to
be revisited when we switch to Python 3.
*/

#ifdef __linux__
#define PyUnicode_DecodeUTF8 UCS2_hack
#endif
