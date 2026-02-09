cdef extern from "gl_init.h":
     cdef char * initLegacyGLAndReturnError( )
     cdef char * initModernGLAndReturnError( )
