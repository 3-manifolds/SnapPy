cdef extern from "gl_init.h":
     char * initLegacyGLAndReturnError();
     char * initModernGLAndReturnError();
