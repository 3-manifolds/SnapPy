
# Important:
# The C code below must occur after including the opengl.h header!
# Otherwise, we don't pick up the USE_GLEW macro in the opengl.h header.

cdef extern from *:
    """
    #define CHECK_GLEW_FOR_FUNCTION(name) if (!name) { return #name; }

    /* Check that GLEW set the gl function pointers */
    char * checkGlewForLegacyOpenGL()
    {
    #ifdef USE_GLEW
        GLenum err = glewInit();
        if (err != 0) {
            return "glewInit failed";
        }

        CHECK_GLEW_FOR_FUNCTION(glMatrixMode);
        CHECK_GLEW_FOR_FUNCTION(glLoadIdentity);
        CHECK_GLEW_FOR_FUNCTION(glOrtho);
        CHECK_GLEW_FOR_FUNCTION(glFrustum);
        CHECK_GLEW_FOR_FUNCTION(glPushMatrix);
        CHECK_GLEW_FOR_FUNCTION(glPopMatrix);
        CHECK_GLEW_FOR_FUNCTION(glTranslatef);
        CHECK_GLEW_FOR_FUNCTION(glRotatef);
        CHECK_GLEW_FOR_FUNCTION(glMultMatrixd);
        CHECK_GLEW_FOR_FUNCTION(glBegin);
        CHECK_GLEW_FOR_FUNCTION(glEnd);
        CHECK_GLEW_FOR_FUNCTION(glNormal3f);
        CHECK_GLEW_FOR_FUNCTION(glVertex3f);
        CHECK_GLEW_FOR_FUNCTION(glColor4fv);
        CHECK_GLEW_FOR_FUNCTION(glNewList);
        CHECK_GLEW_FOR_FUNCTION(glEndList);
        CHECK_GLEW_FOR_FUNCTION(glIsList);
        CHECK_GLEW_FOR_FUNCTION(glCallList);
        CHECK_GLEW_FOR_FUNCTION(glDeleteLists);
        CHECK_GLEW_FOR_FUNCTION(glMaterialf);
        CHECK_GLEW_FOR_FUNCTION(glMaterialfv);
        CHECK_GLEW_FOR_FUNCTION(glColorMaterial);
        CHECK_GLEW_FOR_FUNCTION(glShadeModel);
        CHECK_GLEW_FOR_FUNCTION(glLightModeli);
        CHECK_GLEW_FOR_FUNCTION(glLightModelf);
        CHECK_GLEW_FOR_FUNCTION(glLightModelfv);
        CHECK_GLEW_FOR_FUNCTION(glLightf);
        CHECK_GLEW_FOR_FUNCTION(glLightfv);
        CHECK_GLEW_FOR_FUNCTION(glLineWidth);
        CHECK_GLEW_FOR_FUNCTION(glLineStipple);
        CHECK_GLEW_FOR_FUNCTION(glRasterPos2f);
        CHECK_GLEW_FOR_FUNCTION(glRasterPos3f);
        CHECK_GLEW_FOR_FUNCTION(glBitmap);
        CHECK_GLEW_FOR_FUNCTION(glDrawPixels);
        CHECK_GLEW_FOR_FUNCTION(glGetFloatv);
        CHECK_GLEW_FOR_FUNCTION(glGetDoublev);
        CHECK_GLEW_FOR_FUNCTION(glBlendFunc);
        CHECK_GLEW_FOR_FUNCTION(glFrontFace);
        
        CHECK_GLEW_FOR_FUNCTION(glNormalPointer);
        CHECK_GLEW_FOR_FUNCTION(glVertexPointer);
        CHECK_GLEW_FOR_FUNCTION(glEnableClientState);
        CHECK_GLEW_FOR_FUNCTION(glDrawElements);

        CHECK_GLEW_FOR_FUNCTION(glViewport);
        CHECK_GLEW_FOR_FUNCTION(glClear);
        CHECK_GLEW_FOR_FUNCTION(glClearColor);
        CHECK_GLEW_FOR_FUNCTION(glGetString);
        CHECK_GLEW_FOR_FUNCTION(glGetError);
        CHECK_GLEW_FOR_FUNCTION(glEnable);
        CHECK_GLEW_FOR_FUNCTION(glDisable);
    #endif

        return NULL;
    }

    /* Check that GLEW set the gl function pointers */
    char * checkGlewForModernOpenGL()
    {
    #ifdef USE_GLEW
        GLenum err = glewInit();
        if (err != 0) {
            return "glewInit failed";
        }

        CHECK_GLEW_FOR_FUNCTION(glCreateShader);
        CHECK_GLEW_FOR_FUNCTION(glAttachShader);
        CHECK_GLEW_FOR_FUNCTION(glShaderSource);
        CHECK_GLEW_FOR_FUNCTION(glCompileShader);
        CHECK_GLEW_FOR_FUNCTION(glGetProgramiv);
        CHECK_GLEW_FOR_FUNCTION(glGetShaderiv);
        CHECK_GLEW_FOR_FUNCTION(glGetProgramInfoLog);
        CHECK_GLEW_FOR_FUNCTION(glGetShaderInfoLog);
        CHECK_GLEW_FOR_FUNCTION(glLinkProgram);
        CHECK_GLEW_FOR_FUNCTION(glUseProgram);
        CHECK_GLEW_FOR_FUNCTION(glDeleteShader);
        CHECK_GLEW_FOR_FUNCTION(glDeleteProgram);

        CHECK_GLEW_FOR_FUNCTION(glGetUniformLocation);
        CHECK_GLEW_FOR_FUNCTION(glUniform1i);
        CHECK_GLEW_FOR_FUNCTION(glUniform1f);
        CHECK_GLEW_FOR_FUNCTION(glUniform2i);
        CHECK_GLEW_FOR_FUNCTION(glUniform2f);
        CHECK_GLEW_FOR_FUNCTION(glUniform3i);
        CHECK_GLEW_FOR_FUNCTION(glUniform3f);
        CHECK_GLEW_FOR_FUNCTION(glUniform4i);
        CHECK_GLEW_FOR_FUNCTION(glUniform4f);
        CHECK_GLEW_FOR_FUNCTION(glUniform1iv);
        CHECK_GLEW_FOR_FUNCTION(glUniform1fv);
        CHECK_GLEW_FOR_FUNCTION(glUniform2iv);
        CHECK_GLEW_FOR_FUNCTION(glUniform2fv);
        CHECK_GLEW_FOR_FUNCTION(glUniform3iv);
        CHECK_GLEW_FOR_FUNCTION(glUniform3fv);
        CHECK_GLEW_FOR_FUNCTION(glUniform4iv);
        CHECK_GLEW_FOR_FUNCTION(glUniform4fv);
        CHECK_GLEW_FOR_FUNCTION(glUniformMatrix2fv);
        CHECK_GLEW_FOR_FUNCTION(glUniformMatrix4fv);
        CHECK_GLEW_FOR_FUNCTION(glUniformMatrix2x3fv);
        CHECK_GLEW_FOR_FUNCTION(glUniformMatrix3x2fv);
        
        CHECK_GLEW_FOR_FUNCTION(glGetUniformBlockIndex);
        CHECK_GLEW_FOR_FUNCTION(glUniformBlockBinding);
        CHECK_GLEW_FOR_FUNCTION(glBindBufferBase);
        
        CHECK_GLEW_FOR_FUNCTION(glGenVertexArrays);
        CHECK_GLEW_FOR_FUNCTION(glBindVertexArray);
        CHECK_GLEW_FOR_FUNCTION(glGenBuffers);
        CHECK_GLEW_FOR_FUNCTION(glBindBuffer);
        CHECK_GLEW_FOR_FUNCTION(glBufferData);
        CHECK_GLEW_FOR_FUNCTION(glVertexAttribPointer);
        CHECK_GLEW_FOR_FUNCTION(glDeleteBuffers);
        CHECK_GLEW_FOR_FUNCTION(glDeleteVertexArrays);
        CHECK_GLEW_FOR_FUNCTION(glDrawArrays);

        CHECK_GLEW_FOR_FUNCTION(glViewport);
        CHECK_GLEW_FOR_FUNCTION(glClear);
        CHECK_GLEW_FOR_FUNCTION(glClearColor);
        CHECK_GLEW_FOR_FUNCTION(glGetString);
        CHECK_GLEW_FOR_FUNCTION(glGetError);
        CHECK_GLEW_FOR_FUNCTION(glEnable);
        CHECK_GLEW_FOR_FUNCTION(glDisable);
        CHECK_GLEW_FOR_FUNCTION(glReadPixels);

        CHECK_GLEW_FOR_FUNCTION(glGenTextures);
        CHECK_GLEW_FOR_FUNCTION(glBindTexture);
        CHECK_GLEW_FOR_FUNCTION(glActiveTexture);
        CHECK_GLEW_FOR_FUNCTION(glDeleteTextures);

        CHECK_GLEW_FOR_FUNCTION(glGenFramebuffers);
        CHECK_GLEW_FOR_FUNCTION(glBindFramebuffer);
        CHECK_GLEW_FOR_FUNCTION(glFramebufferTexture2D);
        CHECK_GLEW_FOR_FUNCTION(glCheckFramebufferStatus);
        CHECK_GLEW_FOR_FUNCTION(glDeleteFramebuffers);

    #endif
        return NULL;
    }

    """

    int callGlewInitIfNecessary();
    char * checkGlewForLegacyOpenGL();
    char * checkGlewForModernOpenGL();
