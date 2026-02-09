#ifndef CYOPENGL_GL_INIT
#define CYOPENGL_GL_INIT

/*
 *
 * Methods to load GL functions:
 * initLegacyGLAndReturnError and initModernGLAndReturnError.
 *
 */


#define CHECK_FOR_GL_FUNCTION_PTR(name) if (!name) { return "Missing function " #name; }

char * initLegacyGLAndReturnError()
{
#ifdef GL_HEADERS_FROM_GLEW
    GLenum err = glewInit();
    if (err != 0) {
        return "glewInit failed";
    }

    CHECK_FOR_GL_FUNCTION_PTR(glMatrixMode);
    CHECK_FOR_GL_FUNCTION_PTR(glLoadIdentity);
    CHECK_FOR_GL_FUNCTION_PTR(glOrtho);
    CHECK_FOR_GL_FUNCTION_PTR(glFrustum);
    CHECK_FOR_GL_FUNCTION_PTR(glPushMatrix);
    CHECK_FOR_GL_FUNCTION_PTR(glPopMatrix);
    CHECK_FOR_GL_FUNCTION_PTR(glTranslatef);
    CHECK_FOR_GL_FUNCTION_PTR(glRotatef);
    CHECK_FOR_GL_FUNCTION_PTR(glMultMatrixd);
    CHECK_FOR_GL_FUNCTION_PTR(glBegin);
    CHECK_FOR_GL_FUNCTION_PTR(glEnd);
    CHECK_FOR_GL_FUNCTION_PTR(glNormal3f);
    CHECK_FOR_GL_FUNCTION_PTR(glVertex3f);
    CHECK_FOR_GL_FUNCTION_PTR(glColor4fv);
    CHECK_FOR_GL_FUNCTION_PTR(glNewList);
    CHECK_FOR_GL_FUNCTION_PTR(glEndList);
    CHECK_FOR_GL_FUNCTION_PTR(glIsList);
    CHECK_FOR_GL_FUNCTION_PTR(glCallList);
    CHECK_FOR_GL_FUNCTION_PTR(glDeleteLists);
    CHECK_FOR_GL_FUNCTION_PTR(glMaterialf);
    CHECK_FOR_GL_FUNCTION_PTR(glMaterialfv);
    CHECK_FOR_GL_FUNCTION_PTR(glColorMaterial);
    CHECK_FOR_GL_FUNCTION_PTR(glShadeModel);
    CHECK_FOR_GL_FUNCTION_PTR(glLightModeli);
    CHECK_FOR_GL_FUNCTION_PTR(glLightModelf);
    CHECK_FOR_GL_FUNCTION_PTR(glLightModelfv);
    CHECK_FOR_GL_FUNCTION_PTR(glLightf);
    CHECK_FOR_GL_FUNCTION_PTR(glLightfv);
    CHECK_FOR_GL_FUNCTION_PTR(glLineWidth);
    CHECK_FOR_GL_FUNCTION_PTR(glLineStipple);
    CHECK_FOR_GL_FUNCTION_PTR(glRasterPos2f);
    CHECK_FOR_GL_FUNCTION_PTR(glRasterPos3f);
    CHECK_FOR_GL_FUNCTION_PTR(glBitmap);
    CHECK_FOR_GL_FUNCTION_PTR(glDrawPixels);
    CHECK_FOR_GL_FUNCTION_PTR(glGetFloatv);
    CHECK_FOR_GL_FUNCTION_PTR(glGetDoublev);
    CHECK_FOR_GL_FUNCTION_PTR(glBlendFunc);
    CHECK_FOR_GL_FUNCTION_PTR(glFrontFace);

    CHECK_FOR_GL_FUNCTION_PTR(glNormalPointer);
    CHECK_FOR_GL_FUNCTION_PTR(glVertexPointer);
    CHECK_FOR_GL_FUNCTION_PTR(glEnableClientState);
    CHECK_FOR_GL_FUNCTION_PTR(glDrawElements);

    CHECK_FOR_GL_FUNCTION_PTR(glViewport);
    CHECK_FOR_GL_FUNCTION_PTR(glClear);
    CHECK_FOR_GL_FUNCTION_PTR(glClearColor);
    CHECK_FOR_GL_FUNCTION_PTR(glGetString);
    CHECK_FOR_GL_FUNCTION_PTR(glGetError);
    CHECK_FOR_GL_FUNCTION_PTR(glEnable);
    CHECK_FOR_GL_FUNCTION_PTR(glDisable);
    #endif

        return NULL;
    }

char * initModernGLAndReturnError()
{
#ifdef GL_HEADERS_FROM_GLEW
    GLenum err = glewInit();
    if (err != 0) {
        return "glewInit failed";
    }

    CHECK_FOR_GL_FUNCTION_PTR(glCreateShader);
    CHECK_FOR_GL_FUNCTION_PTR(glAttachShader);
    CHECK_FOR_GL_FUNCTION_PTR(glShaderSource);
    CHECK_FOR_GL_FUNCTION_PTR(glCompileShader);
    CHECK_FOR_GL_FUNCTION_PTR(glGetProgramiv);
    CHECK_FOR_GL_FUNCTION_PTR(glGetShaderiv);
    CHECK_FOR_GL_FUNCTION_PTR(glGetProgramInfoLog);
    CHECK_FOR_GL_FUNCTION_PTR(glGetShaderInfoLog);
    CHECK_FOR_GL_FUNCTION_PTR(glLinkProgram);
    CHECK_FOR_GL_FUNCTION_PTR(glUseProgram);
    CHECK_FOR_GL_FUNCTION_PTR(glDeleteShader);
    CHECK_FOR_GL_FUNCTION_PTR(glDeleteProgram);

    CHECK_FOR_GL_FUNCTION_PTR(glGetUniformLocation);
    CHECK_FOR_GL_FUNCTION_PTR(glUniform1i);
    CHECK_FOR_GL_FUNCTION_PTR(glUniform1f);
    CHECK_FOR_GL_FUNCTION_PTR(glUniform2i);
    CHECK_FOR_GL_FUNCTION_PTR(glUniform2f);
    CHECK_FOR_GL_FUNCTION_PTR(glUniform3i);
    CHECK_FOR_GL_FUNCTION_PTR(glUniform3f);
    CHECK_FOR_GL_FUNCTION_PTR(glUniform4i);
    CHECK_FOR_GL_FUNCTION_PTR(glUniform4f);
    CHECK_FOR_GL_FUNCTION_PTR(glUniform1iv);
    CHECK_FOR_GL_FUNCTION_PTR(glUniform1fv);
    CHECK_FOR_GL_FUNCTION_PTR(glUniform2iv);
    CHECK_FOR_GL_FUNCTION_PTR(glUniform2fv);
    CHECK_FOR_GL_FUNCTION_PTR(glUniform3iv);
    CHECK_FOR_GL_FUNCTION_PTR(glUniform3fv);
    CHECK_FOR_GL_FUNCTION_PTR(glUniform4iv);
    CHECK_FOR_GL_FUNCTION_PTR(glUniform4fv);
    CHECK_FOR_GL_FUNCTION_PTR(glUniformMatrix2fv);
    CHECK_FOR_GL_FUNCTION_PTR(glUniformMatrix4fv);
    CHECK_FOR_GL_FUNCTION_PTR(glUniformMatrix2x3fv);
    CHECK_FOR_GL_FUNCTION_PTR(glUniformMatrix3x2fv);

    CHECK_FOR_GL_FUNCTION_PTR(glGetUniformBlockIndex);
    CHECK_FOR_GL_FUNCTION_PTR(glUniformBlockBinding);
    CHECK_FOR_GL_FUNCTION_PTR(glBindBufferBase);

    CHECK_FOR_GL_FUNCTION_PTR(glGenVertexArrays);
    CHECK_FOR_GL_FUNCTION_PTR(glBindVertexArray);
    CHECK_FOR_GL_FUNCTION_PTR(glGenBuffers);
    CHECK_FOR_GL_FUNCTION_PTR(glBindBuffer);
    CHECK_FOR_GL_FUNCTION_PTR(glBufferData);
    CHECK_FOR_GL_FUNCTION_PTR(glVertexAttribPointer);
    CHECK_FOR_GL_FUNCTION_PTR(glDeleteBuffers);
    CHECK_FOR_GL_FUNCTION_PTR(glDeleteVertexArrays);
    CHECK_FOR_GL_FUNCTION_PTR(glDrawArrays);

    CHECK_FOR_GL_FUNCTION_PTR(glViewport);
    CHECK_FOR_GL_FUNCTION_PTR(glClear);
    CHECK_FOR_GL_FUNCTION_PTR(glClearColor);
    CHECK_FOR_GL_FUNCTION_PTR(glGetString);
    CHECK_FOR_GL_FUNCTION_PTR(glGetError);
    CHECK_FOR_GL_FUNCTION_PTR(glEnable);
    CHECK_FOR_GL_FUNCTION_PTR(glDisable);
    CHECK_FOR_GL_FUNCTION_PTR(glReadPixels);

    CHECK_FOR_GL_FUNCTION_PTR(glGenTextures);
    CHECK_FOR_GL_FUNCTION_PTR(glBindTexture);
    CHECK_FOR_GL_FUNCTION_PTR(glActiveTexture);
    CHECK_FOR_GL_FUNCTION_PTR(glDeleteTextures);

    CHECK_FOR_GL_FUNCTION_PTR(glGenFramebuffers);
    CHECK_FOR_GL_FUNCTION_PTR(glBindFramebuffer);
    CHECK_FOR_GL_FUNCTION_PTR(glFramebufferTexture2D);
    CHECK_FOR_GL_FUNCTION_PTR(glCheckFramebufferStatus);
    CHECK_FOR_GL_FUNCTION_PTR(glDeleteFramebuffers);

#endif
    return NULL;
}

#endif
