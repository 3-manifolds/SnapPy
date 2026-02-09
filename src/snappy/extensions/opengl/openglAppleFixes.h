
/* When Apple backported vertex array objects to OpenGL 2.1 and
   introduced glGenVertexArraysAPPLE, ...
   it also remvoed the OpenGL functions from the header.
   Note that glGenVertexArraysApple, ... give an invalid operation
   OpenGL error in a modern OpenGL context. */

extern void glGenVertexArrays( GLsizei n, GLuint *arrays );
extern void glBindVertexArray( GLuint array );
extern void glDeleteVertexArrays( GLsizei n, const GLuint *arrays );
