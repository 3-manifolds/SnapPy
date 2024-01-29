def clear_gl_errors():
    """
    Clears any previous OpenGL errors.
    """

    while glGetError() != GL_NO_ERROR:
        pass

def print_gl_errors(msg):
    """
    Prints all OpenGL errors using given message.
    """

    while True:
        err = glGetError()
        if err == GL_NO_ERROR:
            return
        if err == GL_INVALID_ENUM:
            k = "GL_INVALID_ENUM"
        elif err == GL_INVALID_VALUE:
            k = "GL_INVALID_VALUE"
        elif err == GL_INVALID_OPERATION:
            k = "GL_INVALID_OPERATION"
        elif err == GL_INVALID_FRAMEBUFFER_OPERATION:
            k = "GL_INVALID_FRAMEBUFFER_OPERATION"
        else:
            k = "GL_ENUM 0x%x" % err
        print("Error %s in %s" % (k, msg))
