cdef class DataBasedTexture:
    cdef GLuint _textureName

    def __cinit__(self):
        self._textureName = 0

    def __init__(self,
                 unsigned int width,
                 unsigned int height,
                 unsigned char[:] rgba_data):

        if rgba_data.size != 4 * width * height:
            raise RuntimeError("Length of rgba_data not matching")

        glGenTextures(1, &self._textureName)
        glActiveTexture(GL_TEXTURE0)
        glBindTexture(GL_TEXTURE_2D, self._textureName)

        glTexImage2D(GL_TEXTURE_2D, 0,
                     GL_RGBA,
                     width,
                     height,
                     0,
                     GL_RGBA,
                     GL_UNSIGNED_BYTE,
                     &rgba_data[0])

        glTexParameteri(GL_TEXTURE_2D,
                        GL_TEXTURE_MIN_FILTER,
                        GL_LINEAR)
        glTexParameteri(GL_TEXTURE_2D,
                        GL_TEXTURE_MAG_FILTER,
                        GL_LINEAR)

        glBindTexture(GL_TEXTURE_2D, 0)

    def bind(self):
        glBindTexture(GL_TEXTURE_2D, self._textureName)

    def unbind(self):
        glBindTexture(GL_TEXTURE_2D, 0)

    def delete_resource(self):
        glDeleteTextures(1, &self._textureName)
        self._textureName = 0
