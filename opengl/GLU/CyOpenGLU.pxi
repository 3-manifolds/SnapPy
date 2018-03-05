cdef extern from "glu.h":
    cdef enum:

#  Extensions
        GLU_EXT_object_space_tess
        GLU_EXT_nurbs_tessellator

#  Boolean
        GLU_FALSE
        GLU_TRUE

#  Version
        GLU_VERSION_1_1
        GLU_VERSION_1_2
        GLU_VERSION_1_3

#  StringName
        GLU_VERSION
        GLU_EXTENSIONS

#  ErrorCode
        GLU_INVALID_ENUM
        GLU_INVALID_VALUE
        GLU_OUT_OF_MEMORY
        GLU_INCOMPATIBLE_GL_VERSION
        GLU_INVALID_OPERATION

#  NurbsDisplay
        GLU_OUTLINE_POLYGON
        GLU_OUTLINE_PATCH

#  NurbsCallback
        GLU_NURBS_ERROR
        GLU_ERROR
        GLU_NURBS_BEGIN
        GLU_NURBS_BEGIN_EXT
        GLU_NURBS_VERTEX
        GLU_NURBS_VERTEX_EXT
        GLU_NURBS_NORMAL
        GLU_NURBS_NORMAL_EXT
        GLU_NURBS_COLOR
        GLU_NURBS_COLOR_EXT
        GLU_NURBS_TEXTURE_COORD
        GLU_NURBS_TEX_COORD_EXT
        GLU_NURBS_END
        GLU_NURBS_END_EXT
        GLU_NURBS_BEGIN_DATA
        GLU_NURBS_BEGIN_DATA_EXT
        GLU_NURBS_VERTEX_DATA
        GLU_NURBS_VERTEX_DATA_EXT
        GLU_NURBS_NORMAL_DATA
        GLU_NURBS_NORMAL_DATA_EXT
        GLU_NURBS_COLOR_DATA
        GLU_NURBS_COLOR_DATA_EXT
        GLU_NURBS_TEXTURE_COORD_DATA
        GLU_NURBS_TEX_COORD_DATA_EXT
        GLU_NURBS_END_DATA
        GLU_NURBS_END_DATA_EXT

#  NurbsError
        GLU_NURBS_ERROR1
        GLU_NURBS_ERROR2
        GLU_NURBS_ERROR3
        GLU_NURBS_ERROR4
        GLU_NURBS_ERROR5
        GLU_NURBS_ERROR6
        GLU_NURBS_ERROR7
        GLU_NURBS_ERROR8
        GLU_NURBS_ERROR9
        GLU_NURBS_ERROR10
        GLU_NURBS_ERROR11
        GLU_NURBS_ERROR12
        GLU_NURBS_ERROR13
        GLU_NURBS_ERROR14
        GLU_NURBS_ERROR15
        GLU_NURBS_ERROR16
        GLU_NURBS_ERROR17
        GLU_NURBS_ERROR18
        GLU_NURBS_ERROR19
        GLU_NURBS_ERROR20
        GLU_NURBS_ERROR21
        GLU_NURBS_ERROR22
        GLU_NURBS_ERROR23
        GLU_NURBS_ERROR24
        GLU_NURBS_ERROR25
        GLU_NURBS_ERROR26
        GLU_NURBS_ERROR27
        GLU_NURBS_ERROR28
        GLU_NURBS_ERROR29
        GLU_NURBS_ERROR30
        GLU_NURBS_ERROR31
        GLU_NURBS_ERROR32
        GLU_NURBS_ERROR33
        GLU_NURBS_ERROR34
        GLU_NURBS_ERROR35
        GLU_NURBS_ERROR36
        GLU_NURBS_ERROR37

#  NurbsProperty
        GLU_AUTO_LOAD_MATRIX
        GLU_CULLING
        GLU_SAMPLING_TOLERANCE
        GLU_DISPLAY_MODE
        GLU_PARAMETRIC_TOLERANCE
        GLU_SAMPLING_METHOD
        GLU_U_STEP
        GLU_V_STEP
        GLU_NURBS_MODE
        GLU_NURBS_MODE_EXT
        GLU_NURBS_TESSELLATOR
        GLU_NURBS_TESSELLATOR_EXT
        GLU_NURBS_RENDERER
        GLU_NURBS_RENDERER_EXT

#  NurbsSampling
        GLU_OBJECT_PARAMETRIC_ERROR
        GLU_OBJECT_PARAMETRIC_ERROR_EXT
        GLU_OBJECT_PATH_LENGTH
        GLU_OBJECT_PATH_LENGTH_EXT
        GLU_PATH_LENGTH
        GLU_PARAMETRIC_ERROR
        GLU_DOMAIN_DISTANCE

#  NurbsTrim
        GLU_MAP1_TRIM_2
        GLU_MAP1_TRIM_3

#  QuadricDrawStyle
        GLU_POINT
        GLU_LINE
        GLU_FILL
        GLU_SILHOUETTE

#  QuadricNormal
        GLU_SMOOTH
        GLU_FLAT
        GLU_NONE

#  QuadricOrientation
        GLU_OUTSIDE
        GLU_INSIDE

#  TessCallback
        GLU_TESS_BEGIN
        GLU_BEGIN
        GLU_TESS_VERTEX
        GLU_VERTEX
        GLU_TESS_END
        GLU_END
        GLU_TESS_ERROR
        GLU_TESS_EDGE_FLAG
        GLU_EDGE_FLAG
        GLU_TESS_COMBINE
        GLU_TESS_BEGIN_DATA
        GLU_TESS_VERTEX_DATA
        GLU_TESS_END_DATA
        GLU_TESS_ERROR_DATA
        GLU_TESS_EDGE_FLAG_DATA
        GLU_TESS_COMBINE_DATA

#  TessContour
        GLU_CW
        GLU_CCW
        GLU_INTERIOR
        GLU_EXTERIOR
        GLU_UNKNOWN

#  TessProperty
        GLU_TESS_WINDING_RULE
        GLU_TESS_BOUNDARY_ONLY
        GLU_TESS_TOLERANCE

#  TessError
        GLU_TESS_ERROR1
        GLU_TESS_ERROR2
        GLU_TESS_ERROR3
        GLU_TESS_ERROR4
        GLU_TESS_ERROR5
        GLU_TESS_ERROR6
        GLU_TESS_ERROR7
        GLU_TESS_ERROR8
        GLU_TESS_MISSING_BEGIN_POLYGON
        GLU_TESS_MISSING_BEGIN_CONTOUR
        GLU_TESS_MISSING_END_POLYGON
        GLU_TESS_MISSING_END_CONTOUR
        GLU_TESS_COORD_TOO_LARGE
        GLU_TESS_NEED_COMBINE_CALLBACK

#  TessWinding
        GLU_TESS_WINDING_ODD
        GLU_TESS_WINDING_NONZERO
        GLU_TESS_WINDING_POSITIVE
        GLU_TESS_WINDING_NEGATIVE
        GLU_TESS_WINDING_ABS_GEQ_TWO

    cdef struct GLUnurbs
    cdef struct GLUquadric
    cdef struct GLUtesselator

    ctypedef GLUnurbs GLUnurbsObj
    ctypedef GLUquadric GLUquadricObj
    ctypedef GLUtesselator GLUtesselatorObj
    ctypedef GLUtesselator GLUtriangulatorObj

    cdef float GLU_TESS_MAX_COORD = 1.0e150

# Internal convenience typedefs
#    ctypedef void (void* _GLUfuncptr)()

    cdef void gluBeginCurve ( GLUnurbs* nurb )
    cdef void gluBeginPolygon ( GLUtesselator* tess )
    cdef void gluBeginSurface ( GLUnurbs* nurb )
    cdef void gluBeginTrim ( GLUnurbs* nurb )
    cdef GLint gluBuild1DMipmapLevels ( GLenum target, GLint internalFormat,
                                        GLsizei width, GLenum format, GLenum type,
                                        GLint level, GLint base, GLint max,
                                        void *data )
    cdef GLint gluBuild1DMipmaps ( GLenum target, GLint internalFormat,
                                   GLsizei width, GLenum format, GLenum type,
                                   void *data )
    cdef GLint gluBuild2DMipmapLevels ( GLenum target, GLint internalFormat,
                                        GLsizei width, GLsizei height,
                                        GLenum format, GLenum type,
                                        GLint level, GLint base, GLint max,
                                        void *data )
    cdef GLint gluBuild2DMipmaps ( GLenum target, GLint internalFormat,
                                   GLsizei width, GLsizei height, GLenum format,
                                   GLenum type, void *data )
    cdef GLint gluBuild3DMipmapLevels ( GLenum target, GLint internalFormat,
                                        GLsizei width, GLsizei height,
                                        GLsizei depth, GLenum format,
                                        GLenum type, GLint level, GLint base,
                                        GLint max, void *data )
    cdef GLint gluBuild3DMipmaps ( GLenum target, GLint internalFormat,
                                   GLsizei width, GLsizei height, GLsizei depth,
                                   GLenum format, GLenum type, void *data )
    cdef GLboolean gluCheckExtension ( GLubyte *extName, GLubyte *extString )
    cdef void gluCylinder ( GLUquadric* quad, GLdouble base, GLdouble top,
                            GLdouble height, GLint slices, GLint stacks )
    cdef void gluDeleteNurbsRenderer ( GLUnurbs* nurb )
    cdef void gluDeleteQuadric ( GLUquadric* quad )
    cdef void gluDeleteTess ( GLUtesselator* tess )
    cdef void gluDisk ( GLUquadric* quad, GLdouble inner, GLdouble outer,
                        GLint slices, GLint loops )
    cdef void gluEndCurve ( GLUnurbs* nurb )
    cdef void gluEndPolygon ( GLUtesselator* tess )
    cdef void gluEndSurface ( GLUnurbs* nurb )
    cdef void gluEndTrim ( GLUnurbs* nurb )
    cdef GLubyte * gluErrorString ( GLenum error )
    cdef void gluGetNurbsProperty ( GLUnurbs* nurb, GLenum property, GLfloat* data )
    cdef GLubyte * gluGetString ( GLenum name )
    cdef void gluGetTessProperty ( GLUtesselator* tess, GLenum which,
                                   GLdouble* data )
    cdef void gluLoadSamplingMatrices ( GLUnurbs* nurb, GLfloat *model,
                                        GLfloat *perspective, GLint *view )
    cdef void gluLookAt ( GLdouble eyeX, GLdouble eyeY, GLdouble eyeZ,
                          GLdouble centerX, GLdouble centerY, GLdouble centerZ,
                          GLdouble upX, GLdouble upY, GLdouble upZ )
    cdef GLUnurbs* gluNewNurbsRenderer ( )
    cdef GLUquadric* gluNewQuadric ( )
    cdef GLUtesselator* gluNewTess ( )
    cdef void gluNextContour ( GLUtesselator* tess, GLenum type )
    cdef void gluNurbsCallback ( GLUnurbs* nurb, GLenum which,
                                 void (*CallbackFunc)() )
    cdef void gluNurbsCallbackData ( GLUnurbs* nurb, GLvoid* userData )
    cdef void gluNurbsCallbackDataEXT ( GLUnurbs* nurb, GLvoid* userData )
    cdef void gluNurbsCurve ( GLUnurbs* nurb, GLint knotCount, GLfloat *knots,
                              GLint stride, GLfloat *control, GLint order,
                              GLenum type )
    cdef void gluNurbsProperty ( GLUnurbs* nurb, GLenum property, GLfloat value )
    cdef void gluNurbsSurface ( GLUnurbs* nurb, GLint sKnotCount, GLfloat* sKnots,
                                GLint tKnotCount, GLfloat* tKnots, GLint sStride,
                                GLint tStride, GLfloat* control, GLint sOrder,
                                GLint tOrder, GLenum type )
    cdef void gluOrtho2D ( GLdouble left, GLdouble right, GLdouble bottom,
                           GLdouble top )
    cdef void gluPartialDisk ( GLUquadric* quad, GLdouble inner, GLdouble outer,
                               GLint slices, GLint loops, GLdouble start,
                               GLdouble sweep )
    cdef void gluPerspective ( GLdouble fovy, GLdouble aspect, GLdouble zNear,
                               GLdouble zFar )
    cdef void gluPickMatrix ( GLdouble x, GLdouble y, GLdouble delX, GLdouble delY,
                              GLint *viewport )
    cdef GLint gluProject ( GLdouble objX, GLdouble objY, GLdouble objZ,
                            GLdouble *model, GLdouble *proj, GLint *view,
                            GLdouble* winX, GLdouble* winY, GLdouble* winZ )
    cdef void gluPwlCurve ( GLUnurbs* nurb, GLint count, GLfloat* data,
                            GLint stride, GLenum type )
    cdef void gluQuadricCallback ( GLUquadric* quad, GLenum which,
                                   void (*CallbackFunc)() )
    cdef void gluQuadricDrawStyle ( GLUquadric* quad, GLenum draw )
    cdef void gluQuadricNormals ( GLUquadric* quad, GLenum normal )
    cdef void gluQuadricOrientation ( GLUquadric* quad, GLenum orientation )
    cdef void gluQuadricTexture ( GLUquadric* quad, GLboolean texture )
    cdef GLint gluScaleImage ( GLenum format, GLsizei wIn, GLsizei hIn,
                               GLenum typeIn, void *dataIn, GLsizei wOut,
                               GLsizei hOut, GLenum typeOut, GLvoid* dataOut )
    cdef void gluSphere ( GLUquadric* quad, GLdouble radius, GLint slices,
                          GLint stacks )
    cdef void gluTessBeginContour ( GLUtesselator* tess )
    cdef void gluTessBeginPolygon ( GLUtesselator* tess, GLvoid* data )
    cdef void gluTessCallback ( GLUtesselator* tess, GLenum which,
                                void (*CallbackFunc)() )
    cdef void gluTessEndContour ( GLUtesselator* tess )
    cdef void gluTessEndPolygon ( GLUtesselator* tess )
    cdef void gluTessNormal ( GLUtesselator* tess, GLdouble valueX,
                              GLdouble valueY, GLdouble valueZ )
    cdef void gluTessProperty ( GLUtesselator* tess, GLenum which, GLdouble data )
    cdef void gluTessVertex ( GLUtesselator* tess, GLdouble *location, GLvoid* data )
    cdef GLint gluUnProject ( GLdouble winX, GLdouble winY, GLdouble winZ,
                              GLdouble *model, GLdouble *proj, GLint *view,
                              GLdouble* objX, GLdouble* objY, GLdouble* objZ )
    cdef GLint gluUnProject4 ( GLdouble winX, GLdouble winY, GLdouble winZ,
                               GLdouble clipW, GLdouble *model, GLdouble *proj,
                               GLint *view, GLdouble nearVal, GLdouble farVal,
                               GLdouble* objX, GLdouble* objY, GLdouble* objZ,
                               GLdouble* objW )
