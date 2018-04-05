cdef extern from "stdlib.h":

#    ctypedef unsigned long size_t
    void* malloc(size_t size)
    void free(void *mem)

cdef extern from "opengl.h":

# Datatypes
    ctypedef unsigned int GLenum
    ctypedef unsigned char GLboolean
    ctypedef unsigned int GLbitfield
    ctypedef void GLvoid
    ctypedef signed char GLbyte
    ctypedef short GLshort
    ctypedef int GLint
    ctypedef unsigned char GLubyte
    ctypedef unsigned short GLushort
    ctypedef unsigned int GLuint
    ctypedef int GLsizei
    ctypedef float GLfloat
    ctypedef float GLclampf
    ctypedef double GLdouble
    ctypedef double GLclampd

    ctypedef void* GLintptr
    ctypedef void* GLsizeiptr

# Constants
    cdef enum:
# Boolean values
        GL_FALSE
        GL_TRUE

# Data types
        GL_BYTE
        GL_UNSIGNED_BYTE
        GL_SHORT
        GL_UNSIGNED_SHORT
        GL_INT
        GL_UNSIGNED_INT
        GL_FLOAT
        GL_2_BYTES
        GL_3_BYTES
        GL_4_BYTES
        GL_DOUBLE

# Primitives

        GL_POINTS
        GL_LINES
        GL_LINE_LOOP
        GL_LINE_STRIP
        GL_TRIANGLES
        GL_TRIANGLE_STRIP
        GL_TRIANGLE_FAN
        GL_QUADS
        GL_QUAD_STRIP
        GL_POLYGON

# Vertex Arrays
        GL_VERTEX_ARRAY
        GL_NORMAL_ARRAY
        GL_COLOR_ARRAY
        GL_INDEX_ARRAY
        GL_TEXTURE_COORD_ARRAY
        GL_EDGE_FLAG_ARRAY
        GL_VERTEX_ARRAY_SIZE
        GL_VERTEX_ARRAY_TYPE
        GL_VERTEX_ARRAY_STRIDE
        GL_NORMAL_ARRAY_TYPE
        GL_NORMAL_ARRAY_STRIDE
        GL_COLOR_ARRAY_SIZE
        GL_COLOR_ARRAY_TYPE
        GL_COLOR_ARRAY_STRIDE
        GL_INDEX_ARRAY_TYPE
        GL_INDEX_ARRAY_STRIDE
        GL_TEXTURE_COORD_ARRAY_SIZE
        GL_TEXTURE_COORD_ARRAY_TYPE
        GL_TEXTURE_COORD_ARRAY_STRIDE
        GL_EDGE_FLAG_ARRAY_STRIDE
        GL_VERTEX_ARRAY_POINTER
        GL_NORMAL_ARRAY_POINTER
        GL_COLOR_ARRAY_POINTER
        GL_INDEX_ARRAY_POINTER
        GL_TEXTURE_COORD_ARRAY_POINTER
        GL_EDGE_FLAG_ARRAY_POINTER
        GL_V2F
        GL_V3F
        GL_C4UB_V2F
        GL_C4UB_V3F
        GL_C3F_V3F
        GL_N3F_V3F
        GL_C4F_N3F_V3F
        GL_T2F_V3F
        GL_T4F_V4F
        GL_T2F_C4UB_V3F
        GL_T2F_C3F_V3F
        GL_T2F_N3F_V3F
        GL_T2F_C4F_N3F_V3F
        GL_T4F_C4F_N3F_V4F

# Matrix Mode
        GL_MATRIX_MODE
        GL_MODELVIEW
        GL_PROJECTION
        GL_TEXTURE

# Points
        GL_POINT_SMOOTH
        GL_POINT_SIZE
        GL_POINT_SIZE_GRANULARITY 
        GL_POINT_SIZE_RANGE

# Lines
        GL_LINE_SMOOTH
        GL_LINE_STIPPLE
        GL_LINE_STIPPLE_PATTERN
        GL_LINE_STIPPLE_REPEAT
        GL_LINE_WIDTH
        GL_LINE_WIDTH_GRANULARITY
        GL_LINE_WIDTH_RANGE

# Polygons
        GL_POINT
        GL_LINE
        GL_FILL
        GL_CW
        GL_CCW
        GL_FRONT
        GL_BACK
        GL_POLYGON_MODE
        GL_POLYGON_SMOOTH
        GL_POLYGON_STIPPLE
        GL_EDGE_FLAG
        GL_CULL_FACE
        GL_CULL_FACE_MODE
        GL_FRONT_FACE
        GL_POLYGON_OFFSET_FACTOR
        GL_POLYGON_OFFSET_UNITS
        GL_POLYGON_OFFSET_POINT
        GL_POLYGON_OFFSET_LINE
        GL_POLYGON_OFFSET_FILL

# Display Lists
        GL_COMPILE
        GL_COMPILE_AND_EXECUTE
        GL_LIST_BASE
        GL_LIST_INDEX
        GL_LIST_MODE

# Depth buffer
        GL_NEVER
        GL_LESS
        GL_EQUAL
        GL_LEQUAL
        GL_GREATER
        GL_NOTEQUAL
        GL_GEQUAL
        GL_ALWAYS
        GL_DEPTH_TEST
        GL_DEPTH_BITS
        GL_DEPTH_CLEAR_VALUE
        GL_DEPTH_FUNC
        GL_DEPTH_RANGE
        GL_DEPTH_WRITEMASK
        GL_DEPTH_COMPONENT

# Lighting
        GL_LIGHTING
        GL_LIGHT0
        GL_LIGHT1
        GL_LIGHT2
        GL_LIGHT3
        GL_LIGHT4
        GL_LIGHT5
        GL_LIGHT6
        GL_LIGHT7
        GL_SPOT_EXPONENT
        GL_SPOT_CUTOFF
        GL_CONSTANT_ATTENUATION
        GL_LINEAR_ATTENUATION
        GL_QUADRATIC_ATTENUATION
        GL_AMBIENT
        GL_DIFFUSE
        GL_SPECULAR
        GL_SHININESS
        GL_EMISSION
        GL_POSITION
        GL_SPOT_DIRECTION
        GL_AMBIENT_AND_DIFFUSE
        GL_COLOR_INDEXES
        GL_LIGHT_MODEL_TWO_SIDE
        GL_LIGHT_MODEL_LOCAL_VIEWER
        GL_LIGHT_MODEL_AMBIENT
        GL_FRONT_AND_BACK
        GL_SHADE_MODEL
        GL_FLAT
        GL_SMOOTH
        GL_COLOR_MATERIAL
        GL_COLOR_MATERIAL_FACE
        GL_COLOR_MATERIAL_PARAMETER
        GL_NORMALIZE

# User clipping planes
        GL_CLIP_PLANE0
        GL_CLIP_PLANE1
        GL_CLIP_PLANE2
        GL_CLIP_PLANE3
        GL_CLIP_PLANE4
        GL_CLIP_PLANE5

# Accumulation buffer
        GL_ACCUM_RED_BITS
        GL_ACCUM_GREEN_BITS
        GL_ACCUM_BLUE_BITS
        GL_ACCUM_ALPHA_BITS
        GL_ACCUM_CLEAR_VALUE
        GL_ACCUM
        GL_ADD
        GL_LOAD
        GL_MULT
        GL_RETURN

# Alpha testing
        GL_ALPHA_TEST
        GL_ALPHA_TEST_REF
        GL_ALPHA_TEST_FUNC

# Blending
        GL_BLEND
        GL_BLEND_SRC
        GL_BLEND_DST
        GL_ZERO
        GL_ONE
        GL_SRC_COLOR
        GL_ONE_MINUS_SRC_COLOR
        GL_SRC_ALPHA
        GL_ONE_MINUS_SRC_ALPHA
        GL_DST_ALPHA
        GL_ONE_MINUS_DST_ALPHA
        GL_DST_COLOR
        GL_ONE_MINUS_DST_COLOR
        GL_SRC_ALPHA_SATURATE

# Render Mode
        GL_FEEDBACK
        GL_RENDER
        GL_SELECT

# Feedback
        GL_2D
        GL_3D
        GL_3D_COLOR
        GL_3D_COLOR_TEXTURE
        GL_4D_COLOR_TEXTURE
        GL_POINT_TOKEN
        GL_LINE_TOKEN
        GL_LINE_RESET_TOKEN
        GL_POLYGON_TOKEN
        GL_BITMAP_TOKEN
        GL_DRAW_PIXEL_TOKEN
        GL_COPY_PIXEL_TOKEN
        GL_PASS_THROUGH_TOKEN
        GL_FEEDBACK_BUFFER_POINTER
        GL_FEEDBACK_BUFFER_SIZE
        GL_FEEDBACK_BUFFER_TYPE

# Selection
        GL_SELECTION_BUFFER_POINTER
        GL_SELECTION_BUFFER_SIZE

# Fog
        GL_FOG
        GL_FOG_MODE
        GL_FOG_DENSITY
        GL_FOG_COLOR
        GL_FOG_INDEX
        GL_FOG_START
        GL_FOG_END
        GL_LINEAR
        GL_EXP
        GL_EXP2

# Logic Ops
        GL_LOGIC_OP
        GL_INDEX_LOGIC_OP
        GL_COLOR_LOGIC_OP
        GL_LOGIC_OP_MODE
        GL_CLEAR
        GL_SET
        GL_COPY
        GL_COPY_INVERTED
        GL_NOOP
        GL_INVERT
        GL_AND
        GL_NAND
        GL_OR
        GL_NOR
        GL_XOR
        GL_EQUIV
        GL_AND_REVERSE
        GL_AND_INVERTED
        GL_OR_REVERSE
        GL_OR_INVERTED

# Stencil
        GL_STENCIL_BITS
        GL_STENCIL_TEST
        GL_STENCIL_CLEAR_VALUE
        GL_STENCIL_FUNC
        GL_STENCIL_VALUE_MASK
        GL_STENCIL_FAIL
        GL_STENCIL_PASS_DEPTH_FAIL
        GL_STENCIL_PASS_DEPTH_PASS
        GL_STENCIL_REF
        GL_STENCIL_WRITEMASK
        GL_STENCIL_INDEX
        GL_KEEP
        GL_REPLACE
        GL_INCR
        GL_DECR

# Buffers, Pixel Drawing/Reading
        GL_NONE
        GL_LEFT
        GL_RIGHT
        GL_FRONT_LEFT
        GL_FRONT_RIGHT
        GL_BACK_LEFT
        GL_BACK_RIGHT
        GL_AUX0
        GL_AUX1
        GL_AUX2
        GL_AUX3
        GL_COLOR_INDEX
        GL_RED
        GL_GREEN
        GL_BLUE
        GL_ALPHA
        GL_LUMINANCE
        GL_LUMINANCE_ALPHA
        GL_ALPHA_BITS
        GL_RED_BITS
        GL_GREEN_BITS
        GL_BLUE_BITS
        GL_INDEX_BITS
        GL_SUBPIXEL_BITS
        GL_AUX_BUFFERS
        GL_READ_BUFFER
        GL_DRAW_BUFFER
        GL_DOUBLEBUFFER
        GL_STEREO
        GL_BITMAP
        GL_COLOR
        GL_DEPTH
        GL_STENCIL
        GL_DITHER
        GL_RGB
        GL_RGBA

# Implementation limits
        GL_MAX_LIST_NESTING
        GL_MAX_EVAL_ORDER
        GL_MAX_LIGHTS
        GL_MAX_CLIP_PLANES
        GL_MAX_TEXTURE_SIZE
        GL_MAX_PIXEL_MAP_TABLE
        GL_MAX_ATTRIB_STACK_DEPTH
        GL_MAX_MODELVIEW_STACK_DEPTH
        GL_MAX_NAME_STACK_DEPTH
        GL_MAX_PROJECTION_STACK_DEPTH
        GL_MAX_TEXTURE_STACK_DEPTH
        GL_MAX_VIEWPORT_DIMS
        GL_MAX_CLIENT_ATTRIB_STACK_DEPTH

# Gets
        GL_ATTRIB_STACK_DEPTH
        GL_CLIENT_ATTRIB_STACK_DEPTH
        GL_COLOR_CLEAR_VALUE
        GL_COLOR_WRITEMASK
        GL_CURRENT_INDEX
        GL_CURRENT_COLOR
        GL_CURRENT_NORMAL
        GL_CURRENT_RASTER_COLOR
        GL_CURRENT_RASTER_DISTANCE
        GL_CURRENT_RASTER_INDEX
        GL_CURRENT_RASTER_POSITION
        GL_CURRENT_RASTER_TEXTURE_COORDS
        GL_CURRENT_RASTER_POSITION_VALID
        GL_CURRENT_TEXTURE_COORDS
        GL_INDEX_CLEAR_VALUE
        GL_INDEX_MODE
        GL_INDEX_WRITEMASK
        GL_MODELVIEW_MATRIX
        GL_MODELVIEW_STACK_DEPTH
        GL_NAME_STACK_DEPTH
        GL_PROJECTION_MATRIX
        GL_PROJECTION_STACK_DEPTH
        GL_RENDER_MODE
        GL_RGBA_MODE
        GL_TEXTURE_MATRIX
        GL_TEXTURE_STACK_DEPTH
        GL_VIEWPORT

# Evaluators
        GL_AUTO_NORMAL
        GL_MAP1_COLOR_4
        GL_MAP1_INDEX
        GL_MAP1_NORMAL
        GL_MAP1_TEXTURE_COORD_1
        GL_MAP1_TEXTURE_COORD_2
        GL_MAP1_TEXTURE_COORD_3
        GL_MAP1_TEXTURE_COORD_4
        GL_MAP1_VERTEX_3
        GL_MAP1_VERTEX_4
        GL_MAP2_COLOR_4
        GL_MAP2_INDEX
        GL_MAP2_NORMAL
        GL_MAP2_TEXTURE_COORD_1
        GL_MAP2_TEXTURE_COORD_2
        GL_MAP2_TEXTURE_COORD_3
        GL_MAP2_TEXTURE_COORD_4
        GL_MAP2_VERTEX_3
        GL_MAP2_VERTEX_4
        GL_MAP1_GRID_DOMAIN
        GL_MAP1_GRID_SEGMENTS
        GL_MAP2_GRID_DOMAIN
        GL_MAP2_GRID_SEGMENTS
        GL_COEFF
        GL_ORDER
        GL_DOMAIN

# Hints
        GL_PERSPECTIVE_CORRECTION_HINT
        GL_POINT_SMOOTH_HINT
        GL_LINE_SMOOTH_HINT
        GL_POLYGON_SMOOTH_HINT
        GL_FOG_HINT
        GL_DONT_CARE
        GL_FASTEST
        GL_NICEST

# Scissor box
        GL_SCISSOR_BOX
        GL_SCISSOR_TEST

# Pixel Mode / Transfer
        GL_MAP_COLOR
        GL_MAP_STENCIL
        GL_INDEX_SHIFT
        GL_INDEX_OFFSET
        GL_RED_SCALE
        GL_RED_BIAS
        GL_GREEN_SCALE
        GL_GREEN_BIAS
        GL_BLUE_SCALE
        GL_BLUE_BIAS
        GL_ALPHA_SCALE
        GL_ALPHA_BIAS
        GL_DEPTH_SCALE
        GL_DEPTH_BIAS
        GL_PIXEL_MAP_S_TO_S_SIZE
        GL_PIXEL_MAP_I_TO_I_SIZE
        GL_PIXEL_MAP_I_TO_R_SIZE
        GL_PIXEL_MAP_I_TO_G_SIZE
        GL_PIXEL_MAP_I_TO_B_SIZE
        GL_PIXEL_MAP_I_TO_A_SIZE
        GL_PIXEL_MAP_R_TO_R_SIZE
        GL_PIXEL_MAP_G_TO_G_SIZE
        GL_PIXEL_MAP_B_TO_B_SIZE
        GL_PIXEL_MAP_A_TO_A_SIZE
        GL_PIXEL_MAP_S_TO_S
        GL_PIXEL_MAP_I_TO_I
        GL_PIXEL_MAP_I_TO_R
        GL_PIXEL_MAP_I_TO_G
        GL_PIXEL_MAP_I_TO_B
        GL_PIXEL_MAP_I_TO_A
        GL_PIXEL_MAP_R_TO_R
        GL_PIXEL_MAP_G_TO_G
        GL_PIXEL_MAP_B_TO_B
        GL_PIXEL_MAP_A_TO_A
        GL_PACK_ALIGNMENT
        GL_PACK_LSB_FIRST
        GL_PACK_ROW_LENGTH
        GL_PACK_SKIP_PIXELS
        GL_PACK_SKIP_ROWS
        GL_PACK_SWAP_BYTES
        GL_UNPACK_ALIGNMENT
        GL_UNPACK_LSB_FIRST
        GL_UNPACK_ROW_LENGTH
        GL_UNPACK_SKIP_PIXELS
        GL_UNPACK_SKIP_ROWS
        GL_UNPACK_SWAP_BYTES
        GL_ZOOM_X
        GL_ZOOM_Y

# Texture mapping
        GL_TEXTURE_ENV
        GL_TEXTURE_ENV_MODE
        GL_TEXTURE_1D
        GL_TEXTURE_2D
        GL_TEXTURE_WRAP_S
        GL_TEXTURE_WRAP_T
        GL_TEXTURE_MAG_FILTER
        GL_TEXTURE_MIN_FILTER
        GL_TEXTURE_ENV_COLOR
        GL_TEXTURE_GEN_S
        GL_TEXTURE_GEN_T
        GL_TEXTURE_GEN_MODE
        GL_TEXTURE_BORDER_COLOR
        GL_TEXTURE_WIDTH
        GL_TEXTURE_HEIGHT
        GL_TEXTURE_BORDER
        GL_TEXTURE_COMPONENTS
        GL_TEXTURE_RED_SIZE
        GL_TEXTURE_GREEN_SIZE
        GL_TEXTURE_BLUE_SIZE
        GL_TEXTURE_ALPHA_SIZE
        GL_TEXTURE_LUMINANCE_SIZE
        GL_TEXTURE_INTENSITY_SIZE
        GL_NEAREST_MIPMAP_NEAREST
        GL_NEAREST_MIPMAP_LINEAR
        GL_LINEAR_MIPMAP_NEAREST
        GL_LINEAR_MIPMAP_LINEAR
        GL_OBJECT_LINEAR
        GL_OBJECT_PLANE
        GL_EYE_LINEAR
        GL_EYE_PLANE
        GL_SPHERE_MAP
        GL_DECAL
        GL_MODULATE
        GL_NEAREST
        GL_REPEAT
        GL_CLAMP
        GL_S
        GL_T
        GL_R
        GL_Q
        GL_TEXTURE_GEN_R
        GL_TEXTURE_GEN_Q

# Utility
        GL_VENDOR
        GL_RENDERER
        GL_VERSION
        GL_EXTENSIONS
        GL_SHADING_LANGUAGE_VERSION

# Errors
        GL_NO_ERROR 
        GL_INVALID_ENUM
        GL_INVALID_VALUE
        GL_INVALID_OPERATION
        GL_STACK_OVERFLOW
        GL_STACK_UNDERFLOW
        GL_OUT_OF_MEMORY

# glPush/PopAttrib bits
        GL_CURRENT_BIT
        GL_POINT_BIT
        GL_LINE_BIT
        GL_POLYGON_BIT
        GL_POLYGON_STIPPLE_BIT
        GL_PIXEL_MODE_BIT
        GL_LIGHTING_BIT
        GL_FOG_BIT
        GL_DEPTH_BUFFER_BIT
        GL_ACCUM_BUFFER_BIT
        GL_STENCIL_BUFFER_BIT
        GL_VIEWPORT_BIT
        GL_TRANSFORM_BIT
        GL_ENABLE_BIT
        GL_COLOR_BUFFER_BIT
        GL_HINT_BIT
        GL_EVAL_BIT
        GL_LIST_BIT
        GL_TEXTURE_BIT
        GL_SCISSOR_BIT
        GL_ALL_ATTRIB_BITS

# Vertex Buffer Objects 
        GL_ARRAY_BUFFER
        GL_ELEMENT_ARRAY_BUFFER
        GL_ARRAY_BUFFER_BINDING
        GL_ELEMENT_ARRAY_BUFFER_BINDING
        GL_VERTEX_ARRAY_BUFFER_BINDING
        GL_NORMAL_ARRAY_BUFFER_BINDING
        GL_COLOR_ARRAY_BUFFER_BINDING
        GL_INDEX_ARRAY_BUFFER_BINDING
        GL_TEXTURE_COORD_ARRAY_BUFFER_BINDING
        GL_EDGE_FLAG_ARRAY_BUFFER_BINDING
        GL_SECONDARY_COLOR_ARRAY_BUFFER_BINDING
        GL_FOG_COORD_ARRAY_BUFFER_BINDING
        GL_WEIGHT_ARRAY_BUFFER_BINDING
        GL_VERTEX_ATTRIB_ARRAY_BUFFER_BINDING
        GL_STREAM_DRAW
        GL_STREAM_READ
        GL_STREAM_COPY
        GL_STATIC_DRAW
        GL_STATIC_READ
        GL_STATIC_COPY
        GL_DYNAMIC_DRAW
        GL_DYNAMIC_READ
        GL_DYNAMIC_COPY
        GL_READ_ONLY
        GL_WRITE_ONLY
        GL_READ_WRITE
        GL_BUFFER_SIZE
        GL_BUFFER_USAGE
        GL_BUFFER_ACCESS
        GL_BUFFER_MAPPED
        GL_BUFFER_MAP_POINTER

# Miscellaneous
    cdef void glClearIndex( GLfloat c )
    cdef void glClearColor( GLclampf red, GLclampf green, GLclampf blue, GLclampf alpha )
    cdef void glClear( GLbitfield mask )
    cdef void glIndexMask( GLuint mask )
    cdef void glColorMask( GLboolean red, GLboolean green, GLboolean blue, GLboolean alpha )
    cdef void glAlphaFunc( GLenum func, GLclampf ref )
    cdef void glBlendFunc( GLenum sfactor, GLenum dfactor )
    cdef void glLogicOp( GLenum opcode )
    cdef void glCullFace( GLenum mode )
    cdef void glFrontFace( GLenum mode )
    cdef void glPointSize( GLfloat size )
    cdef void glLineWidth( GLfloat width )
    cdef void glLineStipple( GLint factor, GLushort pattern )
    cdef void glPolygonMode( GLenum face, GLenum mode )
    cdef void glPolygonOffset( GLfloat factor, GLfloat units )
    cdef void glPolygonStipple( GLubyte *mask )
    cdef void glGetPolygonStipple( GLubyte *mask )
    cdef void glEdgeFlag( GLboolean flag )
    cdef void glEdgeFlagv( GLboolean *flag )
    cdef void glScissor( GLint x, GLint y, GLsizei width, GLsizei height)
    cdef void glClipPlane( GLenum plane,  GLdouble *equation )
    cdef void glGetClipPlane( GLenum plane, GLdouble *equation )
    cdef void glDrawBuffer( GLenum mode )
    cdef void glReadBuffer( GLenum mode )
    cdef void glEnable( GLenum cap )
    cdef void glDisable( GLenum cap )
    cdef GLboolean glIsEnabled( GLenum cap )
    cdef void glEnableClientState( GLenum cap )
    cdef void glDisableClientState( GLenum cap )
    cdef void glGetBooleanv( GLenum pname, GLboolean *params )
    cdef void glGetDoublev( GLenum pname, GLdouble *params )
    cdef void glGetFloatv( GLenum pname, GLfloat *params )
    cdef void glGetIntegerv( GLenum pname, GLint *params )
    cdef void glPushAttrib( GLbitfield mask )
    cdef void glPopAttrib( )
    cdef void glPushClientAttrib( GLbitfield mask )
    cdef void glPopClientAttrib( )
    cdef GLint glRenderMode( GLenum mode )
    cdef GLenum glGetError( )
    cdef  GLubyte * glGetString( GLenum name )
    cdef void glFinish( )
    cdef void glFlush( )
    cdef void glHint( GLenum target, GLenum mode )

# Depth Buffer
    cdef void glClearDepth( GLclampd depth )
    cdef void glDepthFunc( GLenum func )
    cdef void glDepthMask( GLboolean flag )
    cdef void glDepthRange( GLclampd near_val, GLclampd far_val )

# Accumulation Buffer
    cdef void glClearAccum( GLfloat red, GLfloat green, GLfloat blue, GLfloat alpha )
    cdef void glAccum( GLenum op, GLfloat value )

# Transformation
    cdef void glMatrixMode( GLenum mode )
    cdef void glOrtho( GLdouble left, GLdouble right,
                                 GLdouble bottom, GLdouble top,
                                 GLdouble near_val, GLdouble far_val )
    cdef void glFrustum( GLdouble left, GLdouble right,
                                   GLdouble bottom, GLdouble top,
                                   GLdouble near_val, GLdouble far_val )
    cdef void glViewport( GLint x, GLint y,
                                    GLsizei width, GLsizei height )
    cdef void glPushMatrix( )
    cdef void glPopMatrix( )
    cdef void glLoadIdentity( )
    cdef void glLoadMatrixd( GLdouble *m )
    cdef void glLoadMatrixf( GLfloat *m )
    cdef void glMultMatrixd( GLdouble *m )
    cdef void glMultMatrixf( GLfloat *m )
    cdef void glRotated( GLdouble angle, GLdouble x, GLdouble y, GLdouble z )
    cdef void glRotatef( GLfloat angle, GLfloat x, GLfloat y, GLfloat z )
    cdef void glScaled( GLdouble x, GLdouble y, GLdouble z )
    cdef void glScalef( GLfloat x, GLfloat y, GLfloat z )
    cdef void glTranslated( GLdouble x, GLdouble y, GLdouble z )
    cdef void glTranslatef( GLfloat x, GLfloat y, GLfloat z )

# Display Lists
    cdef GLboolean glIsList( GLuint list )
    cdef void glDeleteLists( GLuint list, GLsizei range )
    cdef GLuint glGenLists( GLsizei range )
    cdef void glNewList( GLuint list, GLenum mode )
    cdef void glEndList( )
    cdef void glCallList( GLuint list )
    cdef void glCallLists( GLsizei n, GLenum type,  GLvoid *lists )
    cdef void glListBase( GLuint base )

# Drawing Functions
    cdef void glBegin( GLenum mode )
    cdef void glEnd( )
    cdef void glVertex2d( GLdouble x, GLdouble y )
    cdef void glVertex2f( GLfloat x, GLfloat y )
    cdef void glVertex2i( GLint x, GLint y )
    cdef void glVertex2s( GLshort x, GLshort y )
    cdef void glVertex3d( GLdouble x, GLdouble y, GLdouble z )
    cdef void glVertex3f( GLfloat x, GLfloat y, GLfloat z )
    cdef void glVertex3i( GLint x, GLint y, GLint z )
    cdef void glVertex3s( GLshort x, GLshort y, GLshort z )
    cdef void glVertex4d( GLdouble x, GLdouble y, GLdouble z, GLdouble w )
    cdef void glVertex4f( GLfloat x, GLfloat y, GLfloat z, GLfloat w )
    cdef void glVertex4i( GLint x, GLint y, GLint z, GLint w )
    cdef void glVertex4s( GLshort x, GLshort y, GLshort z, GLshort w )
    cdef void glVertex2dv( GLdouble *v )
    cdef void glVertex2fv( GLfloat *v )
    cdef void glVertex2iv( GLint *v )
    cdef void glVertex2sv( GLshort *v )
    cdef void glVertex3dv( GLdouble *v )
    cdef void glVertex3fv( GLfloat *v )
    cdef void glVertex3iv( GLint *v )
    cdef void glVertex3sv( GLshort *v )
    cdef void glVertex4dv( GLdouble *v )
    cdef void glVertex4fv( GLfloat *v )
    cdef void glVertex4iv( GLint *v )
    cdef void glVertex4sv( GLshort *v )

    cdef void glNormal3b( GLbyte nx, GLbyte ny, GLbyte nz )
    cdef void glNormal3d( GLdouble nx, GLdouble ny, GLdouble nz )
    cdef void glNormal3f( GLfloat nx, GLfloat ny, GLfloat nz )
    cdef void glNormal3i( GLint nx, GLint ny, GLint nz )
    cdef void glNormal3s( GLshort nx, GLshort ny, GLshort nz )
    cdef void glNormal3bv( GLbyte *v )
    cdef void glNormal3dv( GLdouble *v )
    cdef void glNormal3fv( GLfloat *v )
    cdef void glNormal3iv( GLint *v )
    cdef void glNormal3sv( GLshort *v )

    cdef void glIndexd( GLdouble c )
    cdef void glIndexf( GLfloat c )
    cdef void glIndexi( GLint c )
    cdef void glIndexs( GLshort c )
    cdef void glIndexub( GLubyte c )
    cdef void glIndexdv( GLdouble *c )
    cdef void glIndexfv( GLfloat *c )
    cdef void glIndexiv( GLint *c )
    cdef void glIndexsv( GLshort *c )
    cdef void glIndexubv( GLubyte *c )

    cdef void glColor3b( GLbyte red, GLbyte green, GLbyte blue )
    cdef void glColor3d( GLdouble red, GLdouble green, GLdouble blue )
    cdef void glColor3f( GLfloat red, GLfloat green, GLfloat blue )
    cdef void glColor3i( GLint red, GLint green, GLint blue )
    cdef void glColor3s( GLshort red, GLshort green, GLshort blue )
    cdef void glColor3ub( GLubyte red, GLubyte green, GLubyte blue )
    cdef void glColor3ui( GLuint red, GLuint green, GLuint blue )
    cdef void glColor3us( GLushort red, GLushort green, GLushort blue )
    cdef void glColor4b( GLbyte red, GLbyte green, GLbyte blue, GLbyte alpha )
    cdef void glColor4d( GLdouble red, GLdouble green, GLdouble blue, GLdouble alpha )
    cdef void glColor4f( GLfloat red, GLfloat green, GLfloat blue, GLfloat alpha )
    cdef void glColor4i( GLint red, GLint green, GLint blue, GLint alpha )
    cdef void glColor4s( GLshort red, GLshort green, GLshort blue, GLshort alpha )
    cdef void glColor4ub( GLubyte red, GLubyte green, GLubyte blue, GLubyte alpha )
    cdef void glColor4ui( GLuint red, GLuint green, GLuint blue, GLuint alpha )
    cdef void glColor4us( GLushort red, GLushort green, GLushort blue, GLushort alpha )

    cdef void glColor3bv( GLbyte *v )
    cdef void glColor3dv( GLdouble *v )
    cdef void glColor3fv( GLfloat *v )
    cdef void glColor3iv( GLint *v )
    cdef void glColor3sv( GLshort *v )
    cdef void glColor3ubv( GLubyte *v )
    cdef void glColor3uiv( GLuint *v )
    cdef void glColor3usv( GLushort *v )
    cdef void glColor4bv( GLbyte *v )
    cdef void glColor4dv( GLdouble *v )
    cdef void glColor4fv( GLfloat *v )
    cdef void glColor4iv( GLint *v )
    cdef void glColor4sv( GLshort *v )
    cdef void glColor4ubv( GLubyte *v )
    cdef void glColor4uiv( GLuint *v )
    cdef void glColor4usv( GLushort *v )

    cdef void glTexCoord1d( GLdouble s )
    cdef void glTexCoord1f( GLfloat s )
    cdef void glTexCoord1i( GLint s )
    cdef void glTexCoord1s( GLshort s )
    cdef void glTexCoord2d( GLdouble s, GLdouble t )
    cdef void glTexCoord2f( GLfloat s, GLfloat t )
    cdef void glTexCoord2i( GLint s, GLint t )
    cdef void glTexCoord2s( GLshort s, GLshort t )
    cdef void glTexCoord3d( GLdouble s, GLdouble t, GLdouble r )
    cdef void glTexCoord3f( GLfloat s, GLfloat t, GLfloat r )
    cdef void glTexCoord3i( GLint s, GLint t, GLint r )
    cdef void glTexCoord3s( GLshort s, GLshort t, GLshort r )
    cdef void glTexCoord4d( GLdouble s, GLdouble t, GLdouble r, GLdouble q )
    cdef void glTexCoord4f( GLfloat s, GLfloat t, GLfloat r, GLfloat q )
    cdef void glTexCoord4i( GLint s, GLint t, GLint r, GLint q )
    cdef void glTexCoord4s( GLshort s, GLshort t, GLshort r, GLshort q )
    cdef void glTexCoord1dv( GLdouble *v )
    cdef void glTexCoord1fv( GLfloat *v )
    cdef void glTexCoord1iv( GLint *v )
    cdef void glTexCoord1sv( GLshort *v )
    cdef void glTexCoord2dv( GLdouble *v )
    cdef void glTexCoord2fv( GLfloat *v )
    cdef void glTexCoord2iv( GLint *v )
    cdef void glTexCoord2sv( GLshort *v )
    cdef void glTexCoord3dv( GLdouble *v )
    cdef void glTexCoord3fv( GLfloat *v )
    cdef void glTexCoord3iv( GLint *v )
    cdef void glTexCoord3sv( GLshort *v )
    cdef void glTexCoord4dv( GLdouble *v )
    cdef void glTexCoord4fv( GLfloat *v )
    cdef void glTexCoord4iv( GLint *v )
    cdef void glTexCoord4sv( GLshort *v )

    cdef void glRasterPos2d( GLdouble x, GLdouble y )
    cdef void glRasterPos2f( GLfloat x, GLfloat y )
    cdef void glRasterPos2i( GLint x, GLint y )
    cdef void glRasterPos2s( GLshort x, GLshort y )
    cdef void glRasterPos3d( GLdouble x, GLdouble y, GLdouble z )
    cdef void glRasterPos3f( GLfloat x, GLfloat y, GLfloat z )
    cdef void glRasterPos3i( GLint x, GLint y, GLint z )
    cdef void glRasterPos3s( GLshort x, GLshort y, GLshort z )
    cdef void glRasterPos4d( GLdouble x, GLdouble y, GLdouble z, GLdouble w )
    cdef void glRasterPos4f( GLfloat x, GLfloat y, GLfloat z, GLfloat w )
    cdef void glRasterPos4i( GLint x, GLint y, GLint z, GLint w )
    cdef void glRasterPos4s( GLshort x, GLshort y, GLshort z, GLshort w )
    cdef void glRasterPos2dv( GLdouble *v )
    cdef void glRasterPos2fv( GLfloat *v )
    cdef void glRasterPos2iv( GLint *v )
    cdef void glRasterPos2sv( GLshort *v )
    cdef void glRasterPos3dv( GLdouble *v )
    cdef void glRasterPos3fv( GLfloat *v )
    cdef void glRasterPos3iv( GLint *v )
    cdef void glRasterPos3sv( GLshort *v )
    cdef void glRasterPos4dv( GLdouble *v )
    cdef void glRasterPos4fv( GLfloat *v )
    cdef void glRasterPos4iv( GLint *v )
    cdef void glRasterPos4sv( GLshort *v )

    cdef void glRectd( GLdouble x1, GLdouble y1, GLdouble x2, GLdouble y2 )
    cdef void glRectf( GLfloat x1, GLfloat y1, GLfloat x2, GLfloat y2 )
    cdef void glRecti( GLint x1, GLint y1, GLint x2, GLint y2 )
    cdef void glRects( GLshort x1, GLshort y1, GLshort x2, GLshort y2 )

    cdef void glRectdv( GLdouble *v1, GLdouble *v2 )
    cdef void glRectfv( GLfloat *v1, GLfloat *v2 )
    cdef void glRectiv( GLint *v1, GLint *v2 )
    cdef void glRectsv( GLshort *v1, GLshort *v2 )

# Vertex Arrays  (1.1)
    cdef void glVertexPointer( GLint size, GLenum type, GLsizei stride, GLvoid *ptr )
    cdef void glNormalPointer( GLenum type, GLsizei stride, GLvoid *ptr )
    cdef void glColorPointer( GLint size, GLenum type, GLsizei stride, GLvoid *ptr )
    cdef void glIndexPointer( GLenum type, GLsizei stride, GLvoid *ptr )
    cdef void glTexCoordPointer( GLint size, GLenum type, GLsizei stride, GLvoid *ptr )
    cdef void glEdgeFlagPointer( GLsizei stride, GLvoid *ptr )
    cdef void glGetPointerv( GLenum pname, GLvoid **params )
    cdef void glArrayElement( GLint i )
    cdef void glDrawArrays( GLenum mode, GLint first, GLsizei count )
    cdef void glDrawElements( GLenum mode, GLsizei count, GLenum type, GLvoid *indices )
    cdef void glInterleavedArrays( GLenum format, GLsizei stride, GLvoid *pointer )

# Lighting
    cdef void glShadeModel( GLenum mode )
    cdef void glLightf( GLenum light, GLenum pname, GLfloat param )
    cdef void glLighti( GLenum light, GLenum pname, GLint param )
    cdef void glLightfv( GLenum light, GLenum pname, GLfloat *params )
    cdef void glLightiv( GLenum light, GLenum pname, GLint *params )
    cdef void glGetLightfv( GLenum light, GLenum pname, GLfloat *params )
    cdef void glGetLightiv( GLenum light, GLenum pname, GLint *params )
    cdef void glLightModelf( GLenum pname, GLfloat param )
    cdef void glLightModeli( GLenum pname, GLint param )
    cdef void glLightModelfv( GLenum pname, GLfloat *params )
    cdef void glLightModeliv( GLenum pname, GLint *params )
    cdef void glMaterialf( GLenum face, GLenum pname, GLfloat param )
    cdef void glMateriali( GLenum face, GLenum pname, GLint param )
    cdef void glMaterialfv( GLenum face, GLenum pname, GLfloat *params )
    cdef void glMaterialiv( GLenum face, GLenum pname, GLint *params )
    cdef void glGetMaterialfv( GLenum face, GLenum pname, GLfloat *params )
    cdef void glGetMaterialiv( GLenum face, GLenum pname, GLint *params )
    cdef void glColorMaterial( GLenum face, GLenum mode )

# Raster functions
    cdef void glPixelZoom( GLfloat xfactor, GLfloat yfactor )
    cdef void glPixelStoref( GLenum pname, GLfloat param )
    cdef void glPixelStorei( GLenum pname, GLint param )
    cdef void glPixelTransferf( GLenum pname, GLfloat param )
    cdef void glPixelTransferi( GLenum pname, GLint param )
    cdef void glPixelMapfv( GLenum map, GLsizei mapsize, GLfloat *values )
    cdef void glPixelMapuiv( GLenum map, GLsizei mapsize, GLuint *values )
    cdef void glPixelMapusv( GLenum map, GLsizei mapsize, GLushort *values )
    cdef void glGetPixelMapfv( GLenum map, GLfloat *values )
    cdef void glGetPixelMapuiv( GLenum map, GLuint *values )
    cdef void glGetPixelMapusv( GLenum map, GLushort *values )
    cdef void glBitmap( GLsizei width, GLsizei height,
                        GLfloat xorig, GLfloat yorig,
                        GLfloat xmove, GLfloat ymove,
                        GLubyte *bitmap )
    cdef void glReadPixels( GLint x, GLint y,
                            GLsizei width, GLsizei height,
                            GLenum format, GLenum type,
                            GLvoid *pixels )
    cdef void glDrawPixels( GLsizei width, GLsizei height,
                            GLenum format, GLenum type,
                            GLvoid *pixels )
    cdef void glCopyPixels( GLint x, GLint y,
                            GLsizei width, GLsizei height,
                            GLenum type )
# Stenciling
    cdef void glStencilFunc( GLenum func, GLint ref, GLuint mask )
    cdef void glStencilMask( GLuint mask )
    cdef void glStencilOp( GLenum fail, GLenum zfail, GLenum zpass )
    cdef void glClearStencil( GLint s )

# Texture mapping
    cdef void glTexGend( GLenum coord, GLenum pname, GLdouble param )
    cdef void glTexGenf( GLenum coord, GLenum pname, GLfloat param )
    cdef void glTexGeni( GLenum coord, GLenum pname, GLint param )
    cdef void glTexGendv( GLenum coord, GLenum pname, GLdouble *params )
    cdef void glTexGenfv( GLenum coord, GLenum pname, GLfloat *params )
    cdef void glTexGeniv( GLenum coord, GLenum pname, GLint *params )
    cdef void glGetTexGendv( GLenum coord, GLenum pname, GLdouble *params )
    cdef void glGetTexGenfv( GLenum coord, GLenum pname, GLfloat *params )
    cdef void glGetTexGeniv( GLenum coord, GLenum pname, GLint *params )

    cdef void glTexEnvf( GLenum target, GLenum pname, GLfloat param )
    cdef void glTexEnvi( GLenum target, GLenum pname, GLint param )
    cdef void glTexEnvfv( GLenum target, GLenum pname, GLfloat *params )
    cdef void glTexEnviv( GLenum target, GLenum pname, GLint *params )
    cdef void glGetTexEnvfv( GLenum target, GLenum pname, GLfloat *params )
    cdef void glGetTexEnviv( GLenum target, GLenum pname, GLint *params )

    cdef void glTexParameterf( GLenum target, GLenum pname, GLfloat param )
    cdef void glTexParameteri( GLenum target, GLenum pname, GLint param )
    cdef void glTexParameterfv( GLenum target, GLenum pname, GLfloat *params )
    cdef void glTexParameteriv( GLenum target, GLenum pname, GLint *params )
    cdef void glGetTexParameterfv( GLenum target, GLenum pname, GLfloat *params)
    cdef void glGetTexParameteriv( GLenum target, GLenum pname, GLint *params )
    cdef void glGetTexLevelParameterfv( GLenum target, GLint level, GLenum pname, GLfloat *params )
    cdef void glGetTexLevelParameteriv( GLenum target, GLint level, GLenum pname, GLint *params )

    cdef void glTexImage1D( GLenum target, GLint level,
                            GLint internalFormat,
                            GLsizei width, GLint border,
                            GLenum format, GLenum type,
                            GLvoid *pixels )
    cdef void glTexImage2D( GLenum target, GLint level,
                            GLint internalFormat,
                            GLsizei width, GLsizei height,
                            GLint border, GLenum format, GLenum type,
                            GLvoid *pixels )
    cdef void glGetTexImage( GLenum target, GLint level,
                             GLenum format, GLenum type,
                             GLvoid *pixels )

# Vertex Buffer Objects

    cdef void glBindBuffer (GLenum target, GLuint buffer)
    cdef void glDeleteBuffers (GLsizei n, GLuint *buffers)
    cdef void glGenBuffers (GLsizei n, GLuint *buffers)
    cdef GLboolean glIsBuffer (GLuint buffer)
    cdef void glBufferData (GLenum target, GLsizeiptr size, GLvoid *data,
                            GLenum usage)
    cdef void glBufferSubData (GLenum target, GLintptr offset,
                               GLsizeiptr size, GLvoid *data)
    cdef void glGetBufferSubData (GLenum target, GLintptr offset,
                                  GLsizeiptr size, GLvoid *data)
    cdef GLvoid * glMapBuffer (GLenum target, GLenum access)
    cdef GLboolean glUnmapBuffer (GLenum target)
    cdef void glGetBufferParameteriv (GLenum target, GLenum pname, GLint *params)
    cdef void glGetBufferPointerv (GLenum target, GLenum pname, GLvoid **params)

# Pixmap Font
cdef extern from "SnapPyfont.h":
    ctypedef struct SnapPy_glyph:
        int     width
        int     height
        int     bytes_per_pixel 
        char*   pixel_data
    cdef SnapPy_glyph* SnapPy_font[]
