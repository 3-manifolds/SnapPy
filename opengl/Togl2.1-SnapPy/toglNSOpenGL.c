/* $Id: toglNSOpenGL.c,v 1.7 2009/10/22 00:06:41 gregcouch Exp $ */

/* vi:set sw=4 expandtab: */

/* 
 * Togl - a Tk OpenGL widget
 *
 * Copyright (C) 1996-2002  Brian Paul and Ben Bederson
 * Copyright (C) 2005-2009  Greg Couch
 * See the LICENSE file for copyright details.
 */

static NSOpenGLPixelFormat *
togl_pixelFormat(Togl *togl)
{
    NSOpenGLPixelFormatAttribute   attribs[32];
    int     na = 0;
    NSOpenGLPixelFormat *pix;

#if 0
    if (togl->MultisampleFlag && !hasMultisampling) {
        Tcl_SetResult(togl->Interp,
                "multisampling not supported", TCL_STATIC);
        return NULL;
    }
#endif

    attribs[na++] = NSOpenGLPFAMinimumPolicy;
    /* ask for hardware-accelerated onscreen */
    /* This is not needed, and can break virtual machines.
       Accelerated rendering is always preferred.
    attribs[na++] = NSOpenGLPFAAccelerated;
    attribs[na++] = NSOpenGLPFANoRecovery;
    */
    if (togl->RgbaFlag) {
        /* RGB[A] mode */
        attribs[na++] = NSOpenGLPFAColorSize;
	attribs[na++] = togl->RgbaRed + togl->RgbaGreen + togl->RgbaBlue;
	/* NSOpenGL does not take separate red,green,blue sizes. */
        if (togl->AlphaFlag) {
            attribs[na++] = NSOpenGLPFAAlphaSize;
            attribs[na++] = togl->AlphaSize;
        }
    } else {
        /* Color index mode */
        Tcl_SetResult(togl->Interp,
                "Color index mode not supported", TCL_STATIC);
        return NULL;
    }
    if (togl->DepthFlag) {
        attribs[na++] = NSOpenGLPFADepthSize;
        attribs[na++] = togl->DepthSize;
    }
    if (togl->DoubleFlag) {
        attribs[na++] = NSOpenGLPFADoubleBuffer;
    }
    if (togl->StencilFlag) {
        attribs[na++] = NSOpenGLPFAStencilSize;
        attribs[na++] = togl->StencilSize;
    }
    if (togl->AccumFlag) {
        attribs[na++] = NSOpenGLPFAAccumSize;
        attribs[na++] = togl->AccumRed + togl->AccumGreen + togl->AccumBlue + (togl->AlphaFlag ? togl->AccumAlpha : 0);
    }
    if (togl->MultisampleFlag) {
        attribs[na++] = NSOpenGLPFAMultisample;
        attribs[na++] = NSOpenGLPFASampleBuffers;
        attribs[na++] = 1;
        attribs[na++] = NSOpenGLPFASamples;
        attribs[na++] = 2;
    }
    if (togl->AuxNumber != 0) {
        attribs[na++] = NSOpenGLPFAAuxBuffers;
        attribs[na++] = togl->AuxNumber;
    }
    if (togl->Stereo == TOGL_STEREO_NATIVE) {
        attribs[na++] = NSOpenGLPFAStereo;
    }
    if (togl->FullscreenFlag) {
        Tcl_SetResult(togl->Interp,
                "FullScreen mode not supported.", TCL_STATIC);
        return NULL;
    }
    switch(togl->profile) {
    case PROFILE_3_2:
      attribs[na++] = NSOpenGLPFAOpenGLProfile;
      attribs[na++] = NSOpenGLProfileVersion3_2Core;
      break;
    case PROFILE_4_1:
      attribs[na++] = NSOpenGLPFAOpenGLProfile;
      attribs[na++] = NSOpenGLProfileVersion4_1Core;
      break;
    default:
      break;
    }
    attribs[na++] = 0;	/* End of attributes. */

    pix = [[NSOpenGLPixelFormat alloc] initWithAttributes:attribs];
    if (pix == nil) {
        Tcl_SetResult(togl->Interp, "couldn't choose pixel format",
                TCL_STATIC);
        return NULL;
    }
    return pix;
}

static int
togl_describePixelFormat(Togl *togl)
{
    NSOpenGLPixelFormat *pfmt = togl->PixelFormat;

    /* fill in RgbaFlag, DoubleFlag, and Stereo */
    GLint   has_rgba, has_doublebuf, has_depth, has_accum, has_alpha,
            has_stencil, has_stereo, has_multisample;

    GLint   vscr = 0;
    [pfmt getValues:&has_rgba forAttribute:NSOpenGLPFAColorSize forVirtualScreen:vscr];
    [pfmt getValues:&has_doublebuf forAttribute:NSOpenGLPFADoubleBuffer forVirtualScreen:vscr];
    [pfmt getValues:&has_depth forAttribute:NSOpenGLPFADepthSize forVirtualScreen:vscr];
    [pfmt getValues:&has_accum forAttribute:NSOpenGLPFAAccumSize forVirtualScreen:vscr];
    [pfmt getValues:&has_alpha forAttribute:NSOpenGLPFAAlphaSize forVirtualScreen:vscr];
    [pfmt getValues:&has_stencil forAttribute:NSOpenGLPFAStencilSize forVirtualScreen:vscr];
    [pfmt getValues:&has_stereo forAttribute:NSOpenGLPFAStereo forVirtualScreen:vscr];
    [pfmt getValues:&has_multisample forAttribute:NSOpenGLPFASampleBuffers forVirtualScreen:vscr];

    togl->RgbaFlag = (has_rgba != 0);
    togl->DoubleFlag = (has_doublebuf != 0);
    togl->DepthFlag = (has_depth != 0);
    togl->AccumFlag = (has_accum != 0);
    togl->AlphaFlag = (has_alpha != 0);
    togl->StencilFlag = (has_stencil != 0);
    togl->Stereo = (has_stereo ? TOGL_STEREO_NATIVE : TOGL_STEREO_NONE);
    togl->MultisampleFlag = (has_multisample != 0);
    return True;
}

#define isPow2(x) (((x) & ((x) - 1)) == 0)

