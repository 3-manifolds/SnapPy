#!/usr/bin/env python

from Tkinter import Widget, Misc, Tk
from OpenGL.GL import *
from OpenGL.GLU import *
import OpenGL

import os, sys

def glTranslateScene(s, x, y, mousex, mousey):
        glMatrixMode(GL_MODELVIEW)
        mat = glGetDoublev(GL_MODELVIEW_MATRIX)
        glLoadIdentity()
        glTranslatef(s * (x - mousex), s * (mousey - y), 0.0)
        glMultMatrixd(mat)


def glRotateScene(s, xcenter, ycenter, zcenter, x, y, mousex, mousey):
        glMatrixMode(GL_MODELVIEW)
        mat = glGetDoublev(GL_MODELVIEW_MATRIX)
        glLoadIdentity()
        glTranslatef(xcenter, ycenter, zcenter)
        glRotatef(s * (y - mousey), 1., 0., 0.)
        glRotatef(s * (x - mousex), 0., 1., 0.)
        glTranslatef(-xcenter, -ycenter, -zcenter)
        glMultMatrixd(mat)

class RawOpengl(Widget, Misc):
  """
  Widget without any sophisticated bindings
  by Tom Schwaller
  """

  def __init__(self, master, cnf={}, **kw):
    Togl_path = os.path.join( os.path.dirname(__file__),
			      sys.platform + "-tk" + master.getvar("tk_version"))
    master.tk.call('lappend', 'auto_path', Togl_path)
    master.tk.call('package', 'require', 'Togl')

    Widget.__init__(self, master, 'togl', cnf, kw)
    self.root = master
    self.bind('<Map>', self.tkMap)
    self.bind('<Expose>', self.tkExpose)
    self.bind('<Configure>', self.tkExpose)

  def tkRedraw(self, *dummy):
    self.tk.call(self._w, 'makecurrent')
    glPushMatrix()
    self.update_idletasks()
    self.redraw()
    glFlush()
    glPopMatrix()
    self.tk.call(self._w, 'swapbuffers')

  def tkMap(self, *dummy):
    self.tkExpose()

  def tkExpose(self, *dummy):
    self.tkRedraw()

class Opengl(RawOpengl):
  """
  Tkinter bindings for an Opengl widget.
  Mike Hartshorn
  Department of Chemistry
  University of York, UK
  http://www.yorvic.york.ac.uk/~mjh/
  """

  def __init__(self, master=None, help='No help is available.', cnf={}, **kw):
    """
    Create an opengl widget.  Arrange for redraws when the window is
    exposed or when it changes size.
    """

    apply(RawOpengl.__init__, (self, master, cnf), kw)
    self.help_text = help
    self.initialised = 0

    # Current coordinates of the mouse.
    self.xmouse = 0
    self.ymouse = 0

    # Where we are centering.
    self.xcenter = 0.0
    self.ycenter = 0.0
    self.zcenter = 0.0

    # The _back color
    self.r_back = 1.
    self.g_back = 0.
    self.b_back = 1.

    # Where the eye is
    self.distance = 10.0

    # Field of view in y direction
    self.fovy = 30.0

    # Position of clipping planes.
    self.near = 1.0
    self.far = 100.0

    # Is the widget allowed to autospin?
    self.autospin_allowed = 0

    # Is the widget currently autospinning?
    self.autospin = 0

    # Dictionary of key actions (keysym:function) .
    self.key_action = {}

    # Bindings for events.
    self.bind('<Map>', self.tkMap)
    self.bind('<Expose>', self.tkExpose)
    self.bind('<Configure>', self.tkExpose)
    self.bind('<Shift-Button-1>', self.tkHandlePick)
#    self.bind('<Button-1><ButtonRelease-1>', self.tkHandlePick)
#    self.bind('<Button-2>', self.tkRecordMouse)
#    self.bind('<B2-Motion>', self.tkTranslate)
    self.bind('<Button-1>', self.StartRotate)
    self.bind('<B1-Motion>', self.tkRotate)
    self.bind('<ButtonRelease-1>', self.tkAutoSpin)
#    self.bind('<Button-3>', self.tkRecordMouse)
#    self.bind('<B3-Motion>', self.tkScale)
#    self.bind('<KeyPress>', self.tkKeyPress)

  def help(self):
    """Help message for the widget."""

    import tkMessageBox
    tkMessageBox.showinfo('Viewer Help', self.help_text)

  def activate(self):
    """Cause this Opengl widget to be the current destination for
       drawing, and to be the focus of keyboard events."""

    self.tk.call(self._w, 'makecurrent')
    self.focus_set()

  def set_background(self, r, g, b):
    """Change the background colour of the widget."""

    self.r_back = r
    self.g_back = g
    self.b_back = b

    self.tkRedraw()

  def set_centerpoint(self, x, y, z):
    """Set the new center point for the model.
    This is where we are looking."""

    self.xcenter = x
    self.ycenter = y
    self.zcenter = z

    self.tkRedraw()

  def set_eyepoint(self, distance):
    """Set how far the eye is from the position we are looking."""

    self.distance = distance
    self.tkRedraw()

  def reset(self):
    """Reset rotation matrix for this widget."""

    self.autospin = 0
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity()
    self.tkRedraw()

  def tkHandlePick(self, event):
    """Handle a pick on the scene."""

    if hasattr(self, 'pick'):
      # here we need to use glu.UnProject
      # Tk and X have their origin top left, 
      # while Opengl has its origin bottom left.
      # So we need to subtract y from the window height to get
      # the proper pick position for Opengl

      realy = self.winfo_height() - event.y

      p1 = gluUnProject(event.x, realy, 0.)
      p2 = gluUnProject(event.x, realy, 1.)

      if self.pick(self, p1, p2):
	"""If the pick method returns true we redraw the scene."""

	self.tkRedraw()

  def tkRecordMouse(self, event):
    """Record the current mouse position."""

    self.xmouse = event.x
    self.ymouse = event.y

  def StartRotate(self, event):

    # Switch off any autospinning if it was happening
    self.autospin = 0
    self.tkRecordMouse(event)

  def tkScale(self, event):
    """Scale the scene.  Achieved by moving the eye position."""

    scale = 1 - 0.01 * (event.y - self.ymouse)
    self.distance = self.distance * scale
    self.tkRedraw()
    self.tkRecordMouse(event)

  def do_AutoSpin(self):
    s = 0.1
    self.activate()

    glRotateScene(s,
		   self.xcenter, self.ycenter, self.zcenter,
		   self.yspin, self.xspin, 0, 0)
    self.tkRedraw()

    if self.autospin:
      self.after(10, self.do_AutoSpin)

  def tkAutoSpin(self, event):
    """Perform autospin of scene."""

    self.after(16)
    self.update_idletasks()

    # This could be done with one call to pointerxy but I'm not sure
    # it would any quicker as we would have to split up the resulting
    # string and then conv

    x = self.tk.getint(self.tk.call('winfo', 'pointerx', self._w))
    y = self.tk.getint(self.tk.call('winfo', 'pointery', self._w))

    if self.autospin_allowed:
      if x != event.x_root and y != event.y_root:
	self.autospin = 1

	self.yspin = x - event.x_root
	self.xspin = y - event.y_root

	self.after(10, self.do_AutoSpin)

  def tkRotate(self, event):
    """Perform rotation of scene."""

    self.activate()
    glRotateScene(0.5,
		   self.xcenter, self.ycenter, self.zcenter,
		   event.x, event.y, self.xmouse, self.ymouse)
    self.tkRedraw()
    self.tkRecordMouse(event)

  def tkTranslate(self, event):
    """Perform translation of scene."""

    self.activate()
    glTranslateScene(0.05, event.x, event.y, self.xmouse, self.ymouse)
    self.tkRedraw()
    self.tkRecordMouse(event)

  def tkRedraw(self, *dummy):
    """Cause the opengl widget to redraw itself."""

    if not self.initialised: return
    self.activate()

    glPushMatrix()			# Protect our matrix
    self.update_idletasks()
    w = self.winfo_width()
    h = self.winfo_height()
    glViewport(0, 0, w, h)

    # Clear the background and depth buffer.
    glClearColor(self.r_back, self.g_back, self.b_back, 0.)
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity()
    gluPerspective(self.fovy, float(w)/float(h), self.near, self.far)
    gluLookAt(self.xcenter, self.ycenter, self.zcenter + self.distance,
	 self.xcenter, self.ycenter, self.zcenter, 0., 1., 0.)
    glMatrixMode(GL_MODELVIEW);

    # Call objects redraw method.
    self.redraw(self)
    glFlush()				# Tidy up
    glPopMatrix()			# Restore the matrix

    self.tk.call(self._w, 'swapbuffers')

  def tkMap(self, *dummy):
    """Cause the opengl widget to redraw itself."""

    self.tkExpose()

  def tkExpose(self, *dummy):
    """
    Redraw the widget.  Make it active, update tk events, call redraw
    procedure and swap the buffers.  Note: swapbuffers is clever
    enough to only swap double buffered visuals.
    """

    self.activate()
    if not self.initialised:
      self.initialised = 1
    self.tkRedraw()

  def tkKeyPress(self, event):
    """
    Handle keyboard events.
    """

    try:
      self.key_action[event.keysym]()
    except KeyError:
      pass
    if not self.autospin:
      self.tkRedraw()

  def tkPrint(self, file):
    """
    Turn the current scene into PostScript via the feedback buffer.
    """

    self.activate()
