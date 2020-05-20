import sys, os, tkinter as tk
import tkinter
from tkinter import ttk

if sys.platform == 'linux':
    label_pad = (4, 0)
else:
    label_pad = 0


class Slider(ttk.Scale):
    _slider_left_end = 1.0
    _slider_right_end = 100.0

    """
    Subclass of ttk.Scale which makes clicking on the track cause the
    knob to move to the click point.
    """
    def __init__(self, master, left_end, right_end,
                 orient = tkinter.HORIZONTAL):
        ttk.Scale.__init__(self,
                           master = master,
                           from_ = self._slider_left_end,
                           to = self._slider_right_end,
                           value = 2.0,
                           orient = orient,
                           takefocus=0)
        self.left_end = left_end
        self.right_end = right_end
        self.orient = orient
        self.callback = None
        self.configure(command = self._command)
        self.bind('<Button-1>', self._handle_mouse)

    def set_callback(self, callback):
        self.callback = callback

    def set_value(self, value):
        length = self.right_end - self.left_end
        slider_length = self._slider_right_end - self._slider_left_end
        v = (value - self.left_end) / length
        slider_value = v * slider_length + self._slider_left_end
        self.configure(value = slider_value)

    def _command(self, slider_value):
        if self.callback:
            length = self.right_end - self.left_end
            slider_length = self._slider_right_end - self._slider_left_end
            v = (float(slider_value) - self._slider_left_end) / slider_length
            value = v * length + self.left_end

            self.callback(value)

    def _handle_mouse(self, event):
        """
        Override the standard mouse event handler which moves the knob
        by 1.0, regardless of the size of the range, when the mouse
        click is in the trough.  This implements the standard behavior
        -- moving the knob to the click point.
        """
        part = self.identify(event.x, event.y)
        orientation = str(self.cget('orient'))
        if (part == 'Scale.trough'):
            length = self.right_end - self.left_end
            slider_length = self._slider_right_end - self._slider_left_end
            if orientation == tk.HORIZONTAL:
                fraction = float(event.x) / float(self.winfo_width())
            else:
                fraction = float(event.y) / float(self.winfo_height())
            if self.callback:
                self.callback(length * fraction + self.left_end)
            self.set(slider_length * fraction + self._slider_left_end)
            return 'break'

class ZoomSlider(ttk.Frame):
    """
    A compound widget containing a Slider, labels and two buttons that
    expand or contract the range of values by a factor of 2.  The nonneg
    option can be set to ensure that the range remains nonnegative. This
    widget assumes that from < to.
    """
    min_span = 0.05
    max_span = 100

    def __init__(self, master, left_end, right_end, label_text=None, on_change=None):
        self._build_icons()
        ttk.Frame.__init__(self, master)
        self.left_end = left_end
        self.right_end = right_end
        self.on_change = on_change
        self.slider = Slider(self, left_end, right_end)
        self.slider.set_callback(self._slider_callback)

        style_args = ( { 'style':'Toolbutton' }
                       if (sys.platform == 'darwin')
                       else {})

        self.compresser = ttk.Button(self,
                                     image=self.compress_icon,
                                     command=self.zoom_in,
                                     **style_args)
        self.expander = ttk.Button(self,
                                   image=self.expand_icon,
                                   command=self.zoom_out,
                                   **style_args)

        padding_cheat = -6 if (sys.platform == 'darwin') else 0

        if not label_text is None:
            self.title_label = ttk.Label(self, text=label_text, padding=label_pad)
            col_base = 1
            self.columnconfigure(0, weight=0)
            self.title_label.grid(row=0, column=0, sticky=tk.W)
        else:
            self.title_lable = None
            col_base = 0

        self.min_label = ttk.Label(self, text='min', padding=(0, padding_cheat, 0, 0))
        self.max_label = ttk.Label(self, text='max', padding=(0, padding_cheat, 0, 0))
        self.value_label = ttk.Label(self, width=6, padding=(8, 0))
        self.columnconfigure(col_base + 0, weight=0)
        self.columnconfigure(col_base + 1, weight=0)
        self.columnconfigure(col_base + 2, weight=1)
        self.columnconfigure(col_base + 3, weight=0)
        self.compresser.grid(row=0, column=col_base, sticky=tk.S)
        self.expander.grid(row=0, column=col_base + 1, sticky=tk.S)
        self.slider.grid(row=0, column=col_base + 2, sticky=tk.EW, padx=(6, 0))
        self.value_label.grid(row=0, column=col_base + 3, sticky=tk.W)
        self.min_label.grid(row=1, column=col_base + 2, sticky=tk.NW)
        self.max_label.grid(row=1, column=col_base + 2, sticky=tk.NE)
        self.slider.bind('<MouseWheel>', self._handle_wheel)
        self.slider.bind('<Shift-Button-1>', self._reset)

        self.current_value = left_end

        self.callback = None

        self._update_labels()

    def set_callback(self, callback):
        self.callback = callback

    def set_value(self, value):
        frac = 0.8

        l = self.slider.left_end
        r = self.slider.right_end

        length = r - l
        l_p = frac * l + (1.0 - frac) * r
        r_p = frac * r + (1.0 - frac) * l

        if value < l_p:
            self.slider.left_end  = value - (1.0 - frac) * length
            self.slider.right_end = value + frac * length

        if value > r_p:
            self.slider.left_end  = value - frac * length
            self.slider.right_end = value + (1.0 - frac) * length

        self.slider.set_value(value)
        self.current_value = value
        self._update_labels()

    def _build_icons(self):
        if sys.platform == 'darwin':
            try:
                self.compress_icon=tk.Image('nsimage', source='NSExitFullScreenTemplate',
                                            width=18, height=18)
                self.expand_icon=tk.Image('nsimage', source='NSEnterFullScreenTemplate',
                                            width=18, height=18)
            except tk.TclError:
                self.compress_icon=tk.Image('photo', width=18, height=18,
                    file=os.path.join(os.path.dirname(__file__), 'inward18.png'))

                self.expand_icon=tk.Image('photo', width=18, height=18,
                    file=os.path.join(os.path.dirname(__file__), 'outward18.png'))
        else:
            suffix = 'gif' if tkinter.TkVersion < 8.6 else 'png'
            self.compress_icon=tk.Image('photo', width=18, height=18,
                file=os.path.join(os.path.dirname(__file__), 'inward18.' + suffix))

            self.expand_icon=tk.Image('photo', width=18, height=18,
                file=os.path.join(os.path.dirname(__file__), 'outward18.' + suffix))

    def _slider_callback(self, value):
        self.current_value = value
        self._update_labels()
        if self.callback:
            self.callback(value)
        if self.on_change:
            self.on_change()

    def _update_labels(self):
        l = self.slider.left_end
        r = self.slider.right_end

        num_digits = _num_digits(r - l)

        format_str1 = '%%.%df' % (num_digits + 1)
        format_str2 = '%%.%df' %  num_digits

        self.value_label.configure(text = format_str1 % self.current_value)
        self.min_label.configure(text = format_str2 % l)
        self.max_label.configure(text = format_str2 % r)
        if self.on_change:
            self.on_change()

    def _reset(self, event):

        self.slider.left_end = self.left_end
        self.slider.right_end = self.right_end
        self.set_value(self.current_value)

    def _handle_wheel(self, event):
        delta = event.delta
        if delta < 0:
            self.compresser.state(['pressed'])
            while(delta < 0):
                self.zoom_in()
                delta += 1
            self.after(200, lambda :self.compresser.state(['!pressed']))
        elif delta > 0:
            self.expander.state(['pressed'])
            while(delta > 0):
                self.zoom_out()
                delta -= 1
            self.after(200, lambda :self.expander.state(['!pressed']))

    def zoom_in(self):
        l = self.slider.left_end
        r = self.slider.right_end

        if r - l < self.min_span:
            return

        self.slider.left_end  = 0.25 * (3.0 * l + r)
        self.slider.right_end = 0.25 * (3.0 * r + l)

        self.set_value(self.current_value)
        self._update_labels()

    def zoom_out(self):
        l = self.slider.left_end
        r = self.slider.right_end

        if r - l > self.max_span:
            return

        self.slider.left_end  = 0.5 * (3.0 * l - r)
        self.slider.right_end = 0.5 * (3.0 * r - l)
        self.slider.set_value(self.current_value)
        self._update_labels()

def _num_digits(x):
    r = 1
    while x < 1 and r < 7:
        x *= 10.0
        r += 1
    return r
