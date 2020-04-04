import sys, os, tkinter as tk
from tkinter import ttk
        
class ZoomSlider(ttk.Frame):
    """
    A subclass of the ttk Scale widget which adds buttons that expand or
    contract the range of values by a factor of 2.
    """
    min_span = 0.0005
    max_span = 10000
    
    def __init__(self, master, **kwargs):
        self._build_icons()
        # Strip off kwargs that Frames don't understand
        self.nonneg = kwargs.pop('nonneg', False)
        self.orientation = kwargs.pop('orient', tk.HORIZONTAL)
        self.min = min = kwargs.pop('from_', 0)
        self.max = max = kwargs.pop('to', 100)
        initial_value = kwargs.pop('value', (max + min) / 2)
        self.command = kwargs.pop('command', lambda x: None)
        ttk.Frame.__init__(self, master, **kwargs)
        self.original_range = (min, max)
        self.value = value = tk.DoubleVar(self)
        value.trace('w', self._update)
        self.scale = ttk.Scale(self, from_=min, to=max, variable=value,
                                   orient=self.orientation, command=self.set)
        self.compresser = ttk.Button(self, style='Toolbutton', image=self.compress_icon,
                                    command=self.zoom_in)
        self.expander = ttk.Button(self, style='Toolbutton', image=self.expand_icon,
                                    command=self.zoom_out)
        self.min_label = ttk.Label(self, text='min', padding=(0, -6, 0, 0))
        self.max_label = ttk.Label(self, text='max', padding=(0, -6, 0, 0))
        self.value_label = ttk.Label(self, width=6, padding=(8, 0))
        if self.orientation == tk.HORIZONTAL:
            self.columnconfigure(2, weight=1)
            self.compresser.grid(row=0, column=0, sticky=tk.S)
            self.expander.grid(row=0, column=1, sticky=tk.S)
            self.scale.grid(row=0, column=2, sticky=tk.EW+tk.S, padx=(6, 0))
            self.value_label.grid(row=0, column=3, sticky=tk.W)
            self.min_label.grid(row=1, column=2, sticky=tk.NW)
            self.max_label.grid(row=1, column=2, sticky=tk.NE)
        else:
            raise ValueError('Vertical zoom sliders have not been implemented.')
        self.scale.bind('<MouseWheel>', self._handle_wheel)
        self.scale.bind('<Button-1>', self._handle_mouse)
        self.scale.bind('<Shift-Button-1>', self._reset)
        self.set(initial_value)
        self._update()

    def get(self):
        return self.value.get()

    def set(self, value):
        value = float(value)
        if value > self.max:
            value = self.max
        if value < self.min:
            value = self.min
        self.value.set(value)
        self.command(value)

    def configure(self, **kwargs):
        value = kwargs.get('value', None)
        if value:
            self.set(kwargs['value'])
        command = kwargs.get('command', None)
        if command:
            self.command = command

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
            self.compress_icon=tk.Image('photo', width=18, height=18,
                file=os.path.join(os.path.dirname(__file__), 'inward18.png'))
                                            
            self.expand_icon=tk.Image('photo', width=18, height=18,
                file=os.path.join(os.path.dirname(__file__), 'outward18.png'))
            
    def _update(self, *args):
        self.value_label.configure(text='%.5g'%self.value.get())
        self.min_label.configure(text='%.4f'%self.min)
        self.max_label.configure(text='%.4f'%self.max)

    def _handle_mouse(self, event):
        """
        Override the standard mouse event handler which moves the knob
        by 1.0, regardless of the size of the range, when the mouse
        click is in the trough.  this one implements the standard
        behavior -- moving the knob to the click point.
        """
        part = self.scale.identify(event.x, event.y)
        if (part == 'Scale.trough'):
            if self.orientation == tk.HORIZONTAL:
                fraction = float(event.x) / float(self.scale.winfo_width())
            else:
                fraction = event.y / self.scale.winfo_height()
            self.set(self.min + fraction*(self.max - self.min))
            return 'break'

    def _reset(self, event):
        self.min, self.max = self.original_range
        self.set(self.min + self.max / 2)
        self.scale.configure(from_=self.min, to=self.max)

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

    def _zoom(self, new_span):
        current = self.scale.get()
        self.max = current + new_span / 2
        self.min = current - new_span / 2
        if self.nonneg and self.min < 0:
            self.max -= self.min
            self.min = 0
        self.scale.configure(from_=self.min, to=self.max)
        self._update()

    def zoom_in(self):
        new_span = (self.max - self.min) / 2
        if new_span >= self.min_span:
            self._zoom(new_span)

    def zoom_out(self):
        new_span = (self.max - self.min) * 2
        if new_span <= self.max_span:
            self._zoom(new_span)
