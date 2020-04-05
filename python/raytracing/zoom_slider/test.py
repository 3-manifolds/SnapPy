from . import tk, Slider, ZoomSlider

def test_widget(window=None):
    if window is None:
        window = tk.Toplevel()
    window.title('ZoomSliders')
    tk.Label(window, text='', width=70).pack()
    slider = Slider(window, from_=0.0, to=1.0)
    slider.pack(expand=True, fill=tk.X, padx=10)
    zoom_slider = ZoomSlider(window, from_=0.0, to=1.0)
    zoom_slider.pack(expand=True, fill=tk.X)
    pos_slider = ZoomSlider(window, from_=0.0, to=1.0, nonneg=True)
    pos_slider.pack(expand=True, fill=tk.X, side='bottom')

if __name__ == '__main__':
    win = tk.Tk()
    test_widget(win)
    win.mainloop()
