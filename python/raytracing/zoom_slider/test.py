from . import tk, ZoomSlider

def test_widget(window=None):
    if window is None:
        window = tk.Toplevel()
    window.title('ZoomSliders')
    tk.Label(window, text='', width=70).pack()
    slider = ZoomSlider(window, from_=0.0, to=1.0)
    slider.pack(expand=True, fill=tk.X)
    pos_slider = ZoomSlider(window, from_=0.0, to=1.0, nonneg=True)
    pos_slider.pack(expand=True, fill=tk.X, side='bottom')

if __name__ == '__main__':
    win = tk.Tk()
    test_widget(win)
    win.mainloop()
