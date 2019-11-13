try:
    import sage.all
    _within_sage = True
except:
    _within_sage = False

if _within_sage:
    from . import extended
