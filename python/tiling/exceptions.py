class TilingError(RuntimeError):
    pass

class ObjectCloseToCoreCurve(TilingError):
    def __init__(self, obj_name, cusp_index):
        self.obj_name = obj_name
        self.cusp_index = cusp_index
        s = self.obj_name if self.obj_name else "Given geometric object"
        super().__init__(
            "%s is very close to the core curve "
            "of cusp %d and might intersect it." % (s, cusp_index))
