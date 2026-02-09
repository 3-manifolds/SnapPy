class ImageBasedTexture(DataBasedTexture):
    def __init__(self, texture_file):
        w, h, rows, info = png.Reader(texture_file).asRGBA8()
        data = bytearray(4 * w * h)
        for i, row in enumerate(rows):
            data[i * 4 * w : (i + 1) * 4 * w] = row

        super().__init__(w, h, data)
