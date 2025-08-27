def to_byte_str(s):
    if isinstance(s, bytes):
        return s
    return s.encode('utf-8')


