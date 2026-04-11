def to_bytes(s) -> bytes:
    if isinstance(s, bytes):
        return s
    return s.encode()

