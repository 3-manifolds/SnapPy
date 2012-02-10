import array

try:
  unichr
  def encode_torsion(divisors):
      return ''.join([unichr(x) for x in divisors]).encode('utf8')
except: # Python3
  def encode_torsion(divisors):
      return ''.join([chr(x) for x in divisors]).encode('utf8')

def decode_torsion(utf8):
    return [ord(x) for x in utf8.decode('utf8')]

def encode_matrices(matrices):
    """
    Convert a list of 2x2 integer matrices into a sequence of bytes.
    """
    # The tricky point here is converting signed integers to bytes.
    return bytes(array.array('b', sum(sum(matrices,[]),[])).tostring())

# NOTE: tostring is deprecated in python3, but for now
# it does the same thing as tobytes.
    
