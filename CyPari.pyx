# PARI declarations

cdef extern from "pari.h":
     cdef enum:
         t_INT    =  1
         t_REAL   =  2
         t_INTMOD =  3
         t_FRAC   =  4
         t_COMPLEX=  6
         t_PADIC  =  7
         t_QUAD   =  8
         t_POLMOD =  9
         t_POL    =  10
         t_SER    =  11
         t_RFRAC  =  13
         t_QFR    =  15
         t_QFI    =  16
         t_VEC    =  17
         t_COL    =  18
         t_MAT    =  19
         t_LIST   =  20
         t_STR    =  21
         t_VECSMALL= 22

     ctypedef long* GEN
     ctypedef long* pari_sp
     extern void cgiv(GEN x)
     extern GEN cgetg(long length, long type)
     extern GEN matsnf0(GEN x, long flag)
     extern GEN stoi(long x)
     extern long itos(GEN x)
     extern long lg(GEN x)
     extern long signe(GEN x)
     extern void pari_init_opts(size_t parisize, unsigned long maxprime, unsigned long init_opts)
     extern pari_sp avma   # This is the current position on the stack.  

def init_opts(parisize, maxprime, init_opts):
     pari_init_opts(parisize, maxprime, init_opts)
     
def smith_form(M):
    cdef GEN pari_matrix
    cdef GEN pari_vector
    cdef GEN pari_int
    cdef int i, j
    cdef pari_sp av
    global avma
    
    try:
        m, n = M.shape
    except AttributeError:
        try:
             m, n = M.nrows(), M.ncols()
        except AttributeError:
             try:
                  m, n = len(M), len(M[0])
             except:
                  raise ValueError, "This is not a recognized matrix type"
             
    # Memory management in PARI is very primitive.  Here, once we've
    # copied the answer into Python, we don't care about anything PARI
    # created, so we record the current position on the stack and then
    # restore it at the end.
    av = avma 

    # Copy from Python into PARI
    
    pari_matrix = cgetg(n+1, t_MAT)
    for j from 1 <= j <= n:
        pari_matrix[j] = <long>cgetg(m+1, t_COL)
    for i from 1 <= i <= m:
        for j from 1 <= j <= n:
             try:
                  (<GEN*>pari_matrix)[j][i] =  <long>stoi(M[i-1,j-1])
             except TypeError:
                  (<GEN*>pari_matrix)[j][i] =  <long>stoi(M[i-1][j-1])

    # Do the computation
    pari_vector = matsnf0(pari_matrix, 4)

    # Extract the result
    
    result = []
    for i from 1 <= i < lg(pari_vector):
        pari_int = (<GEN*>pari_vector)[i]
        result.append(itos(pari_int))

    # Restore the stack position, trashing all PARI computations that
    # this function did.

    avma = av
    return result
