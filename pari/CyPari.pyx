# This is part of the UCS2 hack.

cdef extern from "UCS2hack.h":
    pass

cdef public UCS2_hack (char *string, Py_ssize_t length, char *errors) :   
    return string

# malloc
cdef extern from "stdlib.h":
    ctypedef unsigned long size_t
    void* malloc(size_t size)
    void free(void* mem)

# setjmp (stolen from sage/libs/pari/decl.pxi)
cdef extern from 'setjmp.h':
    struct __jmp_buf_tag:
        pass
    ctypedef __jmp_buf_tag jmp_buf
    int setjmp (jmp_buf __env)
    int longjmp (jmp_buf __env, int val)
    
# PARI declarations

cdef extern from "paritype.h":
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

cdef extern from "parisys.h":
    pass

cdef extern from "parigen.h":   
    ctypedef long* GEN
    extern long typ(GEN x)
    extern void settyp(GEN z, long n)
    extern long evaltype(GEN x)
    extern long lg(GEN x)
    extern long lgefint(GEN x)
    extern long signe(GEN x)
    extern long setlg(GEN x, long length)
    extern long setlgefint(GEN x, long length)

cdef extern from "paricast.h":
    pass

cdef extern from "paristio.h":
    ctypedef long* pari_sp
    extern pari_sp avma    # This is the current position on the PARI stack.  

cdef extern from "paricom.h":
    extern int INIT_DFTm
    extern long ndec2nlong(long x)
    extern long prec2ndec(long x)
    extern long nbits2nlong(long x)

cdef extern from "paridecl.h":
    extern void cgiv(GEN x)
    extern GEN matsnf0(GEN x, long flag)
    extern GEN dbltor(double x)
    extern GEN gp_read_str(char* s)
    extern void gaffect(GEN x, GEN y)
    extern GEN gcopy(GEN x)
    extern char* GENtostr(GEN obj)
    extern void pari_init_opts(size_t parisize,
                               unsigned long maxprime,
                               unsigned long init_opts)

cdef extern from "parierr.h":
    pass

cdef extern from "mini_pariinl.h":
    extern int TWOPOTBYTES_IN_LONG # originally defined in parigen.h
    extern GEN cgetg(long length, long type)
    extern GEN stoi(long x)
    extern long itos(GEN x)

cdef extern from "paripriv.h":
    ctypedef struct pariout_t:
        char format  # e,f,g
        long fieldw  # 0 (ignored) or field width
        long sigd    # -1 (all) or number of significant digits printed
        int sp       # 0 = suppress whitespace from output
        int prettyp  # output style: raw, prettyprint, etc
        int TeXstyle
    ctypedef struct gp_data:
        unsigned long primelimit
        jmp_buf env
        void *hist
        void *pp
        void *path
        pariout_t *fmt
        unsigned long flags, lim_lines
        char *help, *prompt, *prompt_cont
        void *T
    extern void gen_output(GEN x, pariout_t *T)
    extern char* GENtostr0(GEN x, pariout_t *T, void (*do_out)(GEN, pariout_t*))
    extern gp_data GP_DATA # The PARI environment struct.
     

cdef sizeof_pari_long = 2**TWOPOTBYTES_IN_LONG

try:
    int(5).bit_length()
    bit_length = lambda x : x.bit_length()
except AttributeError:
    def bit_length(self):
        s = bin(self)
        s = s.lstrip('-0b')
        return len(s)

cdef class PariEnvironment:
    """
    Pari represents integers and floats in radix 2^N where N can be 32
    or 64.  The number of radix 2^N "digits" is determined when the number
    is created.  The number of decimal digits used when printing a PARI
    number is globally determined.
 
    This object controls both the default number of radix 2^N "digits"
    for PARI numbers as well as the precision used when printing
    numbers.

    set_decimal_precision(n) will set the default radix 2^N precision
    so that at least n decimal digits are available.

    set_display_precision(n) causes all PARI numbers to be printed to
    n significant digits, regardless of the number of radix 2^N "digits"
    stored by PARI.
    """
    cdef pariout_t pariout
    cdef int PARI_precision

    def __cinit__(self):
        self.pariout.format = 'e'
        self.pariout.fieldw = 0
        self.pariout.sigd = -1
        self.pariout.sp = 0
        self.pariout.prettyp = 0
        self.pariout.TeXstyle = 0
        self.PARI_precision = ndec2nlong(38)

    cdef pariout_t* get_pariout(self):
        return &(self.pariout)
    
    def set_decimal_precision(self, digits):
        N = int(digits)
        if N < 1:
            raise ValueError('Negative precision specified.')
        self.PARI_precision = ndec2nlong(N)

    def set_display_precision(self, digits):
        N = int(digits) if digits > 0 else -1
        self.pariout.sigd = N
        
    def pari_precision(self):
        return self.PARI_precision
    
PARI_environment = PariEnvironment()
    
def init_opts(parisize, maxprime):
    pari_init_opts(parisize, maxprime, INIT_DFTm)

# Smith Normal Form
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
    
cdef class PARInumber:
    """
    Base class for PARI arbitrary precision numbers.  For a floating
    point number z, the number of radix 2^N digits used to represent z
    can be specified with the optional precision argument when
    instantiating.  If not specified, the default precision specified
    in the PARI environment is used.  The precision argument is
    ignored by the PARIint subclass.
    """
    cdef GEN value
    cdef int pari_precision
    
    def __cinit__(self, data=None, precision=None):
        if precision:
            self.pari_precision = precision
        else:
            self.pari_precision = PARI_environment.pari_precision()
        self.__alloc__()

    def __alloc__(self):
        # Subclasses overide this method
        pass

    def __dealloc__(self):
        free(self.value)
        
    def __init__(self, data=0, precision=None):
        cdef pari_sp av
        cdef GEN temp
        global avma
        av = avma
        if isinstance(data, PARInumber):
            self.acquire(data)
        else:
            try:
                py_value = self.__py_convert(data)
            except ValueError:
                raise ValueError, 'Invalid initialization string for PARI number.'
            if isinstance(data, str):
                temp = gp_read_str(data)
                gaffect(temp, self.value)
            else:
                self.__pari_convert(py_value)
        avma = av

    def __py_convert(self, data):
        # Subclasses override this method.
        pass

    def __pari_convert(self, python_obj):
        # Subclasses override this method.
        pass
    
    def __check_coersion(self, other):
        #Subclasses override this method.
        pass
        
    def __repr__(self):
        cdef char* c_buf
        c_buf = GENtostr0(self.value,
                          PariEnvironment.get_pariout(PARI_environment),
                          &gen_output)
        result = str(c_buf)
        free(c_buf)
        return result
        
    cdef acquire(self, PARInumber other):
        if self.__check_coersion(other):
            gaffect(other.value, self.value)
        else:
            raise ValueError('Cannot convert those PARI types.')
        
    def precision(self):
        """
        Return the number of 32- or 64-bit words used to represent this
        number.
        """
        return self.pari_precision

    def decimal_precision(self):
        """
        Return the approximate precision of this number in base 10.
        """
        return prec2ndec(self.pari_precision + 2)
    
cdef class PARIint(PARInumber):
    
    def __cinit__(self, data=0, precision=None):
        if isinstance(data, int) or isinstance(data, long):
            self.pari_precision = nbits2nlong(bit_length(data))
        elif isinstance(data, str):
            self.pari_precision = ndec2nlong(len(data))
        elif isinstance(data, PARInumber):
            self.pari_precision = data.precision()
        else:
            raise ValueError('Inappropriate initialization for PARIint.')
        self.__alloc__()
        
    def __alloc__(self):
        self.value = <GEN>malloc((self.pari_precision+2)*sizeof_pari_long)
        settyp(self.value, t_INT)
        setlg(self.value, self.pari_precision+2)

    def __py_convert(self, data):
        # We could do this directly, but for now: int -> string -> GEN
        return str(int(data))

    def __pari_convert(self, python_obj):
        string_ref = python_obj
        gaffect(gp_read_str(string_ref), self.value)

    def __check_coersion(self, other):
        return isinstance(other, PARIint)

cdef class PARIfloat(PARInumber):
    
    def __alloc__(self):
        self.value = <GEN>malloc((self.pari_precision+2)*sizeof_pari_long)
        settyp(self.value, t_REAL)
        setlg(self.value, self.pari_precision+2)

    def __py_convert(self, data):
        return float(data)

    def __pari_convert(self, python_obj):
        gaffect(dbltor(python_obj), self.value)

    def __check_coersion(self, other):
        return isinstance(other, PARIint) or isinstance(other, PARIfloat)
