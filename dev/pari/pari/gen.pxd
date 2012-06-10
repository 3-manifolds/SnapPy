include 'decl.pxi'

#cimport sage.structure.element 

#cdef class gen(sage.structure.element.RingElement):
cdef class gen:
    cdef GEN g
    cdef object _refers_to
    cdef pari_sp b
    cdef void init(self, GEN g, pari_sp b)
    cdef GEN _gen(self)
    cdef gen new_gen(self, GEN x)
    cdef gen new_gen_noclear(self, GEN x)
    cdef gen add_gens(self, gen right)
    cdef gen pari(self, object x)
    cdef GEN _deepcopy_to_python_heap(self, GEN x, pari_sp* address)
    cdef long get_var(self, v)
    cdef GEN get_nf(self) except NULL

#x#cimport sage.structure.parent_base

#cdef class PariInstance(sage.structure.parent_base.ParentWithBase):
cdef class PariInstance:
    cdef gen PARI_ZERO, PARI_ONE, PARI_TWO
    cdef gen new_gen(self, GEN x)
    cdef gen new_leaf_gen(self, GEN x)
    cdef object new_gen_to_string(self, GEN x)
    cdef gen new_gen_noclear(self, GEN x)
#x#    cdef gen new_gen_from_mpz_t(self, mpz_t value)
#x#    cdef inline GEN _new_GEN_from_mpz_t(self, mpz_t value)
#x#    cdef gen new_gen_from_mpq_t(self, mpq_t value)
#x#    cdef inline GEN _new_GEN_from_mpq_t(self, mpq_t value)
    cdef gen new_gen_from_int(self, int value)
    cdef gen new_t_POL_from_int_star(self, int *vals, int length, long varnum)
#x#    cdef gen new_gen_from_padic(self, long ordp, long relprec, mpz_t prime, mpz_t p_pow, mpz_t unit)
    cdef void clear_stack(self)
    cdef void set_mytop_to_avma(self)
    cdef gen double_to_gen_c(self, double)
    cdef GEN double_to_GEN(self, double)
    cdef GEN deepcopy_to_python_heap(self, GEN x, pari_sp* address)
    cdef gen new_ref(self, GEN g, gen parent)
    cdef gen _empty_vector(self, long n)
    cdef long get_var(self, v)
    cdef GEN toGEN(self, x, int i) except NULL
#x#    cdef GEN _new_GEN_from_mpz_t_matrix(self, mpz_t** B, Py_ssize_t nr, Py_ssize_t nc)
#x#    cdef GEN _new_GEN_from_mpz_t_matrix_rotate90(self, mpz_t** B, Py_ssize_t nr, Py_ssize_t nc)
#x#    cdef gen integer_matrix(self, mpz_t** B, Py_ssize_t nr, Py_ssize_t nc, bint permute_for_hnf)
#x#    cdef GEN _new_GEN_from_mpq_t_matrix(self, mpq_t** B, Py_ssize_t nr, Py_ssize_t nc)
#x#    cdef gen rational_matrix(self, mpq_t** B, Py_ssize_t nr, Py_ssize_t nc)

cdef GEN _Vec_append(GEN v, GEN a, long n)
