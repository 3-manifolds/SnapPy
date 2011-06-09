cdef extern from "twister_main.h":
    int main(int argc, char** argv)
    void* malloc(size_t size)
    void free(void *mem)

def call_twister(args):
    cdef int argc
    cdef char** argv

    argc = len(args)
    argv = <char**>malloc(argc*sizeof(char*))
    for i, a in enumerate(args):
        argv[i] = a
    main(argc, argv)
    free(argv)

