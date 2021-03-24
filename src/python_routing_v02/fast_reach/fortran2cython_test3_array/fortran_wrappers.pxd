cdef extern from "pydiffusive.h":
    void c_diffnw(int *mxncomp_g,
                 int *nrch_g,
                 double *z_ar_g) nogil;
