import cython
import numpy as np
cimport numpy as np

from fortran_wrappers cimport c_diffnw

# TO DO load some example inputs to test the module

@cython.boundscheck(False)
cdef void diffnw(
                int mxncomp_g,
                int nrch_g,
                double[:,::1] z_ar_g
            ):
    
    c_diffnw(
            &mxncomp_g,
            &nrch_g,
            &z_ar_g[0,0]
        )
    
cpdef object compute_diffusive(
                 int mxncomp_g,
                 int nrch_g,
                 double[:,::1] z_ar_g
                ):
    
    diffnw(
         mxncomp_g,
         nrch_g,
         z_ar_g
        )