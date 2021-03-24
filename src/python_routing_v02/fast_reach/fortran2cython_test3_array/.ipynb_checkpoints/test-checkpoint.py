import numpy as np
import diffusive

print("this is only a test")

mxncomp_g = 4
nrch_g = 6
z_ar_g = np.empty([mxncomp_g,nrch_g], dtype = np.double)

print("Hello from Python:")
print(z_ar_g)


print("Hello from Fortran")
A = diffusive.compute_diffusive(
                                mxncomp_g,
                                nrch_g,
                                z_ar_g,
                               )