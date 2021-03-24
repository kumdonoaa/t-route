module mdiffv2_interface

use, intrinsic :: iso_c_binding
use mdiffv2, only: diffnw

implicit none
contains
subroutine c_diffnw(mxncomp_g, nrch_g, z_ar_g) bind(c)

    integer(c_int), intent(in) :: mxncomp_g, nrch_g
    real(c_double), dimension(mxncomp_g, nrch_g), intent(in) :: z_ar_g

    call diffnw(mxncomp_g, nrch_g, z_ar_g)
    
end subroutine c_diffnw
end module mdiffv2_interface