module mdiffv2
contains
    
    subroutine diffnw(mxncomp_g, nrch_g, z_ar_g)
    
        implicit none
        integer, intent(in) :: mxncomp_g, nrch_g
        double precision, dimension(mxncomp_g, nrch_g), intent(in) :: z_ar_g

        print *, z_ar_g
    
    endsubroutine diffnw
endmodule mdiffv2