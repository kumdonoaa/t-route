module var

    implicit none

    integer :: nts_ql_g, nts_ub_g, nts_db_g,  nts_qtrib_g
    integer :: frnw_col
    double precision, dimension(:), allocatable :: timestep_ar_g
    double precision, dimension(:,:), allocatable :: z_ar_g, bo_ar_g, traps_ar_g, tw_ar_g, twcc_ar_g
    double precision, dimension(:,:), allocatable :: mann_ar_g, manncc_ar_g, dx_ar_g, iniq
    double precision, dimension(:,:), allocatable :: dfrnw_g
    double precision, dimension(:,:,:), allocatable :: qlat_g
    double precision, dimension(:,:), allocatable :: ubcd_g
    double precision, dimension(:), allocatable :: dbcd_g
    double precision, dimension(:,:), allocatable :: qtrib_g
    integer :: paradim
    double precision, dimension(:), allocatable  :: para_ar_g
    double precision, dimension(:,:), allocatable  :: so_ar_g
endmodule var
