module var_module

    implicit none
    save

    integer :: ots, option_dsbc
    double precision :: cfl, qus, f, yds, rhs1, rhs2, c11, c12, c21, c22, us, thes
    double precision :: vv
    double precision :: dtini, min_Q

    double precision :: frus2, S_0, y_norm_us, y_crit_us, area_norm, area_crit, minNotSwitchRouting, minNotSwitchRouting2

    !integer, allocatable :: nx1(:), noQSKtable(:)
    integer :: ncomp

    !integer :: ntim, mxncomp
    !integer :: chshp
    !* variables for branching channels
    integer :: nlinks, nupbnds, ndnbnds, NAnum, mxnbrch, nends
    double precision :: theta, thetas, thesinv, alfa2, alfa4, phi
    !double precision :: rhs1, rhs2, c11, c12, c21, c22
    !double precision :: yy, skk, qq  !vv, us,
    double precision :: dxini,lastKnownDiffuDT, skk, minDiffuLm, maxDiffuLm!, tfin ntim
    integer :: applyNaturalSection  ! if 1, then attribute table will be activated, if 0, then rectangular channel will be applied

end module var_module