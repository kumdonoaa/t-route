program hello
    use diffusive

    implicit none

    integer :: mxncomp_g, nrch_g
    double precision :: dtini_g, t0_g, tfin_g, saveInterval_g, saveInterval_ev_g, tfin_qlat_g
    double precision :: dt_ql_g, dt_ub_g, dt_db_g
    integer :: nl_ubcd_ar_g, nts_ql_g, nts_ub_g, nts_db_g, ntss_ev_g  !,ntss_g

    integer :: frnw_col
    double precision :: cfl_g, theta_g, so_llm_g
    integer :: tzeq_flag_g !* 0 for lookup tabale; 1 for using procedures to compute trapz.ch.geo.
    integer :: y_opt_g  !* 1 for normal depth(kinematic); 2 for dept of diffusive wave.

    double precision, dimension(:,:), allocatable :: dfrnw_g

    double precision, dimension(:,:), allocatable :: z_ar_g, bo_ar_g, traps_ar_g, tw_ar_g, twcc_ar_g, So_ar_g
    double precision, dimension(:,:), allocatable :: mann_ar_g, manncc_ar_g, dx_ar_g

    double precision, dimension(:,:,:), allocatable :: qlat_g
    double precision, dimension(:,:), allocatable :: ubcd_g, iniq
    double precision, dimension(:), allocatable :: dbcd_g
    !double precision :: timesDepth_g
    double precision, dimension(:,:,:), allocatable ::  q_ev_g, elv_ev_g
    double precision, dimension(:,:,:), allocatable :: ufhlt_m_g,  ufqlt_m_g, ufhlt_f_g, ufqlt_f_g
    integer :: n, i, j, i1, frj, iseg
    integer :: nhincr_m_g, nhincr_f_g
    double precision :: dmy
    double precision :: timestep_ar_g(10)
    integer :: paradim
    double precision, dimension(:), allocatable :: para_ar_g
    double precision :: dt_qtrib_g
    integer :: nts_qtrib_g
    double precision, dimension(:,:), allocatable :: qtrib_g
    doubleprecision :: wsi, qsi, zsi, tmin

    nrch_g=7 !735 !* TW: 933020089 USGS 02086624
    !nrch_g=29 !* TW:933020161 USGS:02086849
    mxncomp_g=2

    dtini_g= 300.0         !* initial time interval in [sec]
    t0_g= 0.0             !* simulation starting time in [hr]
    !tfin_g= 24.0*(30.0-1.0)+ 22.0 !2.0       !* simulation ending time in [hr] #3670.0
    !tfin_g= 2.0       !* simulation ending time in [hr] #3670.0
    !tfin_g= 408.0 !* TW:933020161 USGS:02086849
    tfin_g= 24 !718.0 !718.0 !* TW: 933020089 USGS 02086624 [hr]
    tfin_qlat_g=24 !718.0 !* TW:933020161 USGS:02086849. ** added as qlat data time period is generally not the same as dsbd data.

    saveInterval_g=dtini_g    !* output publishing time interval in [sec]
    cfl_g= 0.99            !* max. allowable Courant number
    dt_ql_g= 3600.0  !* instead of dtinput #time interval of qlateral input data from wrf-hydro [sec]
    dt_ub_g= dt_ql_g  !* time interval of input data for upstream boundary condition  [sec]
    dt_db_g=  900.0 !* time interval of input data for downstream boundary condition  [sec]
    saveInterval_ev_g= 12*dtini_g !saveInterval_g !for evaluation later [sec]

    !nl_ubcd_ar_g=266749 !* row number of ubcd_ar.txt TW: 933020089 USGS 02086624
    nl_ubcd_ar_g= 336 !267468
    ntss_ev_g= int((tfin_g - t0_g)*3600.0/dtini_g, KIND(ntss_ev_g))+1
    nts_ql_g= int( (tfin_qlat_g - t0_g)*3600.0/dt_ql_g, KIND(nts_ql_g) ) !* based on time of [sec]
    !nts_ql_g= int( (tfin_qlat_g - t0_g)*3600.0/dtini_g+1, KIND(nts_ql_g) ) !* based on time of [sec]
    nts_ub_g= nts_ql_g
    nts_db_g= int( (tfin_g - t0_g)*3600.0/dt_db_g, KIND(nts_db_g) ) !* based on time of [sec]

!    minq_ns_g=0.1 !* absolute minimum discharge (lower limit) for numerical stability [m^3/sec]
!    mindepth_ns_g=0.1 !* minimum water depth(lower limit) for numerical stability [meter]
    so_llm_g=0.00001   !* lower limit of channel bed slope
    frnw_col=8

    !* MC result to tributary (including mainstem upstream boundary) flows
    dt_qtrib_g=60.0 ![sec]
    nts_qtrib_g = int((tfin_g - t0_g) * 3600.0 / dt_qtrib_g)


    open(unit=21, file="./input/frnw_ar.txt")
    !open(unit=22, file="./input/z_bo_traps_tw_twcc_m_mcc_so_dx.txt")
    open(unit=22, file="./input/z_ar.txt")
    open(unit=23, file="./input/bo_ar.txt")
    open(unit=24, file="./input/traps_ar.txt")
    open(unit=25, file="./input/tw_ar.txt")
    open(unit=26, file="./input/twcc_ar.txt")
    open(unit=27, file="./input/mann_ar.txt")
    open(unit=28, file="./input/manncc_ar.txt")
    open(unit=29, file="./input/so_ar.txt")
    open(unit=30, file="./input/dx_ar.txt")
    open(unit=31, file="./input/qlat_ar.txt")
    !open(unit=32, file="./input/ubcd_ar.txt")
    !open(unit=33, file="./input/dbcd_ar.txt")
    !open(unit=33, file="./input/dbcd_ar_datumCorrected.txt")
    open(unit=34, file="./input/qtrib_ar.txt")

    allocate(dfrnw_g(nrch_g, frnw_col))
    do j=1, nrch_g
        read(21,*) (dfrnw_g(j,i),i=1, frnw_col)
    enddo

    allocate( z_ar_g(mxncomp_g, nrch_g), bo_ar_g(mxncomp_g, nrch_g), traps_ar_g(mxncomp_g, nrch_g) )
    allocate( tw_ar_g(mxncomp_g, nrch_g), twcc_ar_g(mxncomp_g, nrch_g) )
    allocate( mann_ar_g(mxncomp_g, nrch_g), manncc_ar_g(mxncomp_g, nrch_g) )
    allocate( so_ar_g(mxncomp_g, nrch_g), dx_ar_g(mxncomp_g, nrch_g) )
    allocate( qlat_g(nts_ql_g, mxncomp_g, nrch_g) )
    allocate( iniq(mxncomp_g, nrch_g) )


    iniq=0.5

    do i=1, mxncomp_g
        read(22,*) (z_ar_g(i,j),j=1,nrch_g)
    enddo
    do i=1, mxncomp_g
        read(23,*) (bo_ar_g(i,j),j=1,nrch_g)
    enddo
    do i=1, mxncomp_g
        read(24,*) (traps_ar_g(i,j),j=1,nrch_g)
    enddo
    do i=1, mxncomp_g
        read(25,*) (tw_ar_g(i,j),j=1,nrch_g)
    enddo
    do i=1, mxncomp_g
        read(26,*) (twcc_ar_g(i,j),j=1,nrch_g)
    enddo
    do i=1, mxncomp_g
        read(27,*) (mann_ar_g(i,j),j=1,nrch_g)
    enddo
    do i=1, mxncomp_g
        read(28,*) (manncc_ar_g(i,j),j=1,nrch_g)
    enddo
    do i=1, mxncomp_g
        read(29,*) (so_ar_g(i,j),j=1,nrch_g)
    enddo
    do i=1, mxncomp_g
        read(30,*) (dx_ar_g(i,j),j=1,nrch_g)
    enddo
    !++------------------------------------------------------+
    !+              Lateral inflow
    !+
    !+   IMPORTANT: it should have a unit of m^2/sec
    !++------------------------------------------------------+
    do j=1, nrch_g
        ncomp=int(dfrnw_g(j,1))
        do i=1, ncomp
            do n=1, nts_ql_g
                !* tarr_ql(n) in [min]
                read(31,*) frj, iseg, dmy, qlat_g(n,i,j) !* qlat_g(n,i,j) in [m^2/sec]
                qlat_g(n,i,j)=0.0
            end do
        end do
    enddo




    !++--------------------------------------------------------------+
    !+              Upstream boundary condition
    !+
    !+ Assume that hydrographs are always available at the upper ends
    !+ of the most upstream links. m^3/sec
    !++--------------------------------------------------------------+
    allocate( ubcd_g(nts_ub_g, nrch_g) )
    do j=1, nrch_g
        do n=1,nts_ub_g
            ubcd_g(n,j)= 0.0 !* initialize
        end do
    end do
!    !* ubcd_g in [m^3/sec]
!    do i1=1, nl_ubcd_ar_g
!        read(32,*) frj, n, dmy, ubcd_g(n,frj)  !* time of tarr_ub in [min]
!    enddo
    !++----------------------------------------------------------------+
    !+          Stage data at the downstream end
    !+
    !+
    !+ when flow conditions of lower end of link j are known.
    !+ Assume that stage time-series data is always available at the
    !+ lower end of the most downstream link.
    !++----------------------------------------------------------------+
    allocate( dbcd_g(nts_db_g) )
    dbcd_g=0.0 !* not used for now.
!    allocate( dbcd_g(nts_db_g) )
!    do n=1, nts_db_g
!        !* time of tarr_db in [min]
!        read(33,*) dmy, dbcd_g(n) !, dbcd_g(n,2) !* dbcd_g(n,1): elevation in [meter] based on NAD88 meter datum that
!                                                            !*               RouteLink uses OR stage data [meter] if available
!                                                            !* dbcd_g(n,2): discharge data [cms] if available.
!    end do
    !++---------------------------------------------------------------------------+
    !+    MC results to tributary(including mainstem upstream boundary) flow data
    !+
    !++---------------------------------------------------------------------------+
    allocate (qtrib_g(nts_qtrib_g, nrch_g))
    do j=1, nrch_g
        do n=1, nts_qtrib_g
            !read(34,*) frj, dmy, qtrib_g(n,j) !*  [m^3/sec]
            qtrib_g(n,j)=0.0
        end do
    enddo
    !* mainstem upstream boundary data
    qsi= 100.0
    wsi= 80.0
    zsi= 0.6
    do n=1, nts_qtrib_g
        tmin= dt_qtrib_g/60.0*real(n-1) !*[min]
        if (n==1) then
            tmin=0.1
        end if
        dmy= qsi*exp(zsi*(2.0 - tmin/wsi - wsi/tmin))/((tmin/wsi)**(3.0/2.0))
        qtrib_g(n,2)= dmy
    enddo



    nhincr_m_g=2
    nhincr_f_g=2
    allocate(ufhlt_m_g(mxncomp_g,nrch_g,nhincr_m_g),  ufqlt_m_g(mxncomp_g,nrch_g,nhincr_m_g))
    allocate(ufhlt_f_g(mxncomp_g,nrch_g,nhincr_f_g),  ufqlt_f_g(mxncomp_g,nrch_g,nhincr_f_g))



    allocate(q_ev_g(ntss_ev_g, mxncomp_g, nrch_g), elv_ev_g(ntss_ev_g, mxncomp_g, nrch_g))


!    call diffnw(dtini_g, t0_g, tfin_g, saveinterval_g, saveinterval_ev_g, dt_ql_g, dt_ub_g, dt_db_g, &
!                nts_ql_g, nts_ub_g, nts_db_g, &
!                mxncomp_g, nrch_g, z_ar_g, bo_ar_g, traps_ar_g, tw_ar_g, twcc_ar_g, &
!                mann_ar_g, manncc_ar_g, So_ar_g, dx_ar_g, &
!                nhincr_m_g, nhincr_f_g, ufhlt_m_g,  ufqlt_m_g, ufhlt_f_g, ufqlt_f_g, &
!                frnw_col, dfrnw_g, qlat_g, ubcd_g, dbcd_g, &
!                cfl_g, theta_g, tzeq_flag_g, y_opt_g, so_llm_g, &
!                ntss_ev_g, q_ev_g, elv_ev_g)
!    call diffnw(dtini_g, t0_g, tfin_g, saveinterval_ev_g, dt_ql_g, dt_ub_g, dt_db_g, &
!                        nts_ql_g, nts_ub_g, nts_db_g, &
!                        mxncomp_g, nrch_g, z_ar_g, bo_ar_g, traps_ar_g, tw_ar_g, twcc_ar_g, &
!                        mann_ar_g, manncc_ar_g, so_ar_g, dx_ar_g, iniq, &
!                        nhincr_m_g, nhincr_f_g, ufhlt_m_g,  ufqlt_m_g, ufhlt_f_g, ufqlt_f_g, &
!                        frnw_col, dfrnw_g, qlat_g, ubcd_g, dbcd_g, &
!                        cfl_g, theta_g, tzeq_flag_g, y_opt_g, so_llm_g, &
!                        ntss_ev_g, q_ev_g, elv_ev_g)

        !*---------------------
        !* time step variables
        !*---------------------
        timestep_ar_g(1) = dtini_g  !*[sec]
        timestep_ar_g(2) = t0_g ![hr]
        timestep_ar_g(3) = tfin_g
        timestep_ar_g(4) = saveinterval_ev_g  !*[sec]
        timestep_ar_g(5) = dt_ql_g
        timestep_ar_g(6) = dt_ub_g
        timestep_ar_g(7) = dt_db_g
        timestep_ar_g(8) = dt_qtrib_g
        !*----------------------------
        !* sensitive model parameters
        !*----------------------------
        paradim=9
        allocate(para_ar_g(paradim))
        para_ar_g(1) = 0.5 !, lower limit of celerity
        para_ar_g(2) = 50.0 !, lower limit of diffusivity
        para_ar_g(3) = 400.0 !, upper limit of diffusivity
        para_ar_g(4) = 0.1  !, lower limit of dimensionless space step, used to reduce ocillation of hydrograph
        para_ar_g(5) = 0.3 !, lower limit of dimensionless time step, used to reduce ociallation of hydrograph
        para_ar_g(6) = 2  !* the number of ghost time node, used to get adequate dischage boundary condition in time.
        para_ar_g(7) =	0.02831 !. lower limit of discharge
        para_ar_g(8) = 0.0001 !, lower limit of channel bed slope
        para_ar_g(9) =  3 !* 1: bisection, 2: simple iterative with contribution from dU/dX, 3: Newton Raphson for water depth

        call diffnw(timestep_ar_g, nts_ql_g, nts_ub_g, nts_db_g, ntss_ev_g, nts_qtrib_g, &
                    mxncomp_g, nrch_g, z_ar_g, bo_ar_g, traps_ar_g, tw_ar_g, twcc_ar_g, &
                    mann_ar_g, manncc_ar_g, so_ar_g, dx_ar_g, iniq, &
                    frnw_col, dfrnw_g, qlat_g, ubcd_g, dbcd_g, qtrib_g, &
                    paradim, para_ar_g, q_ev_g, elv_ev_g)



end program
