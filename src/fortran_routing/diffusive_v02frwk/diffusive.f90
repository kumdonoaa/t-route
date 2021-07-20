module diffusive
    !*-----------------------------------------------------------------------------------------------------
    !*      This diffusive model is developed by Dong Ha Kim at the National Water Center by numerically
    !*      solving a diffusive equation using Crank-Nicholson scheme (CNX) as described in a paper titled
    !*      "Algorithms for solving the diffusive wave flood routing equation" by Moussa R. and et al.
    !*      As described in the paper, a technique called resolution technique solves the Crank-Nicholson
    !*      equations. This C-N scheme requires uniform space steps. Subsequently, each and every reach
    !*      computes a uniform segment length that is actually used in discharge computation. Any variables
    !*      here with "_m_" handles model variables operating in the uniform space domain. The CNX is
    !*      unconditionally stable thus a fixed simulation time step (dtini) is used in this model.
    !*      After computing discharge values by the resolution technique at each simulation time
    !*      step, water depth is accordingly computed always by normal depth.
    !*------------------------------------------------------------------------------------------------------
    implicit none

    double precision, parameter :: grav = 9.81
    double precision, parameter :: TOLERANCE = 1e-8
    integer :: nrch
    doubleprecision :: dtini
    double precision, dimension(:,:), allocatable :: q_g, elv_g, qpx_g, q_pre_g, elv_pre_g
    double precision, dimension(:), allocatable :: tarr_ql, tarr_ub, tarr_db, varr_ql, varr_ub, varr_db
    integer, dimension(:,:), allocatable :: frnw_g
    double precision, dimension(:), allocatable :: qlatj_g
    double precision, dimension(:), allocatable :: qlatj_m_g
    double precision :: z_g, bo_g, traps_g, tw_g, twcc_g, so_g, mann_g, manncc_g
    double precision, dimension(:,:), allocatable :: adjz_ar_g, adjso_ar_g
    double precision, dimension(:,:), allocatable :: dx_ar, bo_ar, traps_ar, tw_ar, twcc_ar, mann_ar, manncc_ar
    integer, dimension(:,:), allocatable :: lcvQ_g, perQmap_g
    double precision, dimension(:), allocatable :: pqrapj_g
    double precision, dimension(:,:), allocatable :: c_g, d_g
    integer :: perdim_g,  mxncomp_m_g
    double precision, dimension(:,:), allocatable ::  qnm1_g,  elvnm1_g, qlatn_g
    integer, dimension(:), allocatable ::  ncomp_m_g
    double precision, dimension(:,:), allocatable :: bo_ar_m_g, traps_ar_m_g, tw_ar_m_g, twcc_ar_m_g
    double precision, dimension(:,:), allocatable :: mann_ar_m_g, manncc_ar_m_g, adjz_ar_m_g, dx_ar_m_g, adjso_ar_m_g
    double precision, dimension(:,:), allocatable :: q_m_g, qn_m_g
    double precision, dimension(:,:), allocatable :: elv_m_g, elvn_m_g, elvnm1_m_g
    double precision, dimension(:,:), allocatable ::  qlatn_m_g
    double precision :: minq_ns_g, mindepth_ns_g  !* minimum Q and depth for numerical stability
    double precision, dimension(:,:), allocatable :: refQ_g  !* reference discharge for computing C and D
    double precision, dimension(:,:), allocatable :: q_drainql_g, mxq_drainql_g
    double precision :: dmytime_g, dmydq_nm_g, dmyq_nm_g, dmyqn_nm_g, dmymxqij_g
    integer, dimension(:,:), allocatable :: fail_flag_g
    integer :: dim_nmLT_g
    double precision, dimension(:,:,:,:), allocatable :: normdepthLT_g, normdepthLT_m_g
    double precision, dimension(:,:), allocatable :: p, qq, r, w

contains

    subroutine diffnw(dtini_g, t0_g, tfin_g, saveinterval_g, saveinterval_ev_g, dt_ql_g, dt_ub_g, dt_db_g, &
                        nts_ql_g, nts_ub_g, nts_db_g, &
                        mxncomp_g, nrch_g, z_ar_g, bo_ar_g, traps_ar_g, tw_ar_g, twcc_ar_g, &
                        mann_ar_g, manncc_ar_g, so_ar_g, dx_ar_g, &
                        nhincr_m_g, nhincr_f_g, ufhlt_m_g,  ufqlt_m_g, ufhlt_f_g, ufqlt_f_g, &
                        frnw_col, dfrnw_g, qlat_g, ubcd_g, dbcd_g, &
                        cfl_g, theta_g, tzeq_flag_g, y_opt_g, so_llm_g, &
                        ntss_ev_g, q_ev_g, elv_ev_g)

        implicit none

        integer, intent(in) :: mxncomp_g, nrch_g
        integer, intent(in) :: nts_ql_g, nts_ub_g, nts_db_g, ntss_ev_g
        integer, intent(in) :: frnw_col
        integer, intent(in) :: nhincr_m_g, nhincr_f_g
        double precision,intent(in) :: dtini_g, t0_g, tfin_g, saveinterval_g, saveinterval_ev_g, dt_ql_g, dt_ub_g, dt_db_g
        double precision, dimension(mxncomp_g, nrch_g), intent(in) :: z_ar_g, bo_ar_g, traps_ar_g, tw_ar_g, twcc_ar_g
        double precision, dimension(mxncomp_g, nrch_g), intent(in) :: mann_ar_g, manncc_ar_g, dx_ar_g
        double precision, dimension(mxncomp_g, nrch_g, nhincr_m_g), intent(in) :: ufhlt_m_g,  ufqlt_m_g
        double precision, dimension(mxncomp_g, nrch_g, nhincr_f_g), intent(in) :: ufhlt_f_g, ufqlt_f_g
        double precision, dimension(nrch_g, frnw_col), intent(in) :: dfrnw_g
        double precision, dimension(nts_ql_g, mxncomp_g, nrch_g), intent(in) :: qlat_g
        double precision, dimension(nts_ub_g, nrch_g), intent(in) :: ubcd_g
        double precision, dimension(nts_db_g), intent(in) :: dbcd_g
        double precision, intent(in) :: cfl_g, theta_g, so_llm_g
        integer, intent(in) :: tzeq_flag_g !* 0 for lookup tabale; 1 for using procedures to compute trapz.ch.geo.
        integer, intent(in) :: y_opt_g  !* 1 for normal depth(kinematic); 2 for dept of diffusive wave.
        double precision, dimension(mxncomp_g, nrch_g), intent(inout) :: so_ar_g
        double precision, dimension(ntss_ev_g, mxncomp_g, nrch_g), intent(out) :: q_ev_g, elv_ev_g
        integer :: ncomp
        integer :: ts, n, i, j, ts_ev
        double precision ::  tc, tf0
        integer  :: nusrch, rch, usrchj, ncomp_usrchj, dsrchj
        double precision :: qjt
        double precision ::   saveInterval_min,  rmnd
        integer ::  ncomp_m, i_m
        double precision ::  initval
        double precision :: hbf, z_tw
        double precision, dimension(:,:,:), allocatable :: q_drainql
        double precision ::  sumdx
        integer :: nseg_m, ncomp_m_usrchj
        double precision, dimension(:,:), allocatable :: dmy2d
        double precision, dimension(:), allocatable :: dmy1d
        double precision :: dmy1, dmy2, dmy3, dxst
        integer :: dimLT, ilt, ict, ict_1, ict_2, ict_3
        double precision :: step, qnorm, h_m, h_f
        integer :: twj
        double precision :: Qxs, y_norm

        nrch=nrch_g
        dtini= dtini_g
        allocate(dx_ar(mxncomp_g, nrch_g), bo_ar(mxncomp_g, nrch_g), traps_ar(mxncomp_g, nrch_g))
        allocate(tw_ar(mxncomp_g, nrch_g), twcc_ar(mxncomp_g, nrch_g))
        allocate(mann_ar(mxncomp_g, nrch_g), manncc_ar(mxncomp_g, nrch_g))
        dx_ar= dx_ar_g
        bo_ar= bo_ar_g
        traps_ar= traps_ar_g
        tw_ar= tw_ar_g
        twcc_ar= twcc_ar_g
        mann_ar= mann_ar_g
        manncc_ar= manncc_ar_g
        allocate(frnw_g(nrch_g,frnw_col))
        frnw_g=int(dfrnw_g)
        minq_ns_g=0.1 !* absolute minimum discharge (lower limit) for numerical stability [m^3/sec]
        mindepth_ns_g=0.1 !* minimum water depth(lower limit) for numerical stability [meter]

        allocate(q_g(mxncomp_g, nrch_g), elv_g(mxncomp_g, nrch_g))
        allocate(q_pre_g(mxncomp_g, nrch_g), elv_pre_g(mxncomp_g, nrch_g))
        allocate(tarr_ql(nts_ql_g), varr_ql(nts_ql_g))
        allocate(tarr_ub(nts_ub_g), varr_ub(nts_ub_g))
        allocate(tarr_db(nts_db_g), varr_db(nts_db_g))

        !open(unit=2, file="./forran diffusive discharge water depth.txt")
        !write(2,"(3A10, 4A20)") "ts", "i", "j", "q_g", "depth", "C", "D"

        !* time step series for lateral flow
        do n=1, nts_ql_g
            tarr_ql(n)= t0_g*60.0 + dt_ql_g*real(n-1,KIND(dt_ql_g))/60.0 !* [min]
        end do
        !* time step series for upstream boundary data
        do n=1, nts_ub_g
            tarr_ub(n)= t0_g*60.0 + dt_ub_g*real(n-1,KIND(dt_ub_g))/60.0 !* [min]
        enddo
        !* when measured data is used, time step series for downstream boundary data
        do n=1, nts_db_g
            tarr_db(n)= t0_g*60.0 + dt_db_g*real(n-1,KIND(dt_db_g))/60.0 !* [min]
        enddo
        !**-----------------------------------------------------------------------------------*
        !*                  Adjust abnormally small channel bottom slope
        !*
        !* Definition of abnormal slope:
        !*  1. slopes less than a chosen lower limit (=so_llm)
        !* Repair method:
        !*  1. take average of slopes of adjacent segments
        !**-----------------------------------------------------------------------------------*
        allocate(adjso_ar_g(mxncomp_g,nrch_g))
        adjso_ar_g= so_ar_g
        !** NOTE:
        !* so_flag=1: activate the adjustment of bottom slope as below
        !* so_flag=0: Use slope as it is
!        so_flag=0
!        if (so_flag==1) then
!            do j=1,nrch_g
!                ncomp=frnw_g(j,1)
!                do i=1, ncomp-1
!                    if (so_ar_g(i,j).lt.so_llm_g) then
!                        !* adjacent upstream segment's slope
!                        so_us=0.0
!                        c_us=0.0
!                        i2=i-1
!                        do while (i2.ge.1)
!                            if (so_ar_g(i2,j).ge.so_llm_g) then
!                                so_us= so_ar_g(i2,j)
!                                c_us=1.0
!                                exit
!                            endif
!                            i2=i2-1
!                        end do
!                        !* adjacent downstream segment's slope
!                        so_ds=0.0
!                        c_ds=0.0
!                        i2=i+1
!                        do while (i2.le.ncomp-1)
!                            if (so_ar_g(i2,j).ge.so_llm_g) then
!                                so_ds= so_ar_g(i2,j)
!                                c_ds=1.0
!                                exit
!                            endif
!                            i2=i2+1
!                        end do
!                        if (c_us+c_ds.gt.0.0) then
!                            adjso= (so_us+so_ds)/(c_us+c_ds)
!                            adjso_ar_g(i,j)= adjso
!                        else
!                            adjso_ar_g(i,j)= so_llm_g
!                        endif
!                    endif
!                enddo
!            enddo
!        endif
        !++-----------------------------------------------------------------------------------+
        !+        adjust altitude one more time for counting for z_us = z_ds + so*dx
        !+
        !+        Refer to p48, RM5
        !++-----------------------------------------------------------------------------------+
        allocate(adjz_ar_g(mxncomp_g, nrch_g))
        do j=nrch_g, 1, -1
            if (frnw_g(j,2)<0.0) then
                twj= j !* TW reach index
                exit
            endif
        enddo
        ncomp=frnw_g(twj,1)
        z_tw= z_ar_g(ncomp, twj)
        do j=nrch_g, 1, -1
            ncomp=frnw_g(j,1)
            if (frnw_g(j,2)<0.0) then
            !* downstream boundary node at TW
                adjz_ar_g(ncomp, j)= z_tw
            else
            !* downstream boundary node at a junction
                dsrchj= frnw_g(j,2)    !* j index of the downstream reach
                adjz_ar_g(ncomp, j)= adjz_ar_g(1, dsrchj)
            endif
            !* For the following nodes within a reach before the bottom node:
            do i=ncomp-1,1,-1
                adjz_ar_g(i, j)= adjz_ar_g(i+1, j) + adjso_ar_g(i,j)*dx_ar_g(i,j)
            end do
        enddo
        !**----------------------------------------------------------------------*
        !*         Uniform flow lookup tables on original space domain
        !*
        !**----------------------------------------------------------------------*
        dim_nmLT_g=80
        dimLT=dim_nmLT_g
        allocate (normdepthLT_g(mxncomp_g,nrch_g,dimLT,2))
        do j=1, nrch_g
            do i=1, frnw_g(j,1)
                bo_g= bo_ar_g(i,j)
                traps_g= traps_ar_g(i,j)
                tw_g= tw_ar_g(i,j)
                twcc_g= twcc_ar_g(i,j)
                so_g= adjso_ar_g(i,j)
                mann_g= mann_ar_g(i,j)
                manncc_g= manncc_ar_g(i,j)

                hbf= (tw_g - bo_g)/(2.0*traps_g) !* bankfull depth
                !* main channel
                step= hbf/real(0.5*dimLT)
                h_m=0.0
                do ilt=1, 0.5*dimLT
                    call ufQ_tmrf(h_m, qnorm)
                    normdepthLT_g(i,j,ilt,1)= qnorm
                    normdepthLT_g(i,j,ilt,2)= h_m
                    h_m= h_m+ step
                enddo
                h_m= hbf
                call ufQ_tmrf(h_m, qnorm)
                ilt=0.5*dimLT
                normdepthLT_g(i,j,ilt,1)= qnorm
                normdepthLT_g(i,j,ilt,2)= h_m
                !* flood plain
                step= 0.2*hbf
                h_f= hbf+step
                ict=0
                ict_1=0.25*(0.5*dimLT)
                ict_2=0.5*(0.5*dimLT)
                ict_3=0.75*(0.5*dimLT)
                do ilt=int(0.5*dimLT)+1, dimLT
                    call ufQ_tmrf(h_f, qnorm)
                    normdepthLT_g(i,j,ilt,1)= qnorm
                    normdepthLT_g(i,j,ilt,2)= h_f
                    ict=ict+1
                    if ((ict.gt.ict_1).and.(ict.le.ict_2)) then
                        step=0.5*hbf
                    elseif ((ict.gt.ict_2).and.(ict.le.ict_3)) then
                        step=1.0*hbf
                    elseif (ict.gt.ict_3) then
                        step=5.0*hbf
                    endif
                    h_f= h_f+ step
                enddo
            end do
        end do
        !**----------------------------------------------------------------------------------------+
        !*                      Uniform domain of channel reach network
        !*
        !*  To apply CN method on uniform space domain, the number of new nodes are
        !*  computed by multiplying the number of current nodes of a reach by at least 1.5 (
        !*  that is, mainly to make the minimum number of current nodes from two to three in the
        !*  new domain). Then, a uniform length of any segment between nodes of a given reach is
        !*  equal to total length divided by the number of new nodes.
        !**----------------------------------------------------------------------------------------+
        allocate (ncomp_m_g(nrch_g), dmy2d(100,nrch_g) )
        mxncomp_m_g=0
        do j=1,nrch_g
            ncomp= frnw_g(j,1)
            sumdx= sum(dx_ar_g(1:ncomp-1,j))
            nseg_m=ceiling(1.5*ncomp)+1
            dxst= sumdx/real(nseg_m)
            ncomp_m_g(j)= nseg_m+1
            if (ncomp_m_g(j).gt.mxncomp_m_g) then
                mxncomp_m_g= ncomp_m_g(j)
            end if
            do i=1, ncomp_m_g(j)
                dmy2d(i,j)=dxst
            enddo
        enddo

        allocate (dx_ar_m_g(mxncomp_m_g,nrch_g))
        allocate (bo_ar_m_g(mxncomp_m_g,nrch_g), traps_ar_m_g(mxncomp_m_g,nrch_g))
        allocate (tw_ar_m_g(mxncomp_m_g,nrch_g), twcc_ar_m_g(mxncomp_m_g,nrch_g))
        allocate (mann_ar_m_g(mxncomp_m_g,nrch_g), manncc_ar_m_g(mxncomp_m_g,nrch_g))
        allocate (adjz_ar_m_g(mxncomp_m_g,nrch_g), adjso_ar_m_g(mxncomp_m_g,nrch_g) )
        !* new channel geometries interpolated by the new node configuration (uniform space domain)
        do j=1,nrch_g
            do i=1,ncomp_m_g(j)
                dx_ar_m_g(i,j)= dmy2d(i,j)
            end do
        end do
        do j=1, nrch_g
            ncomp=frnw_g(j,1)
            ncomp_m= ncomp_m_g(j)
            call mod_chgeo(j, ncomp, ncomp_m)
        enddo
        !* new lateral flow interpolated for new node configuration.
        allocate(qlatn_m_g(mxncomp_m_g,nrch_g))
        allocate(dmy1d(nrch_g))
        qlatn_m_g= 0.0  !* initialize
        do j=1,nrch_g
            ncomp=frnw_g(j,1)
            allocate(qlatj_g(ncomp))
            qlatj_g=0.0
            dmy1=0.0
            do i=1,ncomp
                qlatj_g(i)= qlat_g(1,i,j)
                !* test mass conservation
                if (i<ncomp) then
                    dmy1= dmy1+qlatj_g(i)*dx_ar_g(i,j)
                endif
            enddo

            ncomp_m= ncomp_m_g(j)
            allocate(qlatj_m_g(ncomp_m))
            qlatj_m_g=0.0

            call mod_qlatj(j, ncomp, ncomp_m)

            dmy2=0.0
            do i_m=1, ncomp_m
                qlatn_m_g(i_m,j) = qlatj_m_g(i_m)
                !* test mass conservation
                if (i_m<ncomp_m) then
                    dmy2= dmy2 + qlatj_m_g(i_m)*dx_ar_m_g(i_m,j)
                endif
            enddo
            !* test mass conservation (tested correct)
            dmy3=dmy1-dmy2
            dmy1d(j)=dmy3
            deallocate(qlatj_g, qlatj_m_g)
        enddo
        !**-----------------------------------------------------------------------------------*
        !*     INITIAL CONDITIONS of DISCHARGE either on original or uniform space domain
        !*
        !**-----------------------------------------------------------------------------------*
        allocate(q_m_g(mxncomp_m_g,nrch_g))
        ts=1
        do j=1, nrch_g
        !* For the head node of a reach,
            if (frnw_g(j,3)==0) then !* frnw_g(j,3) shows the number of upstream reaches.
            !* head water reach
               i=1
               q_g(i,j)= ubcd_g(ts,j)       !* on original domain
               q_m_g(i,j)= ubcd_g(ts,j)     !* on uniform domain
            else
            !* at a junction
                !* on original domain
                qjt= 0.0
                nusrch= frnw_g(j,3) !* then number of upstream reaches
                do rch=1, nusrch
                    usrchj= frnw_g(j,3+rch) !* js corresponding to upstream reaches
                    ncomp_usrchj= frnw_g(usrchj,1)
                    qjt= qjt + q_g(ncomp_usrchj, usrchj)
                enddo
                i=1
                q_g(i,j)= qjt
                !* on uniform domain
                qjt= 0.0
                nusrch= frnw_g(j,3) !* then number of upstream reaches
                do rch=1, nusrch
                    usrchj= frnw_g(j,3+rch) !* js corresponding to upstream reaches
                    ncomp_m_usrchj= ncomp_m_g(usrchj)
                    qjt= qjt + q_m_g(ncomp_m_usrchj, usrchj)
                enddo
                i=1
                q_m_g(i,j)= qjt
            endif
        !* For the following nodes within a reach after the head node,
            !* on original domain
            if (abs(q_g(1,j))<minq_ns_g) then
                q_g(1,j)=sign(minq_ns_g, q_g(1,j))
            end if
            do i=2, frnw_g(j,1)
                q_g(i,j) = q_g(i-1,j) + qlat_g(ts,i-1,j)*dx_ar_g(i-1,j) !* qlat_g(n,i,j) in unit of m2/sec
                !* keep lower limit of q for numerical stability
                if (abs(q_g(i,j))<minq_ns_g) then
                    q_g(i,j)=sign(minq_ns_g, q_g(i,j))
                end if
            enddo
            !* on uniform domain
            if (abs(q_m_g(1,j))<minq_ns_g) then
                q_m_g(1,j)=sign(minq_ns_g, q_m_g(1,j))
            end if
            do i=2, ncomp_m_g(j)
                q_m_g(i,j) = q_m_g(i-1,j) + qlatn_m_g(i-1,j)*dx_ar_m_g(i-1,j)
                !* keep lower limit of q for numerical stability
                if (abs(q_m_g(i,j))<minq_ns_g) then
                    q_m_g(i,j)=sign(minq_ns_g,q_m_g(i,j))
                end if
            enddo
        enddo
        !**---------------------------------------------------------------------------------*
        !*                  INITIAL CONDITION of water depth

        !*      assuming depths of all reaches and nodes are equal to TW water depth.
        !**---------------------------------------------------------------------------------*
        allocate(elv_m_g(mxncomp_m_g,nrch_g))
        !* 1. use water depth data at TW node
!        do j=1, nrch_g
!            !* on original domain
!            do i=1,frnw_g(j,1)
!                elv_g(i,j)= dbcd_g(1)+ adjz_ar_g(i,j)
!            enddo
!            !* uniform domain
!            do i=1,ncomp_m_g(j)
!                elv_m_g(i,j)= dbcd_g(1)+ adjz_ar_m_g(i,j)
!            enddo
!        enddo
        !* 2. use normal depth at  TW node
        do j=nrch_g, 1, -1
            if (frnw_g(j,2)<0.0) then
                twj= j !* TW reach index
                exit
            endif
        enddo
        ncomp=frnw_g(twj, 1)
        Qxs= q_g(ncomp, twj)
        y_norm= normdepth(ncomp, twj, Qxs)
        do j=1, nrch_g
            !* on original domain
            do i=1,frnw_g(j,1)
                elv_g(i,j)= y_norm+ adjz_ar_g(i,j)
            enddo
            !* uniform domain
            do i=1,ncomp_m_g(j)
                elv_m_g(i,j)= y_norm+ adjz_ar_m_g(i,j)
            enddo
        enddo
        !**---------------------------------------------------------------------------------*
        !*                 Q by draining lateral flow with zero lagging time.
        !*               This computation is done for checking mass conservation.
        !**---------------------------------------------------------------------------------*
        initval= 0.0 !* initial value after allocating variables
        allocate(q_drainql(nts_ql_g, mxncomp_g, nrch_g))
        open(unit=331, file="./temp_test/q by draining lateral w zero lagging.txt")
        write(331,"(A10, 3A10, A20)") "tc", "ts_ev", "i", "j", "q_drainql[m^3/s]"
        tc = t0_g*60.0  !* t0 is in hour. tc is in minutes
        ts=1            !* simulation time step index in parallel with tc.
	ts_ev=1
        do while ( tc .lt. tfin_g*60.0)
            do j=1, nrch_g
                ncomp=frnw_g(j,1)
            !* For the head node of a reach:
                if (frnw_g(j,3)==0) then !* frnw_g(j,3) indicates the number of upstream reaches.
                !* head water reach
                    do n=1,nts_ub_g
                        varr_ub(n)= ubcd_g(n,j) !* ubcd_g in [m^3/sec]
                    enddo
                    q_drainql(ts,1,j)= intp_y(nts_ub_g, tarr_ub, varr_ub, tc)
                else
                !* at a junction
                    qjt= 0.0
                    nusrch= frnw_g(j,3) !* then number of upstream reaches
                    do rch=1, nusrch
                        usrchj= frnw_g(j,3+rch) !* j index of upstream reaches
                        ncomp_usrchj= frnw_g(usrchj,1)
                        qjt= qjt + q_drainql(ts,ncomp_usrchj, usrchj)
                    enddo
                    q_drainql(ts,1,j)= qjt
                endif
            !* For the following nodes within a reach after the head node:
                !* interpolate qlat_g at the current time tc
                allocate(qlatj_g(ncomp))
                qlatj_g=initval !* initialize elements' dummy values
                do i=1,ncomp-1
                    do n=1,nts_ql_g
                        varr_ql(n)= qlat_g(n,i,j) !* qlat_g(n,i,j) in unit of m2/sec
                    enddo
                    qlatj_g(i)= intp_y(nts_ql_g, tarr_ql, varr_ql, tc)
                enddo
                do i=2, ncomp
                    q_drainql(ts,i,j) = q_drainql(ts,i-1,j) + qlatj_g(i-1)*dx_ar_g(i-1,j)
                enddo
                deallocate(qlatj_g)
            enddo
            !* test
	    saveInterval_min=saveInterval_ev_g/60.0 !* [min]
            rmnd= mod(tc,saveInterval_min)  !* computes the remainder of the division of tc by saveinterval_min
            if ( (rmnd<0.0001).or.(tc==tfin_g*60.0) ) then
                do j=1, nrch_g
                    ncomp= frnw_g(j,1)
                    do i=1, ncomp
                         write(331,"(f10.3, 3I10, F20.5)") tc, ts_ev, i, j, q_drainql(ts,i,j)
                    enddo
                enddo
                ts_ev=ts_ev+1
            endif

            tc= tc + dt_ql_g/60.0 !* [min]
            ts= ts+1
        enddo !* do while ( tc .lt. tfin_g*60.0)

        ts_ev=1 !* output time step index for recording q and elv
        do j=1, nrch_g
            do i=1,frnw_g(j,1)
                q_ev_g(ts_ev, i, j)=q_g(i,j)
                elv_ev_g(ts_ev, i, j)= z_ar_g(i,j)+ elv_g(i,j)-adjz_ar_g(i,j)
            enddo
        end do

        tc = t0_g*60.0  !* t0 is in hour. tc is in minutes
        ts=1    !*simulation time step in parallel with tc.
        ts_ev=ts_ev+1   !* time index for recording output
        allocate(c_g(mxncomp_m_g,nrch_g), d_g(mxncomp_m_g,nrch_g)) !* C and D variables of per matrix equation
        allocate(qn_m_g(mxncomp_m_g,nrch_g))
        allocate(elvn_m_g(mxncomp_m_g,nrch_g), elvnm1_m_g(mxncomp_m_g,nrch_g))
        allocate (pqrapj_g(nrch_g))
        allocate ( qnm1_g(mxncomp_g,nrch_g))
        allocate (elvnm1_g(mxncomp_g,nrch_g))
        allocate (qlatn_g(mxncomp_g,nrch_g))
        elvnm1_m_g=elv_m_g !* store current values

        do while (tc .lt. tfin_g*60.0)
            !**---------------------------------------------------------------------------------------------
            !*   discharge is computed by solving CNX with resolution technique on uniform grid, p91, RM5
            !**---------------------------------------------------------------------------------------------
            !* store q and elv on uniform grid at current time n
            do j=1, nrch_g
                do i=1, ncomp_m_g(j)
                    qn_m_g(i,j)= q_m_g(i,j)
                    elvn_m_g(i,j)= elv_m_g(i,j)
                enddo
            end do
            !* upstream/downstream boundary condition for Q at reach j at time ts+1.
            do j=1, nrch_g
                if (frnw_g(j,3)==0) then !* frnw_g(j,3) indicates the number of upstream reaches.
                !* head water reach
                    do n=1,nts_ub_g
                        varr_ub(n)= ubcd_g(n,j)
                    enddo
                    tf0= tc+dtini_g/60.0
                    q_m_g(1,j)= intp_y(nts_ub_g, tarr_ub, varr_ub, tf0) !* tarr_ub in min.
                endif
            enddo
            !* estimate lateral flow at current time n (=tc) for all the nodes of all reaches.
            !* interpolated lateral flow for uniform space domain, p80,RM5
            do j=1,nrch_g
                ncomp=frnw_g(j,1)
                allocate(qlatj_g(ncomp))
                qlatj_g=0.0
                dmy1=0.0
                do i=1,ncomp
                    do n=1,nts_ql_g
                        varr_ql(n)= qlat_g(n,i,j) !* qlat_g(n,i,j) in unit of m2/sec
                    enddo
                    qlatj_g(i)= intp_y(nts_ql_g, tarr_ql, varr_ql, tc) !* interpolated at current time.
                    !* test mass conservation
                    if (i<ncomp) then
                        dmy1= dmy1+qlatj_g(i)*dx_ar_g(i,j)
                    endif
                enddo

                ncomp_m= ncomp_m_g(j)
                allocate(qlatj_m_g(ncomp_m))
                qlatj_m_g=0.0

                call mod_qlatj(j, ncomp, ncomp_m)

                dmy2=0.0
                do i_m=1, ncomp_m
                    qlatn_m_g(i_m,j) = qlatj_m_g(i_m)
                    !* test mass conservation
                    if (i_m<ncomp_m) then
                        dmy2= dmy2 + qlatj_m_g(i_m)*dx_ar_m_g(i_m,j)
                    endif
                enddo
                !* test mass conservation (tested correct)
                dmy3=dmy1-dmy2
                dmy1d(j)=dmy3
                deallocate(qlatj_g, qlatj_m_g)
            enddo

            call cnx_main

            !* keep lower limit of q for numerical stability
            do j=1,nrch_g
                do i=1,ncomp_m_g(j)
                    if (q_m_g(i,j)<minq_ns_g) then
                        q_m_g(i,j)= minq_ns_g
                    endif
                enddo
            enddo
            !**--------------------------------------------------------------------------------------------
            !* Map back q_m_g on uniform space domain to q_g on original domain  for elevation computation
            !**--------------------------------------------------------------------------------------------
            do j=1,nrch_g
                ncomp=frnw_g(j,1)
                ncomp_m=ncomp_m_g(j)
                call mapback_q(j, ncomp, ncomp_m)
            enddo
            !*---------------------------------------------------------------------------------------
            !* Compute y(=elv) or water depth at ts+1 from downstream reach to upstream reach while
            !* computing from downstream node to upstream node for each reach.
            !*---------------------------------------------------------------------------------------
            do j=nrch_g, 1, -1
                ncomp= frnw_g(j,1)
                !* downstream boundary condition for elv at ts+1
                if (frnw_g(j,2)<0.0) then
                !* downstream boundary node at TW
                !* 1. measured water depth data at TW node
!                    do n=1,nts_db_g
!                        varr_db(n)= dbcd_g(n) + adjz_ar_g(ncomp,j) !* when dbcd is water stage, channel bottom elev is added.
!                    enddo
!                    tf0= tc+dtini_g/60.0
!                    elv_g(ncomp,j)= intp_y(nts_db_g, tarr_db, varr_db, tf0)
                !* 2. use normal depth at  TW node
                    ncomp=frnw_g(j,1)
                    Qxs= q_g(ncomp,j)
                    y_norm= normdepth(ncomp, j, Qxs)
                    elv_g(ncomp,j) = y_norm + adjz_ar_g(ncomp,j)
                else
                !* downstream end node of a reach at a junction
                    dsrchj= frnw_g(j,2)    !* j index for downstream reach
                    elv_g(ncomp,j)= elv_g(1,dsrchj)
                endif

                call elv_calc(ncomp, j) !* -> elv_g
            enddo
            !* interpolate water elevation from original domain to uniform domain (elv_g -> elv_m_g)
            do j=1, nrch_g
                ncomp=frnw_g(j,1)
                ncomp_m= ncomp_m_g(j)
                call mod_elv(j, ncomp, ncomp_m)
            enddo
            !* store results of time n-1 w.r.t. the next time computation
            do j=1, nrch_g
                do i=1, ncomp_m_g(j)
                     elvnm1_m_g(i,j)=elvn_m_g(i,j)
                enddo
            end do

            tc= tc + dtini_g/60.0 !* [min]
            ts= ts+1
                !*test
                !do j=1, nrch_g
                !    i=1
                !    print*, tc, i, j, q_g(i,j), elv_g(i,j)-adjz_ar_g(i,j)
                !end do

            saveInterval_min=saveInterval_ev_g/60.0 !* [min]
            rmnd= mod(tc,saveInterval_min)  !* computes the remainder of the division of tc by saveinterval_min
            if ( (rmnd<0.0001).or.(tc==tfin_g*60.0) ) then
                do j=1, nrch_g
                    ncomp= frnw_g(j,1)
                    do i=1, ncomp
                        q_ev_g(ts_ev, i, j)=q_g(i,j)
                        elv_ev_g(ts_ev, i, j)= z_ar_g(i,j)+ elv_g(i,j)-adjz_ar_g(i,j)
                        !write(2,"(F8.1, 2I10, 2F20.4, 2f20.7)") tc,  i, j, q_g(i,j), elv_g(i,j)-adjz_ar_g(i,j), &
                        !                                            c_g(i,j), d_g(i,j)
                    enddo
                enddo
                ts_ev=ts_ev+1
            endif
        enddo !* end of simulation time loop

        deallocate(dx_ar, bo_ar, traps_ar, tw_ar, twcc_ar, mann_ar, manncc_ar)
        deallocate(frnw_g)
        deallocate(q_g, elv_g, q_pre_g, elv_pre_g)
        deallocate(tarr_ql, varr_ql, tarr_ub, varr_ub, tarr_db, varr_db)
        deallocate(adjso_ar_g, adjz_ar_g)
        deallocate(normdepthLT_g)
        deallocate(ncomp_m_g, dmy2d)
        deallocate(dx_ar_m_g, bo_ar_m_g, traps_ar_m_g, tw_ar_m_g, twcc_ar_m_g)
        deallocate(mann_ar_m_g, manncc_ar_m_g, adjz_ar_m_g, adjso_ar_m_g)
        deallocate(qlatn_m_g, dmy1d, q_m_g, elv_m_g)
        deallocate(q_drainql)
        deallocate(c_g, d_g)
        deallocate(qn_m_g, elvn_m_g, elvnm1_m_g)
        deallocate(pqrapj_g)
        deallocate(qnm1_g, elvnm1_g)
        deallocate(qlatn_g)

    endsubroutine diffnw
    !* ----------------------------------------------------------------------------------
    !*      Interpolate bo, traps, tw, twcc, n, ncc for uniform space domain, p80,RM5
    !*
    !* ----------------------------------------------------------------------------------
    subroutine mod_chgeo(j, ncomp, ncomp_m)
        implicit none
        integer,intent(in) :: j, ncomp, ncomp_m
        integer :: i, i_m, i_s, ifc
        double precision :: x_m, x_mflag, sumx1, sumx2, dx_s, u1, u2, u_m
        !* first node
        bo_ar_m_g(1,j)= bo_ar(1,j)
        traps_ar_m_g(1,j)= traps_ar(1,j)
        tw_ar_m_g(1,j)= tw_ar(1,j)
        twcc_ar_m_g(1,j)= twcc_ar(1,j)
        mann_ar_m_g(1,j)= mann_ar(1,j)
        manncc_ar_m_g(1,j)= manncc_ar(1,j)
        adjz_ar_m_g(1,j)= adjz_ar_g(1,j)
        adjso_ar_m_g(1,j)= adjso_ar_g(1,j)
        !* last node
        bo_ar_m_g(ncomp_m,j)= bo_ar(ncomp,j)
        traps_ar_m_g(ncomp_m,j)= traps_ar(ncomp,j)
        tw_ar_m_g(ncomp_m,j)= tw_ar(ncomp,j)
        twcc_ar_m_g(ncomp_m,j)= twcc_ar(ncomp,j)
        mann_ar_m_g(ncomp_m,j)= mann_ar(ncomp,j)
        manncc_ar_m_g(ncomp_m,j)= manncc_ar(ncomp,j)
        adjz_ar_m_g(ncomp_m,j)= adjz_ar_g(ncomp,j)
        adjso_ar_m_g(ncomp_m,j)= adjso_ar_g(ncomp,j)
        !* in-between nodes
        x_m=0.0
        do i_m=2, ncomp_m-1
            x_m= x_m + dx_ar_m_g(i_m-1,j)
            x_mflag=0
            i=1
            sumx1=0.0
            sumx2=0.0
            do while ((x_mflag.eq.0).and.(i.le.ncomp-1))
                sumx1= sumx1 + dx_ar(i,j)
                sumx2= sumx1 + dx_ar(i+1,j)
                if (x_m.le.dx_ar(1,j)) then
                    i_s=1
                    sumx1=0.0
                    sumx2= dx_ar(1,j)
                    x_mflag=1
                elseif ((x_m.gt.sumx1).and.(x_m.le.sumx2)) then
                    i_s= i+1
                    x_mflag=1
                endif
                i=i+1
            enddo
            dx_s= x_m- sumx1
            do ifc=1,8
                if (ifc==1) then
                    u1= bo_ar(i_s,j)
                    u2= bo_ar(i_s+1,j)
                elseif (ifc==2) then
                    u1= traps_ar(i_s,j)
                    u2= traps_ar(i_s+1,j)
                elseif (ifc==3) then
                    u1= tw_ar(i_s,j)
                    u2= tw_ar(i_s+1,j)
                elseif (ifc==4) then
                    u1= twcc_ar(i_s,j)
                    u2= twcc_ar(i_s+1,j)
                elseif (ifc==5) then
                    u1= mann_ar(i_s,j)
                    u2= mann_ar(i_s+1,j)
                elseif (ifc==6) then
                    u1= manncc_ar(i_s,j)
                    u2= manncc_ar(i_s+1,j)
                elseif (ifc==7) then
                    u1= adjz_ar_g(i_s,j)
                    u2= adjz_ar_g(i_s+1,j)
                elseif (ifc==8) then
                    u1= adjso_ar_g(i_s,j)
                    u2= adjso_ar_g(i_s+1,j)
                endif

                u_m= (u2-u1)*dx_s/dx_ar(i_s,j) + u1

                if (ifc==1) then
                    bo_ar_m_g(i_m,j) = u_m
                elseif (ifc==2) then
                    traps_ar_m_g(i_m,j) = u_m
                elseif (ifc==3) then
                    tw_ar_m_g(i_m,j) = u_m
                elseif (ifc==4) then
                    twcc_ar_m_g(i_m,j) = u_m
                elseif (ifc==5) then
                    mann_ar_m_g(i_m,j) = u_m
                elseif (ifc==6) then
                    manncc_ar_m_g(i_m,j) = u_m
                elseif (ifc==7) then
                    adjz_ar_m_g(i_m,j) = u_m
                elseif (ifc==8) then
                    adjso_ar_m_g(i_m,j) = u_m
                endif
            enddo
        enddo !* do i_m=2, ncomp_m-1
    end subroutine mod_chgeo
    !*--------------------------------------------------------------------------
    !* Compute qlatj for the changed length delxj_m on uniform space domain
    !* Refer to p.2-5, RM5.
    !*--------------------------------------------------------------------------
    subroutine mod_qlatj(j, ncomp, ncomp_m)
        implicit none
        integer,intent(in) :: j, ncomp, ncomp_m
        integer :: i, i_m, i_s, i2, x_mflag
        double precision :: x_m, sumx1, sumx2, sumqlatdx, dx_s, sumqlatx_m, sumqlatdx_m, delxj_m

        delxj_m= dx_ar_m_g(1,j)
        x_m=0.0
        do i_m=1, ncomp_m-1
            x_m= x_m + dx_ar_m_g(i_m,j)
            x_mflag=0
            i=1
            sumx1=0.0
            sumx2=0.0
            do while ((x_mflag.eq.0).and.(i.le.ncomp-1))
                sumx1= sumx1 + dx_ar(i,j)
                sumx2= sumx1 + dx_ar(i+1,j)
                if (x_m.le.dx_ar(1,j)) then
                    i_s=1
                    sumx2= dx_ar(1,j)
                    x_mflag=1
                elseif ((x_m.gt.sumx1).and.(x_m.le.sumx2)) then
                    i_s= i+1
                    x_mflag=1
                endif
                i=i+1
            enddo

            sumqlatdx=0.0
            do i2= 1, i_s
                sumqlatdx= sumqlatdx + qlatj_g(i2)*dx_ar(i2,j)
            enddo
            dx_s= sumx2 - x_m
            !* sum of qlat up to distance x_m from node1(=distance zero)
            sumqlatx_m= sumqlatdx - qlatj_g(i_s)*dx_s
            !* sum of qlat_m up to index i_m - 1
            sumqlatdx_m=0.0
            if (i_m>1) then
                do i2=1, i_m-1
                    sumqlatdx_m = sumqlatdx_m + qlatj_m_g(i2)*delxj_m
                enddo
            endif
            qlatj_m_g(i_m)= (sumqlatx_m - sumqlatdx_m)/delxj_m
        enddo !* do i_m=1, ncomp_m-1
    end subroutine mod_qlatj
    !* ----------------------------------------------------------------------------------*
    !*      Map back q on uniform space domain to original domain (q_m_g -> q_g)
    !*
    !* ----------------------------------------------------------------------------------*
    subroutine mapback_q(j, ncomp, ncomp_m)
        implicit none
        integer,intent(in) :: j, ncomp, ncomp_m
        integer :: xflag, i, i_m, i_m_s
        double precision :: x_o, sumx1_m, sumx2_m, dx_s, u1, u2, u_o, delxj_m

        delxj_m= dx_ar_m_g(1,j)
        q_g(1,j)= q_m_g(1,j)
        q_g(ncomp,j)= q_m_g(ncomp_m,j)
        x_o= 0
        do i=2, ncomp-1
            x_o= x_o + dx_ar(i-1,j) !* x distance from upstream end node of original segment lengths.
            xflag=0
            i_m=0
            sumx1_m=0.0
            sumx2_m=0.0
            do while ((xflag.eq.0).and.(i_m.lt.ncomp_m))
                sumx1_m= delxj_m*real(i_m,KIND(delxj_m))
                sumx2_m= delxj_m*real(i_m+1,KIND(delxj_m))
                if ((x_o.gt.sumx1_m).and.(x_o.le.sumx2_m)) then
                    i_m_s= i_m+1
                    xflag=1
                endif
                i_m= i_m + 1
            enddo
            dx_s=  x_o - sumx1_m
            u1= q_m_g(i_m_s,j)
            u2= q_m_g(i_m_s+1,j)
            u_o= (u2-u1)*dx_s/delxj_m + u1
            q_g(i,j) = u_o
        enddo !* i=2, ncomp-1
    end subroutine mapback_q
    !* --------------------------------------------------------------------------*
    !*      Interpolate elevation for uniform space domain
    !*
    !* --------------------------------------------------------------------------*
    subroutine mod_elv(j, ncomp, ncomp_m)
        implicit none
        integer,intent(in) :: j, ncomp, ncomp_m
        integer :: i, i_m, i_s
        double precision :: x_m, x_mflag, sumx1, sumx2, dx_s, u1, u2, u_m

        elv_m_g(1,j)= elv_g(1,j)            !* first node
        elv_m_g(ncomp_m,j)= elv_g(ncomp,j)  !* last node
        !* in-between nodes
        x_m=0.0
        do i_m=2, ncomp_m-1
            x_m= x_m + dx_ar_m_g(i_m-1,j)
            x_mflag=0
            i=1
            sumx1=0.0
            sumx2=0.0
            do while ((x_mflag.eq.0).and.(i.le.ncomp-1))
                sumx1= sumx1 + dx_ar(i,j)
                sumx2= sumx1 + dx_ar(i+1,j)
                if (x_m.le.dx_ar(1,j)) then
                    i_s=1
                    sumx1=0.0
                    sumx2= dx_ar(1,j)
                    x_mflag=1
                elseif ((x_m.gt.sumx1).and.(x_m.le.sumx2)) then
                    i_s= i+1
                    x_mflag=1
                endif
                i=i+1
            enddo
            dx_s= x_m- sumx1
            u1= elv_g(i_s,j)
            u2= elv_g(i_s+1,j)
            u_m= (u2-u1)*dx_s/dx_ar(i_s,j) + u1
            elv_m_g(i_m,j) = u_m
        enddo
    end subroutine mod_elv
    !**------------------------------------------------------------------------------
    !*         CNX and resolution technique on uniform grid, p106~109,RM5

    !**------------------------------------------------------------------------------
    subroutine cnx_main
        implicit none
        allocate (p(mxncomp_m_g, nrch), qq(mxncomp_m_g, nrch), r(mxncomp_m_g, nrch))
        allocate (w(mxncomp_m_g, nrch))
        !* initialize
        p=0.0
        qq=0.0
        r=0.0
        w=0.0
        !* Compute C and D at (i,j) and n, p27_10~11
        call C_D_calc
        !* Compute p, q, r, and w coefficients for entire network at time n
        call p_q_r_w_calc
        !* Compute Q using resolution technique
        call Q_calc_resoltech
        !call Q_calc_Thomas
        deallocate(p, qq, r, w)
     endsubroutine cnx_main
    !**------------------------------------------------------------------------------
    !*                          Celerity and Diffusivity

    !*                              p88 ~ 89, RM5
    !**------------------------------------------------------------------------------
    subroutine C_D_calc
        implicit none
        integer :: i, j, ncomp_m, time_flag
        double precision :: bi_n, bip1_n
        double precision :: dxx1, dmy
        double precision :: dBdxi
        double precision :: refQ, sumc, avec, sumd, aved
        double precision :: depth_n, depth_nm1, ki_n, ki_nm1, dKdyi

        do j=1, nrch
            ncomp_m= ncomp_m_g(j)
            do i=1, ncomp_m
                !* dB/dx at (i,j) and time n
                time_flag=0 !* for time n
                call twwi(i, j, time_flag, bi_n)
                if (i<ncomp_m) then
                    call twwi(i+1, j, time_flag, bip1_n)
                    dxx1= dx_ar_m_g(i,j)
                    dBdxi=(bip1_n-bi_n)/dxx1
                else
                    dBdxi=0.0
                endif
                !* dK/dy at (i,j)
                depth_n= elvn_m_g(i,j) - adjz_ar_m_g(i,j)
                depth_nm1= elvnm1_m_g(i,j) - adjz_ar_m_g(i,j)
                if (depth_n==depth_nm1) then
                    depth_n=depth_n+ mindepth_ns_g
                endif
                ki_n= convey(i, j, depth_n)
                ki_nm1= convey(i, j, depth_nm1)
                dKdyi= (ki_n - ki_nm1)/(depth_n - depth_nm1)
                refQ= qn_m_g(i,j) !* reference discharge
                !* Diffusivity at (i,j) and time n
                d_g(i,j)= ki_n**2.0/(2.0*bi_n*abs(refQ))
                !* Celerity at (i,j) and time n
                dmy= ki_n**2.0/(2.0*bi_n*refQ)
                c_g(i,j)= dKdyi*refQ/(ki_n*bi_n) + dmy*dBdxi/bi_n
            enddo
            !* use diffusivity and celerity values averaged along all the nodes of a reach
            sumc=  sum(c_g(1:ncomp_m,j))
            avec= sumc/real(ncomp_m)
            sumd= sum(d_g(1:ncomp_m,j))
            aved= sumd/real(ncomp_m)
            do i=1, ncomp_m
                c_g(i,j)=avec
                d_g(i,j)=aved
            enddo
        enddo
    endsubroutine C_D_calc
    !**-------------------------------------------------------------------*
    !*                   p, q, r, w coefficients of CNX

    !*                          p90, RM5
    !**-------------------------------------------------------------------*
     subroutine p_q_r_w_calc
        implicit none
        integer :: i, j, ncomp_m, im1, ip1
        double precision :: dxst, sumc, avec, sumd, aved, hhi, ggi, ppr, qpr, rpr

        do j=1, nrch
            ncomp_m= ncomp_m_g(j)
            dxst=dx_ar_m_g(1,j)    !*uniform segment length for each j reach
            if (ncomp_m.ge.3) then
                do i=2, ncomp_m-1
                    im1= i-1
                    ip1= i+1
                    sumc=  sum(c_g(im1:ip1,j))
                    avec= sumc/real(ip1-im1+1)
                    sumd=  sum(d_g(im1:ip1,j))
                    aved= sumd/real(ip1-im1+1)

                    hhi= avec*dtini/dxst
                    ggi= aved*dtini/(dxst**2.0)

                    p(i,j)= -hhi/4.0 - ggi/2.0
                    qq(i,j)= 1.0 + ggi
                    r(i,j)= hhi/4.0 - ggi/2.0

                    ppr= hhi/4.0 + ggi/2.0
                    qpr= 1.0 - ggi
                    rpr= -hhi/4.0 + ggi/2.0

                    w(i,j)= ppr*qn_m_g(i-1,j) + qpr*qn_m_g(i,j) + rpr*qn_m_g(i+1,j)
                enddo
            endif
        enddo
    endsubroutine p_q_r_w_calc
    !**----------------------------------------------------------------------------------------------*
    !*              Q computation using resolution technique of CNX

    !* source: 1. Moussa at el. "Algorithms for solving the diffusive wave flood routing equation",
    !*              Hydrological Processes. vol 10, 105-123 (1996)
    !**----------------------------------------------------------------------------------------------*
     subroutine Q_calc_resoltech
        implicit none
        integer :: i, j, ncomp_m, nusrch, rch, usrchj, ncomp_usrchj
        integer :: i1, i2, neg_method
        double precision :: dxst
        double precision :: qjt, usposQ, dsposQ
        double precision, dimension(:), allocatable :: x, y

        do j=1, nrch
            ncomp_m= ncomp_m_g(j)
            dxst=dx_ar_m_g(1,j)  !*uniform segment length for each j reach
            !* compute Q at top node and time n+1 (** Q for headwater reach is already computed in the calling subroutine)
            if (frnw_g(j,3)>0) then
                !* At junction.
                qjt= 0.0
                nusrch= frnw_g(j,3)  !* the number of upstream reaches
                do rch=1, nusrch
                    usrchj= frnw_g(j,3+rch)  !* j index corresponding to upstream reaches
                    ncomp_usrchj= ncomp_m_g(usrchj)
                    qjt= qjt + q_m_g(ncomp_usrchj, usrchj)
                enddo
                !* top node of the junction.
                q_m_g(1,j)= qjt
            endif

            if (ncomp_m.ge.3) then
                !* x and y computation in backward-sweep
                allocate(x(ncomp_m), y(ncomp_m))
                x(ncomp_m)= 1.0
                y(ncomp_m)= 0.0
                do i=ncomp_m-1, 2, -1
                    x(i)= -p(i,j)/(qq(i,j) + r(i,j)*x(i+1))
                    y(i)= (w(i,j)- r(i,j)*y(i+1))/(qq(i,j)+r(i,j)*x(i+1))
                end do
                !* Q computation at time n+1 in forward-sweep
                do i=2, ncomp_m
                    q_m_g(i,j)= x(i)*q_m_g(i-1,j) + y(i)
                end do

                deallocate(x, y)
            else
                i=2
                q_m_g(i,j)= q_m_g(i-1,j)
            end if
            !* overlay qlateral flow at time n on top of computed Q at n+1
            do i=2, ncomp_m
                q_m_g(i,j)= q_m_g(i,j) + qlatn_m_g(i-1,j)*dxst
            enddo
            !* negative discharge treatment
            neg_method=1
            if (neg_method==1) then
            !* method 1: interpolate between neighboring positive discharge values
                do i=2, ncomp_m
                    if (q_m_g(i,j)<0.0) then
                        !* find upstream positive discharge
                        i1=i
                        do while(q_m_g(i1,j)<minq_ns_g)
                            i1=i1-1
                            if (i1<1) then
                                exit
                            end if
                        end do
                        if (i1.ge.1) then
                            usposQ= q_m_g(i1,j)
                        else
                            usposQ= minq_ns_g
                            i1=1
                        endif
                         !* find downstream positive discharge
                        i2=i
                        do while(q_m_g(i2,j)<minq_ns_g)
                            i2=i2+1
                            if (i2>ncomp_m) then
                                exit
                            end if
                        end do
                        if (i2.le.ncomp_m) then
                            dsposQ= q_m_g(i2,j)
                        else
                            dsposQ= minq_ns_g
                            i2=ncomp_m
                        endif
                        q_m_g(i,j)= (dsposQ-usposQ)*real(i-i1)/real(i2-i1) + usposQ
                    endif
                enddo
            elseif (neg_method==2) then
            !* method 2: replace by minimum discharge value
                do i=2, ncomp_m
                    if (q_m_g(i,j).le.0.0) then
                        q_m_g(i,j)=minq_ns_g
                    endif
                enddo
            endif
        enddo
     endsubroutine Q_calc_resoltech
    !**---------------------------------------------------------------
    !*      compute area only for a given i and j
    !*
    !**---------------------------------------------------------------
    subroutine areai(i, j, time_flag, ai)
        implicit none
        integer, intent(in) :: i, j, time_flag
        double precision, intent(out) :: ai
        double precision :: depth

        bo_g= bo_ar_m_g(i,j)
        traps_g= traps_ar_m_g(i,j)
        tw_g= tw_ar_m_g(i,j)
        twcc_g= twcc_ar_m_g(i,j)
        if (time_flag==-1) then
        !* time n-1
            depth= elvnm1_m_g(i, j)-adjz_ar_m_g(i,j)
            call areacalc(depth, ai)
        elseif (time_flag==0) then
        !* time n
            depth= elvn_m_g(i, j)-adjz_ar_m_g(i,j)
            call areacalc(depth, ai)
        endif
    endsubroutine areai
    !**---------------------------------------------------------------
    !*      compute water top width for a given i and j
    !*
    !**---------------------------------------------------------------
    subroutine twwi(i, j, time_flag, bi)
        implicit none
        integer, intent(in) :: i, j, time_flag
        double precision, intent(out) :: bi
        double precision :: depth

        bo_g= bo_ar_m_g(i,j)
        traps_g= traps_ar_m_g(i,j)
        tw_g= tw_ar_m_g(i,j)
        twcc_g= twcc_ar_m_g(i,j)

        if (time_flag==-1) then
        !* time n-1
            depth= elvnm1_m_g(i, j)-adjz_ar_m_g(i,j)
            bi= topwidth(depth)
        elseif (time_flag==0) then
        !* time n
            depth= elvn_m_g(i, j)-adjz_ar_m_g(i,j)
            bi= topwidth(depth)
        endif
    endsubroutine twwi
    !**---------------------------------------------------------------
    !*      compute conveyance for a given depth at i,j,n
    !*
    !**---------------------------------------------------------------
    doubleprecision function convey(i, j, depthi)
        implicit none
        integer, intent(in) :: i, j
        double precision, intent(in) :: depthi
        double precision :: areai, hydR, conv
        !* conveyance at (i,j,n)
        bo_g= bo_ar_m_g(i,j)
        traps_g= traps_ar_m_g(i,j)
        tw_g= tw_ar_m_g(i,j)
        twcc_g= twcc_ar_m_g(i,j)
        mann_g= mann_ar_m_g(i,j)
        manncc_g= manncc_ar_m_g(i,j)
        call areacalc(depthi, areai)
        call hydRcalc(depthi, areai, hydR)
        call Kcalc(depthi, areai, hydR, conv)
        convey= conv
    endfunction convey
    !**--------------------------------------------------------------
    !*          Computation of top width of water surface
    !**--------------------------------------------------------------
    doubleprecision function topwidth(hxs)
        implicit none
        doubleprecision, intent(in) :: hxs
        doubleprecision :: bwd, sslp, hbf, tw, twcc
        bwd= bo_g !*bottom width
        sslp= traps_g !*trapezoidal
        tw= tw_g
        twcc=twcc_g
        if (sslp==0.0) then
        !* rectangular channel
            topwidth= bwd
        else
            hbf= (tw - bwd)/(2.0*sslp) !* bankfull depth
            if (hxs.le.hbf) then
                !* trapezoidal channel inbank flow
                topwidth= bwd + 2.0*sslp*hxs
            else
                !*overbank flow on rect. floodplains
                topwidth= twcc
            end if
        endif
    endfunction topwidth
    !*----------------------------------------------------------------------------------------
    !*      Compute water elevation from downstream to upstream based on normal depth

    !*-----------------------------------------------------------------------------------------
    subroutine elv_calc(ncomp, j)
        implicit none
        integer,intent(in) :: ncomp, j
        integer :: i
        doubleprecision :: qi, yi_tp1

        do i=ncomp-1,1,-1
            qi= q_g(i,j)
            yi_tp1= normdepth(i, j, qi)
            !* make sure water depth larger than arbitrarily chosen minimum positive water depth.
            if (yi_tp1.lt.mindepth_ns_g) then
                yi_tp1= mindepth_ns_g
            end if
            z_g= adjz_ar_g(i,j)
            elv_g(i,j)= yi_tp1 + z_g
        end do
    end subroutine elv_calc
    !**----------------------------------------------------------------
    !*      compute normal depth using Bisection method

    !*-----------------------------------------------------------------
    doubleprecision function normdepth(i, j, Qxs)
        implicit none
        integer, intent(in) :: i, j
        double precision, intent(in) :: Qxs
        double precision, dimension(:), allocatable :: xarr, yarr
        integer :: ilt, dimLT

        dimLT= dim_nmLT_g
        allocate(xarr(dimLT), yarr(dimLT))
        do ilt=1, dimLT
            xarr(ilt)= normdepthLT_g(i,j,ilt,1) !* uniform discharge
            yarr(ilt)= normdepthLT_g(i,j,ilt,2) !* normal depth
        enddo
        normdepth= intp_y(dimLT, xarr, yarr, Qxs)
        deallocate(xarr,yarr)
    endfunction normdepth
    !*--------------------------------------------
    !           Interpolate y for given x
    !
    !*--------------------------------------------
    doubleprecision function intp_y(nrow, xarr, yarr, x)
        implicit none
        integer, intent(in) :: nrow
        doubleprecision, dimension(nrow), intent(in) :: xarr, yarr
        doubleprecision, intent(in) :: x
        integer :: irow
        doubleprecision :: x1, y1, x2, y2, y

        irow= locate(xarr, x)
        if (irow.eq.0) irow= 1
        if (irow.eq.nrow) irow= nrow-1
        x1= xarr(irow); y1= yarr(irow)
        x2= xarr(irow+1); y2= yarr(irow+1)
        y= LInterpol(x1,y1,x2,y2,x)
        intp_y = y
    end function intp_y
    !*------------------------------------------------------------------------------------
    !                       Locate function in f90, p.1045,NR f90
    !
    !   klo=max(min(locate(xa,x),n-1),1) In the Fortran 77 version of splint,
    !   there is in-line code to find the location in the table by bisection. Here
    !   we prefer an explicit call to locate, which performs the bisection. On
    !   some massively multiprocessor (MMP) machines, one might substitute a different,
    !   more parallel algorithm (see next note).
    !*------------------------------------------------------------------------------------
    integer function locate(xx,x)
        implicit none
        doubleprecision, dimension(:), intent(in) :: xx
        doubleprecision, intent(in) :: x
        !* Given an array xx(1:N), and given a value x, returns a value j such that x is between
        !* xx(j) and xx(j + 1). xx must be monotonic, either increasing or decreasing.
        !* j = 0 or j = N is returned to indicate that x is out of range.
        integer :: n,jl,jm,ju
        logical :: ascnd

        n=size(xx)
        ascnd = (xx(n) >= xx(1))  !* True if ascending order of table, false otherwise.
        jl=0    !* Initialize lower
        ju=n+1  !* and upper limits.
        do
            if (ju-jl <= 1) exit    !* Repeat until this condition is satisfied.
            jm=(ju+jl)/2            !* Compute a midpoint,
            if (ascnd .eqv. (x >= xx(jm))) then
                jl=jm               !* and replace either the lower limit
            else
                ju=jm               !* or the upper limit, as appropriate.
            end if
        end do

        if (x == xx(1)) then        !* Then set the output, being careful with the endpoints.
            locate=1
        else if (x == xx(n)) then
            locate=n-1
        else
            locate=jl
        end if
    end function locate
    !*--------------------------------------------------
    !*                 Linear Interpolation
    !
    !*--------------------------------------------------
    doubleprecision function LInterpol(x1,y1,x2,y2,x)
        implicit none
        doubleprecision, intent(in) :: x1, y1, x2, y2, x
        !* interpolate y for the given x
        LInterpol= (y2-y1)/(x2-x1)*(x-x1)+y1
    end function LInterpol
    !*----------------------------------------------------------------------
    !*     Uniform flow discharge for trapz. main and rect. floodplain
    !*     p.102, RM4
    !*----------------------------------------------------------------------
    subroutine ufQ_tmrf(hxs, ufQ)
        implicit none
        double precision, intent(in) :: hxs
        double precision, intent(out) :: ufQ
        double precision :: Axs, hydR, WP, Twd, TCCwd, So
        double precision :: bwd, sslp, hbf, ufQ1, ufQ2, ufQ3
        bwd= bo_g
        sslp= traps_g
        Twd= tw_g       !* top width of main channel
        TCCwd= twcc_g   !* top width of compound channel
        So= so_g
        hbf= (Twd - bwd)/(2.0*sslp) !* bankfull depth
        if (hxs.le.hbf) then
        !* trapezoidal channel inbank flow
            Axs=(bwd + sslp*hxs)*hxs
            WP= bwd + 2.0*hxs*((1.0+(sslp**2.0))**0.5)
            hydR=Axs/WP
            ufQ1=0.0
            ufQ3=0.0
            ufQ2= (1.0/mann_g)*Axs*(hydR**(2.0/3.0))*(So**0.5)
        else
        !*overbank flow on rect. floodplains
            !* subsection 1 and 3, p102,RM4
            Axs= (hxs-hbf)*(TCCwd -Twd)/2.0
            WP= (TCCwd - Twd)/2.0 + (hxs-hbf)
            hydR= Axs/WP
            ufQ1= (1.0/manncc_g)*Axs*(hydR**(2.0/3.0))*(So**0.5)
            ufQ3= ufQ1
            !* subsection 2, p102,RM4
            Axs= (bwd + sslp*hbf)*hbf + (hxs-hbf)*Twd
            WP= bwd + 2.0*hbf*((1.0 + sslp**2.0)**0.5)
            hydR= Axs/WP
            ufQ2= (1.0/mann_g)*Axs*(hydR**(2.0/3.0))*(So**0.5)
        end if

        ufQ= ufQ1+ufQ2+ufQ3
    end subroutine ufQ_tmrf
    !**------------------------------------------------------------------------------
    !*  Computation of area of various channel x-sections with depth as an argument
    !*
    !**------------------------------------------------------------------------------
    subroutine areacalc(hxs, Axs)
        implicit none
        double precision, intent(in) :: hxs
        double precision, intent(out) :: Axs
        double precision :: bwd, sslp, hbf, tw, twcc
        bwd= bo_g !*bottom width
        sslp= traps_g !*trapezoidal
        tw= tw_g
        twcc=twcc_g
        if (sslp==0.0) then
            !* rectangular channel
            Axs=bwd*hxs
        else
            hbf= (tw - bwd)/(2.0*sslp) !* bankfull depth
            if (hxs.le.hbf) then
                !* trapezoidal channel inbank flow
                Axs=(bwd + sslp*hxs)*hxs
            else
                !*overbank flow on rect. floodplains
                Axs=(bwd + sslp*hbf)*hbf + twcc*(hxs-hbf)
            end if
        endif
    end subroutine areacalc
    !**-----------------------------------------------------
    !*      Computation of hydraulic radius R (=A/P)
    !**-----------------------------------------------------
    subroutine hydRcalc(hxs, Axs, hydR)
        implicit none
        double precision, intent(in) ::  hxs, Axs
        double precision, intent(out) :: hydR
        double precision :: bxs, ssxs, twxs, twccxs
        double precision :: hbf
        bxs=bo_g
        ssxs=traps_g
        twxs= tw_g
        twccxs= twcc_g
        if (ssxs==0.0) then
            !* rectangular channel
            hydR= Axs/(bxs + 2.0*hxs)
        else
            hbf= (twxs - bxs)/(2.0*ssxs)
            if (hxs<=hbf) then
            !* inbank flow in trapezoidal channel
                hydR=Axs/(bxs + 2.0*hxs*((1.0 + ssxs**2.0)**0.5))
            else
            !* compound channel having trapezoidal main and rectangular floodplains, p80-,RM1
                hydR= Axs/(bxs + 2.0*hbf*(1.0 + ssxs**2.0)**0.5 + twccxs - twxs + 2.0*(hxs-hbf))
            endif
        endif
    end subroutine hydRcalc
    !**-----------------------------------------------------
    !* Computation of conveyance K (=1/N * A * R^(2/3))
    !**-----------------------------------------------------
    subroutine Kcalc(hxs, Axs, hydR, cnvey)
        implicit none
        double precision, intent(in) :: hxs, Axs, hydR
        double precision, intent(out) :: cnvey
        double precision :: bwd, sslp, Twd, TCCwd
        double precision :: subA, subP, hbf, TwCCi, K0, K1, K2
        bwd= bo_g
        sslp= traps_g
        Twd=tw_g !* top width of main channel
        TCCwd=twcc_g !* top width of compound channel
        hbf= (Twd - bwd)/(2.0*sslp)!* bankfull hxs
        if ((sslp==0.0).or.(hxs<=hbf)) then
        !* inbank flow in rectangular or trapezoidal channel
            cnvey=(1.0/mann_g)*Axs*(hydR**(2.0/3.0))
        else
        !* overbank flow in compound channel having trapezoidal main channel with rectangular floodplains, p84-2~p84-2-1, RM1
            !* conveyance in the main channel, K0
            subA = (bwd + sslp*hbf)*hbf + (hxs - hbf)*Twd
            subP = bwd + 2.0*hbf*(1.0+ sslp**2.0)**0.5
            K0 = (1.0/mann_g)*subA**(5.0/3.0)/subP**(2.0/3.0)
            !* conveyance in the left floodplain (assuming symmetric), K1
            TwCCi=(TCCwd - Twd)/2.0
            subA=(hxs - hbf)*TwCCi
            subP=TwCCi +  (hxs - hbf)
            K1 = (1.0/manncc_g)*subA**(5.0/3.0)/subP**(2.0/3.0)
            !* conveyance in the left floodplain (assuming symmetric), K1
            TwCCi=(TCCwd - Twd)/2.0
            subA=(hxs - hbf)*TwCCi
            subP=TwCCi +  (hxs - hbf)
            K2 = (1.0/manncc_g)*subA**(5.0/3.0)/subP**(2.0/3.0)

            cnvey=K0+K1+K2
        endif
    endsubroutine Kcalc
end module diffusive
