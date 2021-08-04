module arrays_module

    implicit none
    save

    double precision, allocatable :: area(:), bo(:,:) !, y(:, :), q(:, :)
    !double precision, allocatable :: areap(:), qp(:), z(:), dqp(:)
    double precision, allocatable :: av11(:), av12(:), av21(:), av22(:)
    double precision, allocatable ::  ci1(:), ci2(:)!, dqc(:), dap(:), dac(:)
    double precision, allocatable :: aso(:,:), f1(:), f2(:), depth(:)
    double precision, allocatable :: g11inv(:), g12inv(:), g21inv(:), g22inv(:)
    double precision, allocatable :: b11(:), b12(:), b21(:), b22(:)
    double precision, allocatable :: eps2(:), eps4(:), d1(:), d2(:), u(:), c(:)
    double precision, allocatable :: co(:), gso(:,:), dbdx(:,:)!,sk(:)
    double precision, allocatable :: dx(:,:), froud(:), courant(:)
    double precision, allocatable :: dt(:)

	!**arrays for branching channel application
    integer, allocatable :: ndep(:), uslinks(:,:), dslink(:), instrdflag(:,:), nx1(:)
    !double precision, allocatable :: y(:, :, :), q(:, :, :), qlat(:,:,:), bo(:, :), traps(:,:) !,area(:,:), areafnal(:,:,:),
    double precision, allocatable :: areap(:, :), qp(:, :, :), z(:, :), sk(:, :)            ! change 20210628


    double precision, allocatable :: dqp(:,:), dap(:,:), dqc(:,:), dac(:,:)


    double precision, allocatable :: celerity(:,:),velocity(:,:), diffusivity(:,:), diffusivity2(:), celerity2(:)

    double precision, allocatable :: eei(:), ffi(:), exi(:), fxi(:), qpx(:,:), qcx(:)

    double precision, allocatable :: USBoundary(:,:,:), DSBoundary(:,:,:)
    integer, allocatable :: upBoundTableEntry(:), downBoundTableEntry(:) !, normalDepth(:,:)
! change for unsteady flow
    double precision, allocatable :: pere(:,:),dpda(:)

    double precision, allocatable :: oldQ(:,:), newQ(:,:,:), oldArea(:,:), newArea(:,:)
	double precision, allocatable :: oldY(:,:), newY(:,:), normalDepthAtNodes(:,:)      ! change 20210628
    double precision, allocatable :: added_Q(:,:,:)                                                                                         ! change 20210713

    integer, allocatable :: ityp(:), latFlowLocations(:,:), dataInEachLatFlow(:,:), latFlowType(:,:), latFlowXsecs(:,:)
    double precision, allocatable :: lateralFlowTable(:,:,:,:), lateralFlow(:,:,:)

    ! for additional lateral flow of the structures
    integer, allocatable :: latFlowLocations2(:,:), dataInEachLatFlow2(:,:), latFlowType2(:,:), latFlowXsecs2(:,:), noLatFlow2(:)
    double precision, allocatable :: lateralFlowTable2(:,:,:,:), lateralFlow2(:,:)

    double precision, allocatable :: dimensionless_Cr(:,:), dimensionless_Fo(:,:), dimensionless_Fi(:,:)
    double precision, allocatable :: dimensionless_Di(:,:), dimensionless_Fc(:,:), dimensionless_D(:,:)

    double precision, allocatable :: ini_y(:), ini_q(:)
    double precision, allocatable :: ini_q_repeat(:,:), ini_E(:,:), ini_F(:,:)

    integer, allocatable :: Q_sk_tableEntry(:,:), noLatFlow(:), noQSKtable(:)
    double precision, allocatable :: eachQSKtableNodeRange(:,:,:), Q_sk_Table(:,:,:,:)

    double precision, allocatable :: lowerLimitCount(:), higherLimitCount(:), volRemain(:,:)

    character(len=128), allocatable :: downstream_path(:), xSection_path(:), manning_strickler_path(:), upstream_path(:),dx_path(:)
    character(len=128), allocatable :: QSKtablePath(:), lateralFlow_path(:), lateralFlow_path2(:)
    character(len=128), allocatable :: bankLocation_path(:)
    double precision, allocatable :: leftBank(:,:), rightBank(:,:), skLeft(:,:), skMain(:,:), skRight(:,:)

    integer, allocatable :: currentROutingDiffusive(:), notSwitchRouting(:)
    double precision :: minDx, maxCelerity,maxCelDx

    integer, allocatable :: currentRoutingNormal(:,:), routingNotChanged(:,:)

contains

    ! Allocate storage for all of the arrays in this module based on the number
    ! of time steps and spatial points
!    subroutine setup_arrays(num_time, num_points, maxTableEntry1, maxTableEntry2, totalLatFlow, totalQSKtable, totalChannels)
!        implicit none
!        ! Input
!        integer, intent(in) :: num_time, num_points, maxTableEntry1, maxTableEntry2, totalLatFlow, totalQSKtable, totalChannels
!        allocate(area(num_points))
!! change for unsteady flow
!        allocate(bo(num_points,totalChannels))
!        allocate(pere(num_points,totalChannels))
!        allocate(dpda(num_points))
!        !allocate(normalDepth(totalChannels))    !! this parameter indicates which channel will have full diffusive or partial diffusive routing.
!        allocate(normalDepthAtNodes(num_points,totalChannels))
!        normalDepthAtNodes = 0.
!        allocate(areap(num_points,totalChannels))
!        allocate(qp(num_points,num_time,totalChannels))
!        allocate(ini_q_repeat(num_points,totalChannels))
!        allocate(ini_E(num_points,totalChannels))
!        allocate(ini_F(num_points,totalChannels))
!        allocate(z(num_points,totalChannels))
!        z=0. ! initialization
!        allocate(av11(num_points))
!        allocate(av12(num_points))
!        allocate(av21(num_points))
!        allocate(av22(num_points))
!        allocate(dqp(num_points,totalChannels))
!        allocate(dqc(num_points,totalChannels))
!        allocate(dap(num_points,totalChannels))
!        allocate(dac(num_points,totalChannels))
!        allocate(ci1(num_points))
!        allocate(ci2(num_points))
!        allocate(aso(num_points,totalChannels))
!        aso = 0.
!        allocate(depth(num_points))
!        allocate(f1(num_points))
!        allocate(f2(num_points))
!        allocate(g11inv(num_points))
!        allocate(g12inv(num_points))
!        allocate(g21inv(num_points))
!        allocate(g22inv(num_points))
!        allocate(b11(num_points))
!        allocate(b12(num_points))
!        allocate(b21(num_points))
!        allocate(b22(num_points))
!        allocate(eps2(num_points))
!        allocate(eps4(num_points))
!        allocate(d1(num_points))
!        allocate(d2(num_points))
!        allocate(u(num_points))
!        allocate(c(num_points))
!        allocate(sk(num_points,totalChannels))
!        allocate(leftBank(num_points,totalChannels))
!        allocate(rightBank(num_points,totalChannels))
!        allocate(skLeft(num_points,totalChannels))
!        allocate(skMain(num_points,totalChannels))
!        allocate(skRight(num_points,totalChannels))
!        allocate(co(num_points))
!        allocate(gso(num_points,totalChannels))
!        gso = 0.
!        allocate(dbdx(num_points,totalChannels))
!        dbdx = 0.
!        allocate(dt(num_points))
!        allocate(ityp(num_points))
!        allocate(dx(num_points-1,totalChannels))
!        allocate(volRemain(num_points-1,totalChannels))
!        allocate(froud(num_points))
!        allocate(Q_sk_Table(2, maxTableEntry1, totalQSKtable,totalChannels))
!        allocate(Q_sk_tableEntry(totalQSKtable,totalChannels))
!        allocate(USBoundary(2, maxTableEntry2,totalChannels))
!        allocate(DSBoundary(2, maxTableEntry2,totalChannels))
!        allocate(ndep(totalChannels))
!        allocate(dslink(totalChannels))
!        !allocate(uslinks(num_points,totalChannels))
!        allocate(instrdflag(totalChannels,2))
!        allocate(courant(num_points-1))
!        allocate(oldQ(num_points, totalChannels))
!        allocate(newQ(num_points, num_time, totalChannels))     ! change 20210628
!        allocate(added_Q(num_points, num_time, totalChannels))  ! change 20210713
!        allocate(oldArea(num_points, totalChannels))
!        allocate(newArea(num_points, totalChannels))
!        allocate(oldY(num_points, totalChannels))
!        allocate(newY(num_points, totalChannels))
!        oldQ = -999; oldY = -999; newQ = -999; newY = -999
!        allocate(lateralFlowTable(2, maxTableEntry2, totalLatFlow, totalChannels))
!        allocate(dataInEachLatFlow(totalLatFlow, totalChannels))
!        !allocate(lateralFlow(num_points, totalChannels)) ! change 20210311
!        allocate(lateralFlow(num_points, num_time, totalChannels)) ! change 20210707
!        allocate(celerity(num_points, totalChannels))
!        allocate(velocity(num_points, totalChannels))
!        allocate(diffusivity(num_points, totalChannels))
!        allocate(celerity2(num_points))
!        allocate(diffusivity2(num_points))
!
!        allocate(eei(num_points))
!        allocate(ffi(num_points))
!        allocate(exi(num_points))
!        allocate(fxi(num_points))
!        allocate(qpx(num_points, totalChannels))
!        allocate(qcx(num_points))
!
!        allocate(dimensionless_Cr(num_points-1,totalChannels))
!        allocate(dimensionless_Fo(num_points-1,totalChannels))
!        allocate(dimensionless_Fi(num_points-1,totalChannels))
!        allocate(dimensionless_Di(num_points-1,totalChannels))
!        allocate(dimensionless_Fc(num_points-1,totalChannels))
!        allocate(dimensionless_D(num_points-1,totalChannels))
!        dimensionless_Cr = -999; dimensionless_Fo = -999; dimensionless_Fi = -999
!        dimensionless_Di = -999; dimensionless_Fc = -999; dimensionless_D = -999
!
!        allocate(lowerLimitCount(totalChannels))
!        allocate(higherLimitCount(totalChannels))
!        allocate(currentRoutingNormal(num_points-1,totalChannels))
!        allocate(routingNotChanged(num_points-1,totalChannels))
!
!    end subroutine setup_arrays

end module arrays_module
