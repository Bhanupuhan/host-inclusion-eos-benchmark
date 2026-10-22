! =====================================================================================================
! UMAT  umat_eos_total.f90  -  TOTAL-FORM equation of state for olivine-in-diamond (path independent)
! Author: Biswabhanu Puhan.  Part of host-inclusion-eos-benchmark (see README.md, CHANGELOG section).
!
! Replaces an earlier incremental (hypoelastic) UMAT (not included).  Stress is computed from the TOTAL
! strain and the CURRENT temperature, never by adding increments:
!     sigma = -P_EoS(V, T) I + 2 G dev(eps),     V = V_ref * exp(tr eps),   G constant
! where eps is the total strain measured from the reference state (P_ref, T_ref) at which the analysis
! starts (the entrapment state, set in the input with *Initial Conditions, type=STRESS / TEMPERATURE).
! At fixed (eps, T), stress is independent of increment history. Global equilibrium still requires convergence.
! Checked in Abaqus for the 3D sphere (abaqus/inputs/sphere.inp, built by scripts/make_total.py):
!   A one step without trapping step, B trapping step + 100 increments, C trapping step + 1 increment
!   -> P_inc = 0.8409 GPa in all three (EosFit, non-linear isotropic: 0.842 GPa).
!
! EoS parameters: EosFit files olivine_angel2018 / diamond_angel2015 (Angel et al. 2018; Angel et al. 2015),
! reference T 298 K:
!   isothermal part: Birch-Murnaghan 3rd order  P_BM3(V; V0, K0, K')
!   DIAMOND (id 1): + Holland-Powell Einstein thermal pressure (independent of V)
!       Pth = alpha0 K0 thetaE/xi0 [ 1/(exp(thetaE/T)-1) - 1/(exp(thetaE/298)-1) ],  xi0 = u0^2 e^u0/(e^u0-1)^2
!   OLIVINE (id 2): + Mie-Grueneisen-Debye thermal pressure
!       Pth = gamma(V)/V [E(T) - E(298)],  E = 3 n R T D3(theta/T),
!       gamma = gamma0 (V/V0)^q,  theta = theta0 exp[(gamma0 - gamma)/q]
!
! UNITS: model in TPa (stress), micron; the EoS works in GPa internally.
! PROPS (new: the last two entries give the reference state where the strain is zero):
!   diamond (constants=10): G, 1, K0, K', alpha0, thetaE, gamma0, q, P_ref, T_ref
!            (gamma0 and q are read but NOT used for diamond: its thermal pressure does not depend on V)
!   olivine (constants=11): G, 2, K0, K', V0, thetaD, n, gamma0, q, P_ref, T_ref
!   G, K0, P_ref in TPa (as the other model inputs); alpha0 in 1/K; V0 in cm3/mol; temperatures in K.
! STATEV: (1) P (TPa)  (2) T (K)  (3) V/V0.   *Depvar 3 (or 2: then only 1,2 are written).
! NOT done: SSE/SPD/SCD, DDSDDT, RPL, DRPLDE and DRPLDT are not assigned (no coupled heat transfer);
! P is not clipped; the EoS is evaluated at any V (no range check). G is constant (no P or T dependence).
! K_T in DDSDDE is a central finite difference of P_EoS (relative step 1e-6).
! =====================================================================================================
module eos_total_mod
    implicit none
    real(8), parameter :: RGAS = 8.314462618d0, T298 = 298.d0
    real(8), save :: xref_c(2) = -1.d0, pref_c(2) = -1.d99, tref_c(2) = -1.d99
contains
    real(8) function debye3(x)
        real(8), intent(in) :: x
        real(8), parameter :: gp(20) = (/ -0.9931285991850949d0, -0.9639719272779138d0, -0.9122344282513259d0, &
            -0.8391169718222188d0, -0.7463319064601508d0, -0.6360536807265150d0, -0.5108670019508271d0, &
            -0.3737060887154196d0, -0.2277858511416451d0, -0.0765265211334973d0, 0.0765265211334973d0, &
             0.2277858511416451d0,  0.3737060887154196d0,  0.5108670019508271d0, 0.6360536807265150d0, &
             0.7463319064601508d0,  0.8391169718222188d0,  0.9122344282513259d0, 0.9639719272779138d0, &
             0.9931285991850949d0 /)
        real(8), parameter :: gw(20) = (/ 0.0176140071391521d0, 0.0406014298003869d0, 0.0626720483341091d0, &
             0.0832767415767048d0, 0.1019301198172404d0, 0.1181945319615184d0, 0.1316886384491766d0, &
             0.1420961093183821d0, 0.1491729864726037d0, 0.1527533871307258d0, 0.1527533871307258d0, &
             0.1491729864726037d0, 0.1420961093183821d0, 0.1316886384491766d0, 0.1181945319615184d0, &
             0.1019301198172404d0, 0.0832767415767048d0, 0.0626720483341091d0, 0.0406014298003869d0, &
             0.0176140071391521d0 /)
        real(8) :: s, t
        integer :: i
        if (x .lt. 1.d-10) then
            debye3 = 1.d0
            return
        end if
        s = 0.d0
        do i = 1, 20
            t = 0.5d0*x*(gp(i) + 1.d0)
            s = s + gw(i)*t**3/(exp(t) - 1.d0)
        end do
        debye3 = 3.d0/x**3*s*0.5d0*x
    end function debye3

    ! pressure (GPa) at x = V/V0 and T; p(:) = PROPS (moduli in TPa, converted here)
    real(8) function peos(id, x, T, p)
        integer, intent(in) :: id
        real(8), intent(in) :: x, T, p(*)
        real(8) :: K0, Kp, f, pbm, u0, xi0, g, th, e1, e0
        K0 = 1000.d0*p(3)
        Kp = p(4)
        f = 0.5d0*(x**(-2.d0/3.d0) - 1.d0)
        pbm = 3.d0*K0*f*(1.d0 + 2.d0*f)**2.5d0*(1.d0 + 1.5d0*(Kp - 4.d0)*f)
        if (id .eq. 1) then
            u0 = p(6)/T298
            xi0 = u0**2*exp(u0)/(exp(u0) - 1.d0)**2
            peos = pbm + p(5)*K0*p(6)/xi0*(1.d0/(exp(p(6)/T) - 1.d0) - 1.d0/(exp(p(6)/T298) - 1.d0))
        else
            g = p(8)*x**p(9)
            th = p(6)*exp((p(8) - g)/p(9))
            e1 = 3.d0*p(7)*RGAS*T*debye3(th/T)
            e0 = 3.d0*p(7)*RGAS*T298*debye3(th/T298)
            peos = pbm + g/(x*p(5))*(e1 - e0)/1000.d0          ! J/cm3 = MPa -> GPa
        end if
    end function peos

    ! isothermal bulk modulus K_T = -V dP/dV (GPa)
    real(8) function kteos(id, x, T, p)
        integer, intent(in) :: id
        real(8), intent(in) :: x, T, p(*)
        real(8), parameter :: h = 1.d-6
        kteos = -x*(peos(id, x*(1.d0 + h), T, p) - peos(id, x*(1.d0 - h), T, p))/(2.d0*h*x)
    end function kteos

    ! V/V0 at the reference state (P_ref in GPa, T_ref): Newton on P(x) = P_ref
    real(8) function xref_of(id, pref, tref, p)
        integer, intent(in) :: id
        real(8), intent(in) :: pref, tref, p(*)
        real(8) :: x, dx
        integer :: it
        x = 1.d0
        do it = 1, 60
            dx = (peos(id, x, tref, p) - pref)*x/kteos(id, x, tref, p)
            x = x + dx
            if (abs(dx) .lt. 1.d-15) exit
        end do
        xref_of = x
    end function xref_of
end module eos_total_mod


subroutine UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD, &
    RPL,DDSDDT,DRPLDE,DRPLDT, &
    STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME, &
    NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT, &
    CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
    use eos_total_mod
    implicit none
    character*80 :: CMNAME
    integer :: NDI,NSHR,NTENS,NSTATV,NPROPS,NOEL,NPT,LAYER,KSPT,KINC
    integer :: JSTEP(4)
    real(8) :: STRESS(NTENS),STATEV(NSTATV)
    real(8) :: DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS)
    real(8) :: STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1)
    real(8) :: PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
    real(8) :: SSE,SPD,SCD,RPL,DRPLDT,DTIME,TEMP,DTEMP,PNEWDT,CELENT

    integer :: id, i, j
    real(8) :: G, pref, tref, x0, eps(6), tr, x, T, P, K, lam

    id = nint(PROPS(2))
    if ((id .ne. 1 .and. id .ne. 2) .or. (id .eq. 1 .and. NPROPS .lt. 10) .or. (id .eq. 2 .and. NPROPS .lt. 11)) then
        write(6,*) 'UMAT_EOS_TOTAL: wrong PROPS (need 10 for diamond, 11 for olivine incl. P_ref, T_ref)', NPROPS
        call XIT
    end if
    G = PROPS(1)                                   ! TPa
    pref = 1000.d0*PROPS(NPROPS - 1)               ! GPa
    tref = PROPS(NPROPS)

    ! reference volume (cached per mineral; recomputed only if P_ref/T_ref change)
    ! Thread safety is not validated. Cache does not include the complete EoS parameter set.
    ! Use one parameter set per mineral and serial execution for initial verification.
    x0 = xref_c(id)
    if (pref .ne. pref_c(id) .or. tref .ne. tref_c(id) .or. x0 .le. 0.d0) then
        x0 = xref_of(id, pref, tref, PROPS)
        xref_c(id) = x0
        pref_c(id) = pref
        tref_c(id) = tref
    end if

    eps = 0.d0
    do i = 1, NTENS
        eps(i) = STRAN(i) + DSTRAN(i)
    end do
    tr = 0.d0
    do i = 1, NDI
        tr = tr + eps(i)
    end do
    x = x0*exp(tr)
    T = TEMP + DTEMP
    P = 1.d-3*peos(id, x, T, PROPS)                ! TPa
    K = 1.d-3*kteos(id, x, T, PROPS)               ! TPa

    do i = 1, NDI
        STRESS(i) = -P + 2.d0*G*(eps(i) - tr/3.d0)
    end do
    do i = NDI + 1, NTENS
        STRESS(i) = G*eps(i)                       ! engineering shear strain
    end do

    lam = K - 2.d0*G/3.d0
    DDSDDE = 0.d0
    do i = 1, NDI
        do j = 1, NDI
            DDSDDE(i,j) = lam
        end do
        DDSDDE(i,i) = lam + 2.d0*G
    end do
    do i = NDI + 1, NTENS
        DDSDDE(i,i) = G
    end do

    if (NSTATV .ge. 1) STATEV(1) = P
    if (NSTATV .ge. 2) STATEV(2) = T
    if (NSTATV .ge. 3) STATEV(3) = x
    return
end subroutine UMAT
