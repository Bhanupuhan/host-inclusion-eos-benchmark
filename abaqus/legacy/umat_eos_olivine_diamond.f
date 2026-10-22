C =====================================================================
C UMAT: incremental (hypoelastic) EoS law for olivine-in-diamond
C
C *** SUPERSEDED by ../umat_eos_total.f90. Kept to document why. ***
C *** Comments corrected 2026-10 for host-inclusion-eos-benchmark;  ***
C *** code unchanged (see legacy/README.md).                        ***
C
C This replaced the lookup-table approach (*Elastic, dependencies + USDFLD)
C with K and alpha computed in the subroutine.
C
C WHAT IT DOES (explicit update, K and alpha at the START of the increment):
C   sigma_new = sigma_old + C(K,G) : (d_eps - alpha/3 dT I)
C DIAMOND: K = K0 + K'P (Murnaghan-type); alpha = alpha0 xi(T)/xi(298)
C          (Einstein); no thermal pressure, no K(T)
C OLIVINE: K = K0 + K'P; alpha = gamma Cv/(K V), V from Murnaghan V(P),
C          Debye Cv; no thermal pressure, no K(T)
C G is constant for both phases.
C NOTE: this law is PATH DEPENDENT and INCREMENT-SIZE DEPENDENT
C (3D sphere: 0.858 GPa without vs 0.780 GPa with the trapping step).
C P stored in STATEV(1) is clipped to >= 0. SSE is not computed.
C UNITS: model units (TPa for G, K0, P); no conversion except x1e12
C inside the olivine alpha.
C
C =====================================================================
C
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     4 JSTEP(4)
C
C     ---- MATERIAL PROPERTIES (passed via *USER MATERIAL) ----
C     PROPS(1)  = G (shear modulus, TPa) -- CONSTANT
C     PROPS(2)  = mineral_id (1=diamond, 2=olivine)
C     
C     For diamond (mineral_id=1):
C       PROPS(3) = K0 = 0.444 TPa (444 GPa)
C       PROPS(4) = Kp = 4.0
C       PROPS(5) = alpha0 = 0.26720e-5 K^-1
C       PROPS(6) = theta_E = 1500.0 K
C       PROPS(7) = gamma0 = 1.0   (read, not used)
C       PROPS(8) = q = 0.0        (read, not used)
C
C     For olivine (mineral_id=2):
C       PROPS(3) = K0 = 0.126302 TPa (126.302 GPa)
C       PROPS(4) = Kp = 4.54207
C       PROPS(5) = V0 = 43.89 cm3/mol
C       PROPS(6) = theta_D0 = 643.364 K
C       PROPS(7) = n = 7
C       PROPS(8) = gamma0 = 1.04424
C       PROPS(9) = q = 1.87509
C
C     STATE VARIABLES:
C       STATEV(1) = current hydrostatic pressure (TPa), clipped >= 0
C       STATEV(2) = current temperature (K)
C
      DOUBLE PRECISION G_CONST, K_CURRENT, ALPHA_CURRENT
      DOUBLE PRECISION P_CURRENT, T_CURRENT, T_NEW
      DOUBLE PRECISION DVOL, DEPS_TH, DEPS_MECH_VOL
      DOUBLE PRECISION ELAM, EMU
      INTEGER MINERAL_ID, I, J
C     (DVOL, DEPS_MECH_VOL are declared but not used)
C
C     ---- Get constant shear modulus ----
      G_CONST = PROPS(1)     ! TPa (model units)
      MINERAL_ID = NINT(PROPS(2))
C
C     ---- Current state ----
      P_CURRENT = STATEV(1)
      T_CURRENT = TEMP        ! Abaqus passes current temperature
      T_NEW     = TEMP + DTEMP
C
C     ---- K(P) and alpha(P,T) at the START of the increment ----
      IF (MINERAL_ID .EQ. 1) THEN
C       Diamond: K = K0 + K'P, Einstein alpha(T)
        CALL DIAMOND_EOS(P_CURRENT, T_CURRENT,
     1                   PROPS(3), PROPS(4), PROPS(5),
     2                   PROPS(6), PROPS(7), PROPS(8),
     3                   K_CURRENT, ALPHA_CURRENT)
      ELSE IF (MINERAL_ID .EQ. 2) THEN
C       Olivine: K = K0 + K'P, alpha = gamma Cv/(K V)
        CALL OLIVINE_EOS(P_CURRENT, T_CURRENT,
     1                   PROPS(3), PROPS(4), PROPS(5),
     2                   PROPS(6), NINT(PROPS(7)), PROPS(8), PROPS(9),
     3                   K_CURRENT, ALPHA_CURRENT)
      END IF
C
C     (no unit conversion: PROPS and STRESS are in model units, TPa)
C
C     ---- Build tangent stiffness DDSDDE ----
C     Isotropic: C_ijkl = (K - 2G/3)*delta_ij*delta_kl + 2G*I_ijkl
      ELAM = K_CURRENT - 2.0D0*G_CONST/3.0D0  ! Lame's lambda
      EMU  = G_CONST                            ! Lame's mu = G
C
C     Initialize
      DO I = 1, NTENS
        DO J = 1, NTENS
          DDSDDE(I,J) = 0.0D0
        END DO
      END DO
C
C     Normal components
      DO I = 1, NDI
        DO J = 1, NDI
          DDSDDE(I,J) = ELAM
        END DO
        DDSDDE(I,I) = ELAM + 2.0D0*EMU
      END DO
C
C     Shear components
      DO I = NDI+1, NTENS
        DDSDDE(I,I) = EMU
      END DO
C
C     ---- Thermal strain increment ----
C     alpha is computed in the subroutine, not from a table
C     Volumetric thermal strain increment = alpha_vol * dT
C     Linear thermal strain increment = alpha_vol/3 * dT
      DEPS_TH = ALPHA_CURRENT / 3.0D0 * DTEMP  ! linear component
C
C     ---- Update stress ----
C     sigma_new = sigma_old + C : (d_epsilon - d_epsilon_thermal)
      DO I = 1, NDI
        DO J = 1, NDI
          STRESS(I) = STRESS(I) + DDSDDE(I,J)*(DSTRAN(J) - DEPS_TH)
        END DO
      END DO
      DO I = NDI+1, NTENS
        STRESS(I) = STRESS(I) + DDSDDE(I,I)*DSTRAN(I)
      END DO
C
C     ---- Update state variables ----
      P_CURRENT = -(STRESS(1)+STRESS(2)+STRESS(3))/3.0D0
      IF (P_CURRENT .LT. 0.0D0) P_CURRENT = 0.0D0
      STATEV(1) = P_CURRENT
      STATEV(2) = T_NEW
C
      RETURN
      END
C
C =====================================================================
C DIAMOND_EOS: K and alpha at given P,T
C K = K0 + K'P (Murnaghan-type); alpha = alpha0 xi(T)/xi(298) (Einstein)
C GAMMA0, Q_PAR, NE_REF, NE_T, F_BM, V_RATIO are not used
C =====================================================================
      SUBROUTINE DIAMOND_EOS(P, T, K0, KP, ALPHA0, THETA_E,
     1                        GAMMA0, Q_PAR, K_OUT, ALPHA_OUT)
      IMPLICIT NONE
      DOUBLE PRECISION P, T, K0, KP, ALPHA0, THETA_E
      DOUBLE PRECISION GAMMA0, Q_PAR, K_OUT, ALPHA_OUT
      DOUBLE PRECISION TREF, U0, U_T, XI0, NE_REF, NE_T
      DOUBLE PRECISION F_BM, V_RATIO
C
      TREF = 298.0D0
C
C     ---- Bulk modulus at current P (no T dependence) ----
C     K(P) = K0 + Kp*P  (Murnaghan-type linear law, NOT BM3)
      K_OUT = K0 + KP * P
C
C     ---- Thermal expansion from HP Einstein model ----
C     alpha(T) = alpha0 * xi(T)/xi(Tref)
C     where xi(T) = (theta/T)^2 * exp(theta/T) / (exp(theta/T)-1)^2
      U0 = THETA_E / TREF
      U_T = THETA_E / T
      XI0 = U0**2 * DEXP(U0) / (DEXP(U0) - 1.0D0)**2
      ALPHA_OUT = ALPHA0 * (U_T**2 * DEXP(U_T) / 
     1            (DEXP(U_T) - 1.0D0)**2) / XI0
C
C     Note: alpha_out is VOLUMETRIC, instantaneous
C     The UMAT integrates deps_th = alpha(T_start)*dT explicitly per
C     increment; the result depends on the increment size and only
C     approaches the integral of alpha(T)dT for small increments
C
      RETURN
      END
C
C =====================================================================
C OLIVINE_EOS: K and alpha at given P,T
C K = K0 + K'P; alpha = gamma(V) Cv(T,V)/(K(P) V), Murnaghan V(P)
C (D3_T, D3_TREF, CV_TREF are declared but not used)
C =====================================================================
      SUBROUTINE OLIVINE_EOS(P, T, K0, KP, V0_CM3, THETA0,
     1                        N_ATOMS, GAMMA0, Q_PAR,
     2                        K_OUT, ALPHA_OUT)
      IMPLICIT NONE
      DOUBLE PRECISION P, T, K0, KP, V0_CM3, THETA0
      INTEGER N_ATOMS
      DOUBLE PRECISION GAMMA0, Q_PAR, K_OUT, ALPHA_OUT
      DOUBLE PRECISION TREF, R_GAS
      DOUBLE PRECISION GAMMA_V, THETA_V, V_RATIO
      DOUBLE PRECISION D3_T, D3_TREF, CV_T, CV_TREF
      DOUBLE PRECISION DEBYE_D3
      EXTERNAL DEBYE_D3
C
      TREF = 298.0D0
      R_GAS = 8.314462618D0
C
C     ---- Bulk modulus at current P ----
      K_OUT = K0 + KP * P
C
C     ---- Thermal expansion from MGD ----
C     Need V/V0 at current P. From Murnaghan:
C     V/V0 = (1 + Kp*P/K0)^(-1/Kp)
      V_RATIO = (1.0D0 + KP*P/K0)**(-1.0D0/KP)
C
C     gamma(V) = gamma0 * (V/V0)^q
      GAMMA_V = GAMMA0 * V_RATIO**Q_PAR
C
C     theta(V) = theta0 * exp[(gamma0 - gamma(V))/q]
      IF (DABS(Q_PAR) .GT. 1.0D-10) THEN
        THETA_V = THETA0 * DEXP((GAMMA0 - GAMMA_V)/Q_PAR)
      ELSE
        THETA_V = THETA0 * V_RATIO**(-GAMMA0)
      END IF
C
C     Cv from Debye model: Cv = 3nR * 4D3(x) - 3x/(e^x - 1)
C     Actually, alpha = gamma * Cv / (K * V)
C     Cv = 3nR * [4*D3(theta/T) - 3*(theta/T)/(exp(theta/T)-1)]
C     Simpler: Cv = d(3nRT*D3)/dT
C     alpha_vol = gamma * Cv / (K * V)
C
C     For the thermal expansion:
C     alpha = gamma(V) * Cv(T,V) / (K(P) * V)
C     where Cv = 3*n*R * [4*D3(x) - 3*x/(e^x - 1)], x = theta/T
C
      CALL COMPUTE_CV_DEBYE(T, THETA_V, N_ATOMS, R_GAS, CV_T)
C
C     alpha_vol = gamma * Cv / (K * V_molar_m3)
C     V_molar_m3 = V0_cm3 * V_ratio * 1e-6
      ALPHA_OUT = GAMMA_V * CV_T / 
     1            (K_OUT * 1.0D12 * V0_CM3 * V_RATIO * 1.0D-6)
C
C     ALPHA_OUT is volumetric, instantaneous, in K^-1
C
      RETURN
      END
C
C =====================================================================
C COMPUTE_CV_DEBYE: Heat capacity from Debye model
C Cv = 3nR * [4*D3(x) - 3x/(e^x - 1)],  x = theta/T
C =====================================================================
      SUBROUTINE COMPUTE_CV_DEBYE(T, THETA, N_ATOMS, R_GAS, CV)
      IMPLICIT NONE
      DOUBLE PRECISION T, THETA, R_GAS, CV
      INTEGER N_ATOMS
      DOUBLE PRECISION X, D3_VAL
      DOUBLE PRECISION DEBYE_D3
      EXTERNAL DEBYE_D3
C
      IF (T .LT. 1.0D0) THEN
        CV = 0.0D0
        RETURN
      END IF
C
      X = THETA / T
      D3_VAL = DEBYE_D3(X)
C
C     Cv = 3nR * [4*D3(x) - 3*x/(exp(x)-1)]
      CV = 3.0D0 * DBLE(N_ATOMS) * R_GAS * 
     1     (4.0D0*D3_VAL - 3.0D0*X/(DEXP(X) - 1.0D0))
C
C     CV in J/(mol*K)
      RETURN
      END
C
C =====================================================================
C DEBYE_D3: Debye function D3(x) = (3/x^3) * int_0^x t^3/(e^t-1) dt
C Numerical integration using Gauss-Legendre quadrature
C =====================================================================
      DOUBLE PRECISION FUNCTION DEBYE_D3(X)
      IMPLICIT NONE
      DOUBLE PRECISION X
      DOUBLE PRECISION SUM_VAL, T_VAL, W_VAL, HALF_X
C     (W_VAL is declared but not used)
      INTEGER I, N_GAUSS
      PARAMETER (N_GAUSS = 20)
C     Gauss-Legendre points and weights on [-1,1] (20-point)
      DOUBLE PRECISION GP(20), GW(20)
      DATA GP / -0.9931285991850949D0, -0.9639719272779138D0,
     1          -0.9122344282513259D0, -0.8391169718222188D0,
     2          -0.7463319064601508D0, -0.6360536807265150D0,
     3          -0.5108670019508271D0, -0.3737060887154196D0,
     4          -0.2277858511416451D0, -0.0765265211334973D0,
     5           0.0765265211334973D0,  0.2277858511416451D0,
     6           0.3737060887154196D0,  0.5108670019508271D0,
     7           0.6360536807265150D0,  0.7463319064601508D0,
     8           0.8391169718222188D0,  0.9122344282513259D0,
     9           0.9639719272779138D0,  0.9931285991850949D0 /
      DATA GW /  0.0176140071391521D0,  0.0406014298003869D0,
     1           0.0626720483341091D0,  0.0832767415767048D0,
     2           0.1019301198172404D0,  0.1181945319615184D0,
     3           0.1316886384491766D0,  0.1420961093183821D0,
     4           0.1491729864726037D0,  0.1527533871307258D0,
     5           0.1527533871307258D0,  0.1491729864726037D0,
     6           0.1420961093183821D0,  0.1316886384491766D0,
     7           0.1181945319615184D0,  0.1019301198172404D0,
     8           0.0832767415767048D0,  0.0626720483341091D0,
     9           0.0406014298003869D0,  0.0176140071391521D0 /
C
      IF (X .LT. 1.0D-12) THEN
        DEBYE_D3 = 1.0D0
        RETURN
      END IF
      IF (X .GT. 150.0D0) THEN
        DEBYE_D3 = 0.0D0
        RETURN
      END IF
C
C     Transform integral from [0,x] to [-1,1]
      HALF_X = X / 2.0D0
      SUM_VAL = 0.0D0
      DO I = 1, N_GAUSS
        T_VAL = HALF_X * (GP(I) + 1.0D0)  ! map [-1,1] to [0,x]
        IF (T_VAL .GT. 1.0D-12) THEN
          SUM_VAL = SUM_VAL + GW(I) * T_VAL**3 / (DEXP(T_VAL)-1.0D0)
        END IF
      END DO
      SUM_VAL = SUM_VAL * HALF_X  ! Jacobian of transformation
C
      DEBYE_D3 = 3.0D0 / X**3 * SUM_VAL
C
      RETURN
      END
