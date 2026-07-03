!==============================================================================
! test_puck.f90  –  Comprehensive unit tests for the Puck failure criterion
!                   Extracted from UMAT_fix_anstrengung_260624.for
!
! Compile: ifx test_puck.f90 -o test_puck.exe
! Run:     test_puck.exe
!
! Test groups:
!   Fiber Fracture (FF):  Tests 1-4   (tension, compression, biaxial, post-failure)
!   IFF Mode A:           Tests 5-9   (matrix tension, shear, out-of-plane, combined)
!   IFF Mode B:           Test  10    (inclined compression + shear at theta=0)
!   IFF Mode C:           Test  11    (transverse compression, self-consistency)
!   Edge cases:           Tests 12-14 (sub-critical, independence, zero stress)
!==============================================================================
program test_puck

    implicit none

    real(8), parameter :: Pi = 3.14159265358979d+0

    ! ---- Material properties (typical Carbon / Epoxy) ----
    real(8) :: Xt, Xc, Yt, Yc, Zt, Zc, S12, S13, S23
    real(8) :: E11, E11f, ANU12, ANU12f, MGF, THETAF
    real(8) :: EpsXt, EpsXc, EpsYt, EpsYc, EpsZt, EpsZc
    real(8) :: GamS12, GamS13, GamS23

    real(8) :: STRESS(6), UPSTRAN(6), e(6), angles_out(6)
    integer :: failure_id
    integer :: n_pass, n_fail
    ! ---- Damage variables ----
    real(8) :: dmg(6)
    integer :: fflags(6), damage_id, KINC
    real(8) :: beta_ft, beta_fc, beta_mt, beta_mc, beta_s
    real(8) :: a_ft, a_fc, a_mt, a_mc, a_s, n_ft, n_fc, n_mt, n_mc, n_s
    real(8) :: G_ft, G_fc, G_mt, G_mc, G_IIC, le, let, alpha, beta_p, E22, E33, G12, G13, G23


    ! Derived Puck parameters (for computing expected values)
    real(8) :: Pnt_ref, Rnt_ref, Pn1_ref, coeff_MGF
    real(8) :: sigma2_B, tau12_B, e2_comb_exp, e1_biax_exp

    ! ---- Material setup ----
    Xt     = 2000.0d0    ! MPa – fiber tensile strength
    Xc     = 1200.0d0    ! MPa – fiber compressive strength
    Yt     =   50.0d0    ! MPa – matrix tensile strength (= Rn)
    Yc     =  250.0d0    ! MPa – matrix compressive strength
    Zt     =   50.0d0    ! MPa – out-of-plane tensile strength
    Zc     =  250.0d0    ! MPa
    S12    =   70.0d0    ! MPa – in-plane shear strength (= Rn1)
    S13    =   70.0d0    ! MPa – out-of-plane shear strength 13
    S23    =   50.0d0    ! MPa – transverse shear strength (= Rn in plane)
    E11    = 135000.0d0  ! MPa – fiber direction modulus
    E11f   = 230000.0d0  ! MPa – pure fiber modulus
    ANU12  = 0.27d0      ! composite Poisson ratio
    ANU12f = 0.20d0      ! fiber Poisson ratio
    MGF    = 1.1d0       ! Puck magnification factor
    THETAF = 53.0d0 * Pi / 180.0d0   ! fracture angle (rad)

    EpsXt  = 0.015d0;  EpsXc  = 0.010d0
    EpsYt  = 0.005d0;  EpsYc  = 0.025d0
    EpsZt  = 0.005d0;  EpsZc  = 0.025d0
    GamS12 = 0.020d0;  GamS13 = 0.020d0;  GamS23 = 0.020d0

    beta_ft = 0.1d0; beta_fc = 0.1d0; beta_mt = 0.1d0; beta_mc = 0.1d0; beta_s = 0.1d0
    a_ft = 1.0d0; a_fc = 1.0d0; a_mt = 1.0d0; a_mc = 1.0d0; a_s = 1.0d0
    n_ft = 1.0d0; n_fc = 1.0d0; n_mt = 1.0d0; n_mc = 1.0d0; n_s = 1.0d0
    G_ft=0d0; G_fc=0d0; G_mt=0d0; G_mc=0d0; G_IIC=0d0; le=0d0; let=0d0; alpha=0d0; beta_p=0d0
    E22=1d0; E33=1d0; G12=1d0; G13=1d0; G23=1d0
    KINC = 1


    failure_id = 7   ! Puck
    n_pass = 0;  n_fail = 0

    ! Precompute Puck parameters analytically (for expected-value formulas)
    Pnt_ref  = -1.0d0 / tan(2.0d0*THETAF)
    Rnt_ref  = Yc / (2.0d0*tan(THETAF))
    Pn1_ref  = Pnt_ref * (S12 / Rnt_ref)
    coeff_MGF = ANU12 - ANU12f*MGF*(E11/E11f)   ! used in FF formula

    write(*,'(/,A)') "================================================================"
    write(*,'(A)')   "  Puck Failure Criterion  –  Comprehensive Unit Tests"
    write(*,'(A,F7.1,A,F7.1,A,F7.1,A,F7.1,A,F7.1)') &
        "  Xt=",Xt,"  Xc=",Xc,"  Yt=",Yt,"  Yc=",Yc,"  S12=",S12
    write(*,'(A,F7.2,A,F6.3,A,F6.3)') &
        "  THETAF=",THETAF*180.d0/Pi,"deg  Pnt=",Pnt_ref,"  Rnt=",Rnt_ref
    write(*,'(A,/)') "================================================================"

    !==================================================================
    !  GROUP 1: FIBER FRACTURE (FF)
    !==================================================================
    write(*,'(A,/)') "--- Group 1: Fiber Fracture (FF) ---"

    ! Test 1: Pure fiber tension at exactly Xt
    STRESS = 0.0d0;  UPSTRAN = 0.0d0;  e = 0.0d0;  angles_out = 0.0d0
    STRESS(1) = Xt
    call failure_calc(failure_id,STRESS,UPSTRAN,e,Xt,Xc,Yt,Yc,Zt,Zc,S12,S13,S23, &
        EpsXt,EpsXc,EpsYt,EpsYc,EpsZt,EpsZc,GamS12,GamS13,GamS23, &
        THETAF,MGF,ANU12,ANU12f,E11,E11f,angles_out)
    write(*,'(A)') "Test 1: Pure fiber tension  (sigma_1 = Xt)"
    write(*,'(A,F9.5)') "  e(1)=",e(1)
    call check("e(1) = 1.0",         e(1), 1.0d0, 0.001d0, n_pass, n_fail)
    call check("e(2) = 0 (no IFF)",  e(2), 0.0d0, 0.001d0, n_pass, n_fail)
    write(*,*)

    ! Test 2: Pure fiber compression at exactly -Xc
    !         [Key: verifies sign fix – must return positive index]
    STRESS = 0.0d0;  e = 0.0d0;  angles_out = 0.0d0
    STRESS(1) = -Xc
    call failure_calc(failure_id,STRESS,UPSTRAN,e,Xt,Xc,Yt,Yc,Zt,Zc,S12,S13,S23, &
        EpsXt,EpsXc,EpsYt,EpsYc,EpsZt,EpsZc,GamS12,GamS13,GamS23, &
        THETAF,MGF,ANU12,ANU12f,E11,E11f,angles_out)
    write(*,'(A)') "Test 2: Pure fiber compression  (sigma_1 = -Xc)  [sign-fix verify]"
    write(*,'(A,F9.5)') "  e(1)=",e(1)
    call check("e(1) = 1.0 (positive)", e(1), 1.0d0, 0.001d0, n_pass, n_fail)
    write(*,*)

    ! Test 3: Biaxial – fiber tension + transverse stress (MGF cross-coupling)
    !         A = sigma_1 - coeff_MGF*(sigma_2+sigma_3)
    !         sigma_2 chosen so that coeff_MGF*sigma_2 = 0.05*Xt => e(1)=0.95
    STRESS = 0.0d0;  e = 0.0d0;  angles_out = 0.0d0
    STRESS(1) = Xt
    STRESS(2) = 0.05d0*Xt / coeff_MGF   ! ~710 MPa for these params
    e1_biax_exp = (Xt - coeff_MGF*(STRESS(2))) / Xt   ! = 0.95 exactly
    call failure_calc(failure_id,STRESS,UPSTRAN,e,Xt,Xc,Yt,Yc,Zt,Zc,S12,S13,S23, &
        EpsXt,EpsXc,EpsYt,EpsYc,EpsZt,EpsZc,GamS12,GamS13,GamS23, &
        THETAF,MGF,ANU12,ANU12f,E11,E11f,angles_out)
    write(*,'(A,F6.1,A)') "Test 3: Biaxial FF  (sigma_1=Xt, sigma_2=",STRESS(2)," MPa)  [MGF coupling]"
    write(*,'(A,F9.5,A,F9.5)') "  e(1)=",e(1),"  expected=",e1_biax_exp
    call check("e(1) = 0.95 (biaxial MGF)", e(1), e1_biax_exp, 0.001d0, n_pass, n_fail)
    write(*,*)

    ! Test 4: Post-failure – double fiber tension (e > 1)
    STRESS = 0.0d0;  e = 0.0d0;  angles_out = 0.0d0
    STRESS(1) = 2.0d0*Xt
    call failure_calc(failure_id,STRESS,UPSTRAN,e,Xt,Xc,Yt,Yc,Zt,Zc,S12,S13,S23, &
        EpsXt,EpsXc,EpsYt,EpsYc,EpsZt,EpsZc,GamS12,GamS13,GamS23, &
        THETAF,MGF,ANU12,ANU12f,E11,E11f,angles_out)
    write(*,'(A)') "Test 4: Post-failure fiber tension  (sigma_1 = 2*Xt)  [e > 1]"
    write(*,'(A,F9.5)') "  e(1)=",e(1)
    call check("e(1) = 2.0", e(1), 2.0d0, 0.001d0, n_pass, n_fail)
    write(*,*)

    !==================================================================
    !  GROUP 2: IFF MODE A
    !==================================================================
    write(*,'(A,/)') "--- Group 2: IFF Mode A ---"

    ! Test 5: Pure matrix tension at Yt  (standard Mode A at theta=0)
    STRESS = 0.0d0;  e = 0.0d0;  angles_out = 0.0d0
    STRESS(2) = Yt
    call failure_calc(failure_id,STRESS,UPSTRAN,e,Xt,Xc,Yt,Yc,Zt,Zc,S12,S13,S23, &
        EpsXt,EpsXc,EpsYt,EpsYc,EpsZt,EpsZc,GamS12,GamS13,GamS23, &
        THETAF,MGF,ANU12,ANU12f,E11,E11f,angles_out)
    write(*,'(A)') "Test 5: Pure matrix tension  (sigma_2=Yt)  – Mode A, theta=0"
    write(*,'(A,F9.5,A,F7.2,A)') "  e(2)=",e(2),"  angle=",angles_out(2)," deg"
    call check("e(2) = 1.0",       e(2),        1.0d0, 0.001d0, n_pass, n_fail)
    call check("angle = 0 deg",    angles_out(2), 0.0d0, 1.0d0,  n_pass, n_fail)
    write(*,*)

    ! Test 6: Pure in-plane shear at S12  (Mode A at theta=0, sigma_n=0)
    STRESS = 0.0d0;  e = 0.0d0;  angles_out = 0.0d0
    STRESS(4) = S12
    call failure_calc(failure_id,STRESS,UPSTRAN,e,Xt,Xc,Yt,Yc,Zt,Zc,S12,S13,S23, &
        EpsXt,EpsXc,EpsYt,EpsYc,EpsZt,EpsZc,GamS12,GamS13,GamS23, &
        THETAF,MGF,ANU12,ANU12f,E11,E11f,angles_out)
    write(*,'(A)') "Test 6: Pure in-plane shear  (tau_12=S12)  – Mode A, theta=0"
    write(*,'(A,F9.5,A,F7.2,A)') "  e(2)=",e(2),"  angle=",angles_out(2)," deg"
    call check("e(2) = 1.0",    e(2),        1.0d0, 0.001d0, n_pass, n_fail)
    call check("angle = 0 deg", angles_out(2), 0.0d0, 1.0d0,  n_pass, n_fail)
    write(*,*)

    ! Test 7: Out-of-plane shear tau_13 = S13  (Mode A, fracture at theta=90)
    !         TN1 = tau_12*cos + tau_13*sin -> max at theta=90 -> TN1=tau_13
    STRESS = 0.0d0;  e = 0.0d0;  angles_out = 0.0d0
    STRESS(5) = S13
    call failure_calc(failure_id,STRESS,UPSTRAN,e,Xt,Xc,Yt,Yc,Zt,Zc,S12,S13,S23, &
        EpsXt,EpsXc,EpsYt,EpsYc,EpsZt,EpsZc,GamS12,GamS13,GamS23, &
        THETAF,MGF,ANU12,ANU12f,E11,E11f,angles_out)
    write(*,'(A)') "Test 7: Out-of-plane shear  (tau_13=S13)  – Mode A, theta=90"
    write(*,'(A,F9.5,A,F7.2,A)') "  e(2)=",e(2),"  angle=",angles_out(2)," deg"
    call check(    "e(2) = 1.0",         e(2),        1.0d0,  0.001d0, n_pass, n_fail)
    call check_sym("angle = +-90 deg",   angles_out(2), 90.0d0, 1.0d0,  n_pass, n_fail)
    write(*,*)

    ! Test 8: Transverse shear tau_23 = S23  (Mode A, fracture at theta=45)
    !         SFP = 2*tau23*sin*cos -> max at theta=45: SFP = tau23
    !         IFF = SFP/Rn = S23/Yt = 50/50 = 1.0
    STRESS = 0.0d0;  e = 0.0d0;  angles_out = 0.0d0
    STRESS(6) = S23
    call failure_calc(failure_id,STRESS,UPSTRAN,e,Xt,Xc,Yt,Yc,Zt,Zc,S12,S13,S23, &
        EpsXt,EpsXc,EpsYt,EpsYc,EpsZt,EpsZc,GamS12,GamS13,GamS23, &
        THETAF,MGF,ANU12,ANU12f,E11,E11f,angles_out)
    write(*,'(A)') "Test 8: Transverse shear  (tau_23=S23=50MPa)  – Mode A, theta=45"
    write(*,'(A,F9.5,A,F7.2,A)') "  e(2)=",e(2),"  angle=",angles_out(2)," deg"
    call check("e(2) = 1.0",     e(2),        1.0d0, 0.001d0, n_pass, n_fail)
    call check("angle = 45 deg", angles_out(2), 45.0d0, 1.0d0, n_pass, n_fail)
    write(*,*)

    ! Test 9: Out-of-plane normal sigma_3 = Zt  (Mode A at theta=90, SFP=sigma_3)
    STRESS = 0.0d0;  e = 0.0d0;  angles_out = 0.0d0
    STRESS(3) = Zt
    call failure_calc(failure_id,STRESS,UPSTRAN,e,Xt,Xc,Yt,Yc,Zt,Zc,S12,S13,S23, &
        EpsXt,EpsXc,EpsYt,EpsYc,EpsZt,EpsZc,GamS12,GamS13,GamS23, &
        THETAF,MGF,ANU12,ANU12f,E11,E11f,angles_out)
    write(*,'(A)') "Test 9: Out-of-plane normal  (sigma_3=Zt)  – Mode A, theta=90"
    write(*,'(A,F9.5,A,F7.2,A)') "  e(2)=",e(2),"  angle=",angles_out(2)," deg"
    call check(    "e(2) = 1.0",         e(2),        1.0d0,  0.001d0, n_pass, n_fail)
    call check_sym("angle = +-90 deg",   angles_out(2), 90.0d0, 1.0d0,  n_pass, n_fail)
    write(*,*)

    !==================================================================
    !  GROUP 3: IFF MODE B
    !==================================================================
    write(*,'(A,/)') "--- Group 3: IFF Mode B ---"

    ! Test 10: Compressive sigma_2 + shear tau_12 → failure at theta=0 (Mode B)
    !          Critical tau_12 = Rn1 - Pn1*sigma_2 = S12 + Pn1*|sigma_2|
    !          At theta=0: SFP=sigma_2<0, J=90, IFF = tau12/(S12-Pn1*sigma_2) = 1.0
    sigma2_B  = -50.0d0
    tau12_B   = S12 - Pn1_ref*sigma2_B   ! = S12 + Pn1*50
    STRESS = 0.0d0;  e = 0.0d0;  angles_out = 0.0d0
    STRESS(2) = sigma2_B
    STRESS(4) = tau12_B
    call failure_calc(failure_id,STRESS,UPSTRAN,e,Xt,Xc,Yt,Yc,Zt,Zc,S12,S13,S23, &
        EpsXt,EpsXc,EpsYt,EpsYc,EpsZt,EpsZc,GamS12,GamS13,GamS23, &
        THETAF,MGF,ANU12,ANU12f,E11,E11f,angles_out)
    write(*,'(A,F6.1,A,F6.1,A)') "Test 10: Mode B  (sigma_2=",sigma2_B,"  tau_12=",tau12_B," MPa)"
    write(*,'(A)') "         [IFF maximised at theta=0 under compression+shear]"
    write(*,'(A,F9.5,A,F9.5,A,F7.2,A)') &
        "  e(2)=",e(2),"  e(4)=",e(4),"  angle_B=",angles_out(4)," deg"
    call check("e(2) = 1.0",            e(2),        1.0d0, 0.01d0, n_pass, n_fail)
    call check("e(4) = 1.0 (Mode B)",   e(4),        1.0d0, 0.01d0, n_pass, n_fail)
    call check("Mode B angle = 0 deg",  angles_out(4), 0.0d0, 0.1d0, n_pass, n_fail)
    write(*,*)

    !==================================================================
    !  GROUP 4: IFF MODE C
    !==================================================================
    write(*,'(A,/)') "--- Group 4: IFF Mode C ---"

    ! Test 11: Pure transverse compression sigma_2 = -Yc  (self-consistency)
    !          With Pnt=-1/tan(2*THETAF) and Rnt=Yc/(2*tan(THETAF)):
    !          IFF = 1.0 exactly at theta = +-THETAF
    !          Search finds negative mirror angle first.
    STRESS = 0.0d0;  e = 0.0d0;  angles_out = 0.0d0
    STRESS(2) = -Yc
    call failure_calc(failure_id,STRESS,UPSTRAN,e,Xt,Xc,Yt,Yc,Zt,Zc,S12,S13,S23, &
        EpsXt,EpsXc,EpsYt,EpsYc,EpsZt,EpsZc,GamS12,GamS13,GamS23, &
        THETAF,MGF,ANU12,ANU12f,E11,E11f,angles_out)
    write(*,'(A)') "Test 11: Pure transverse compression  (sigma_2=-Yc)  – Mode C"
    write(*,'(A)') "         [Self-consistency: IFF=1.0 at +-THETAF]"
    write(*,'(A,F9.5,A,F9.5,A,F7.2,A)') &
        "  e(2)=",e(2),"  e(5)=",e(5),"  angle_C=",angles_out(5)," deg"
    call check("e(2) = 1.0",              e(2),        1.0d0,  0.001d0, n_pass, n_fail)
    call check("angle_C = -53 deg",       angles_out(5),-53.0d0, 2.0d0,  n_pass, n_fail)
    write(*,*)

    !==================================================================
    !  GROUP 5: COMBINED AND EDGE CASES
    !==================================================================
    write(*,'(A,/)') "--- Group 5: Combined Stress & Edge Cases ---"

    ! Test 12: Combined Mode A – sigma_2=Yt AND tau_12=S12 simultaneously
    !          Both at failure → IFF > 1 (failure already exceeded)
    !          e2_expected = SQRT(1 + (S12/(S12 - Pn1*Yt))^2) from Mode A formula
    e2_comb_exp = SQRT(1.0d0 + (S12/(S12 - Pn1_ref*Yt))**2)
    STRESS = 0.0d0;  e = 0.0d0;  angles_out = 0.0d0
    STRESS(2) = Yt
    STRESS(4) = S12
    call failure_calc(failure_id,STRESS,UPSTRAN,e,Xt,Xc,Yt,Yc,Zt,Zc,S12,S13,S23, &
        EpsXt,EpsXc,EpsYt,EpsYc,EpsZt,EpsZc,GamS12,GamS13,GamS23, &
        THETAF,MGF,ANU12,ANU12f,E11,E11f,angles_out)
    write(*,'(A)') "Test 12: Combined Mode A  (sigma_2=Yt, tau_12=S12)  [e > 1]"
    write(*,'(A,F9.5,A,F9.5)') "  e(2)=",e(2),"  expected=",e2_comb_exp
    call check("e(2) > 1 (combined)", e(2), e2_comb_exp, 0.02d0, n_pass, n_fail)
    write(*,*)

    ! Test 13: Sub-critical – half matrix tension, check linearity of Mode A
    !          sigma_2 = 0.5*Yt -> e(2)=0.5, e(3)=0.5
    STRESS = 0.0d0;  e = 0.0d0;  angles_out = 0.0d0
    STRESS(2) = 0.5d0*Yt
    call failure_calc(failure_id,STRESS,UPSTRAN,e,Xt,Xc,Yt,Yc,Zt,Zc,S12,S13,S23, &
        EpsXt,EpsXc,EpsYt,EpsYc,EpsZt,EpsZc,GamS12,GamS13,GamS23, &
        THETAF,MGF,ANU12,ANU12f,E11,E11f,angles_out)
    write(*,'(A)') "Test 13: Sub-critical Mode A  (sigma_2=0.5*Yt)  [linearity check]"
    write(*,'(A,F9.5,A,F9.5)') "  e(2)=",e(2),"  e(3)=",e(3)
    call check("e(2) = 0.5", e(2), 0.5d0, 0.001d0, n_pass, n_fail)
    call check("e(3) = 0.5", e(3), 0.5d0, 0.001d0, n_pass, n_fail)
    write(*,*)

    ! Test 14: Zero stress – all failure indices must be exactly 0
    !          [Ensures no garbage values / correct initialization]
    STRESS = 0.0d0;  e = 0.0d0;  angles_out = 0.0d0
    call failure_calc(failure_id,STRESS,UPSTRAN,e,Xt,Xc,Yt,Yc,Zt,Zc,S12,S13,S23, &
        EpsXt,EpsXc,EpsYt,EpsYc,EpsZt,EpsZc,GamS12,GamS13,GamS23, &
        THETAF,MGF,ANU12,ANU12f,E11,E11f,angles_out)
    write(*,'(A)') "Test 14: Zero stress  – all e must be 0  [initialization]"
    write(*,'(A,F8.5,A,F8.5)') "  e(1)=",e(1),"  e(2)=",e(2)
    call check("e(1) = 0 (no stress)", e(1), 0.0d0, 1.0d-10, n_pass, n_fail)
    call check("e(2) = 0 (no stress)", e(2), 0.0d0, 1.0d-10, n_pass, n_fail)
    write(*,*)


    !==================================================================
    !  GROUP 6: DAMAGE DEGRADATION (damage_calc)
    !==================================================================
    write(*,'(A,/)') "--- Group 6: Damage Degradation (damage_calc) ---"
    
    ! Test 15: Instantaneous damage (damage_id=1), Fiber Tension
    damage_id = 1
    STRESS = 0.0d0; UPSTRAN = 0.0d0; e = 0.0d0; dmg = 0.0d0; fflags = 0
    e(1) = 1.5d0    ! Exceeds failure
    STRESS(1) = Xt  ! Tension
    call damage_calc(failure_id,damage_id,STRESS,UPSTRAN,dmg,e,fflags,KINC,beta_ft,beta_fc,beta_mt,beta_mc,beta_s,&
        EpsXt,EpsXc,EpsYt,EpsYc,EpsZt,EpsZc,GamS12,GamS13,GamS23,a_ft,a_fc,a_mt,a_mc,a_s,n_ft,n_fc,n_mt,n_mc,n_s,&
        G_ft,G_fc,G_mt,G_mc,G_IIC,le,let,alpha,beta_p,E11,E22,E33,G12,G13,G23,Xt,Xc,Yt,Yc,Zt,Zc,S12,S13,S23)
    write(*,'(A)') "Test 15: Instantaneous Damage - Fiber Tension (e1=1.5)"
    write(*,'(A,F9.5)') "  dmg(1)=",dmg(1)
    call check("dmg(1) = 0.9 (1.0-beta)", dmg(1), 0.99999999d0-beta_ft, 0.0001d0, n_pass, n_fail)
    call check("fflags(1) = KINC", real(fflags(1),8), real(KINC,8), 0.01d0, n_pass, n_fail)
    write(*,*)

    ! Test 16: Instantaneous damage (damage_id=1), Matrix Tension (Puck induced shear damage)
    STRESS = 0.0d0; UPSTRAN = 0.0d0; e = 0.0d0; dmg = 0.0d0; fflags = 0
    e(2) = 1.2d0    ! Matrix failure
    STRESS(2) = Yt
    call damage_calc(failure_id,damage_id,STRESS,UPSTRAN,dmg,e,fflags,KINC,beta_ft,beta_fc,beta_mt,beta_mc,beta_s,&
        EpsXt,EpsXc,EpsYt,EpsYc,EpsZt,EpsZc,GamS12,GamS13,GamS23,a_ft,a_fc,a_mt,a_mc,a_s,n_ft,n_fc,n_mt,n_mc,n_s,&
        G_ft,G_fc,G_mt,G_mc,G_IIC,le,let,alpha,beta_p,E11,E22,E33,G12,G13,G23,Xt,Xc,Yt,Yc,Zt,Zc,S12,S13,S23)
    write(*,'(A)') "Test 16: Instantaneous Damage - Matrix Tension (Puck induced shear)"
    write(*,'(A,F9.5,A,F9.5,A,F9.5,A,F9.5)') "  dmg(2)=",dmg(2),"  dmg(3)=",dmg(3),"  dmg(4)=",dmg(4),"  dmg(6)=",dmg(6)
    call check("dmg(2) = 0.9", dmg(2), 0.99999999d0-beta_mt, 0.0001d0, n_pass, n_fail)
    call check("dmg(3) = 0.9", dmg(3), 0.99999999d0-beta_mt, 0.0001d0, n_pass, n_fail)
    call check("dmg(4) = 0.9", dmg(4), 0.99999999d0-beta_s, 0.0001d0, n_pass, n_fail)
    call check("dmg(6) = 0.9", dmg(6), 0.99999999d0-beta_s, 0.0001d0, n_pass, n_fail)
    write(*,*)

    ! Test 17: Constant Stress Degradation (damage_id=4), Fiber Compression
    damage_id = 4
    STRESS = 0.0d0; UPSTRAN = 0.0d0; e = 0.0d0; dmg = 1.0d0; fflags = 0  ! dmg initially 1.0
    e(1) = 2.0d0
    STRESS(1) = -Xc
    call damage_calc(failure_id,damage_id,STRESS,UPSTRAN,dmg,e,fflags,KINC,beta_ft,beta_fc,beta_mt,beta_mc,beta_s,&
        EpsXt,EpsXc,EpsYt,EpsYc,EpsZt,EpsZc,GamS12,GamS13,GamS23,a_ft,a_fc,a_mt,a_mc,a_s,n_ft,n_fc,n_mt,n_mc,n_s,&
        G_ft,G_fc,G_mt,G_mc,G_IIC,le,let,alpha,beta_p,E11,E22,E33,G12,G13,G23,Xt,Xc,Yt,Yc,Zt,Zc,S12,S13,S23)
    write(*,'(A)') "Test 17: Constant Stress Degradation - Fiber Compression (e1=2.0)"
    write(*,'(A,F9.5)') "  dmg(1)=",dmg(1)
    call check("dmg(1) = 1.0/2.0 = 0.5", dmg(1), 0.5d0, 0.0001d0, n_pass, n_fail)
    call check("fflags(1) = -KINC", real(fflags(1),8), real(-KINC,8), 0.01d0, n_pass, n_fail)
    write(*,*)

    !==========================================================
    ! Summary
    !==========================================================
    write(*,'(A)') "================================================================"
    write(*,'(A,I3,A,I3,A)') "  Results:  ", n_pass, " PASS  /  ", n_fail, " FAIL"
    if (n_fail == 0) then
        write(*,'(A)') "  >>> ALL TESTS PASSED <<<"
    else
        write(*,'(A)') "  >>> SOME TESTS FAILED – check output above <<<"
    end if
    write(*,'(A,/)') "================================================================"

contains

    subroutine check(label, val, expected, tol, np, nf)
        character(len=*), intent(in) :: label
        real(8), intent(in)          :: val, expected, tol
        integer, intent(inout)       :: np, nf
        if (abs(val - expected) <= tol) then
            write(*,'(A,A)') "  [PASS] ", trim(label)
            np = np + 1
        else
            write(*,'(A,A,A,F12.6,A,F12.6,A,F12.6)') &
                "  [FAIL] ", trim(label), &
                "   got=", val, "  expected=", expected, "  tol=+/-", tol
            nf = nf + 1
        end if
    end subroutine check

    ! check_sym: accepts +/-|expected| (for angles with symmetric solutions)
    subroutine check_sym(label, val, expected_abs, tol, np, nf)
        character(len=*), intent(in) :: label
        real(8), intent(in)          :: val, expected_abs, tol
        integer, intent(inout)       :: np, nf
        if (abs(abs(val) - abs(expected_abs)) <= tol) then
            write(*,'(A,A,A,F7.2,A)') "  [PASS] ", trim(label), "  (got ", val, " deg)"
            np = np + 1
        else
            write(*,'(A,A,A,F10.4,A,F10.4,A,F10.4)') &
                "  [FAIL] ", trim(label), &
                "   got=", val, "  expected=+/-", expected_abs, "  tol=+/-", tol
            nf = nf + 1
        end if
    end subroutine check_sym

end program test_puck


!==============================================================================
! failure_calc – Puck criterion (Case 7)
! Mirror of UMAT_fix_anstrengung_260624.for with all fixes applied:
!   1. FF sign fix: e(1) always positive
!   2. Pnt = -1/tan(2*THETAF)  [was -1/(2*tan), factor-2 fix]
!   3. Rnt = Yc/(2*tan(THETAF))  [self-consistent with Pnt]
!   4. Two-pass angle search (37 + ~11 evaluations)
!   5. Constants Pnt/Rnt/Pn1 outside search loop
!==============================================================================
subroutine failure_calc(failure_id,STRESS,UPSTRAN,e,Xt,Xc,Yt,Yc,Zt,Zc,S12,S13,S23, &
    EpsXt,EpsXc,EpsYt,EpsYc,EpsZt,EpsZc,GamS12,GamS13,GamS23, &
    THETAF,MGF,ANU12,ANU12f,E11,E11f,angles_out)

    implicit none
    real*8, dimension(6) :: e, STRESS, UPSTRAN, Xphi, angles_out
    real*8 :: Xt,Xc,Yt,Yc,Zt,Zc,S12,S13,S23,Phi,E11,THETAF,MGF,ANU12,ANU12f,E11f
    real*8 :: EpsXt,EpsXc,EpsYt,EpsYc,EpsZt,EpsZc,GamS12,GamS13,GamS23
    real*8 :: F1,F2,F3,F11,F22,F33,F44,F55,F66,F12,F13,F23
    real*8 :: THETA,Rn,Rn1,Rnt,SFP,TN1,TNT,Pnt,Pn1,THETAMAX,IFF,A
    real*8 :: THETA_A, THETA_B, THETA_C
    real*8, parameter :: Pi = 3.14159265358979d+0
    integer :: failure_id, IMAX, I, J, J_COARSE, J_LOW, J_HIGH

    select case (failure_id)

    case(7)
    !--- Puck: Fiber Fracture ---
    ! A: effective fiber stress accounting for biaxial misalignment (MGF)
    A = STRESS(1) - (ANU12 - ANU12f*MGF*(E11/E11f))*(STRESS(2)+STRESS(3))

    if (A .GE. 0.0d0) then
        e(1) = (1.0d0/Xt)*A          ! tensile FF index
    else
        e(1) = (-1.0d0/Xc)*A         ! compressive FF index (sign-fixed: always positive)
    end if

    !--- Puck: Inter-Fiber Fracture (IFF) ---
    ! Fixed: Pnt=-1/tan(2*THETAF), not -1/(2*tan(2*THETAF))
    ! These two formulas together ensure: max IFF = 1.0 at theta=+-THETAF under -Yc
    Rn  = Yt
    Rn1 = S12
    Pnt = -1.0d0 / tan(2.0d0*THETAF)
    Rnt = Yc / (2.0d0*tan(THETAF))
    Pn1 = Pnt*(Rn1/Rnt)

    e(2) = 0.0d0;  e(3) = 0.0d0;  e(4) = 0.0d0
    e(5) = 0.0d0;  e(6) = 0.0d0
    THETAMAX = 0.0d0
    THETA_A = 0.0d0;  THETA_B = 0.0d0;  THETA_C = 0.0d0

    ! === Pass 1: Coarse (5-degree) – 37 evaluations ===
    J_COARSE = 0
    do J = 0, 180, 5
        THETA = -Pi/2.0d0 + J*(Pi/180.0d0)
        SFP = STRESS(2)*cos(THETA)**2 + STRESS(3)*sin(THETA)**2 &
            + 2.0d0*STRESS(6)*sin(THETA)*cos(THETA)
        TNT = -STRESS(2)*sin(THETA)*cos(THETA) + STRESS(3)*sin(THETA)*cos(THETA) &
            + STRESS(6)*(cos(THETA)**2 - sin(THETA)**2)
        TN1 = STRESS(4)*cos(THETA) + STRESS(5)*sin(THETA)

        if (SFP .GE. 0.0d0) then
            IFF = SQRT((SFP/Rn)**2 + (TN1/(Rn1-Pn1*SFP))**2 + (TNT/(Rnt-Pnt*SFP))**2)
            if (IFF .GT. e(3)) then;  e(3) = IFF;  THETA_A = THETA;  end if
        else
            IFF = SQRT((TN1/(Rn1-Pn1*SFP))**2 + (TNT/(Rnt-Pnt*SFP))**2)
            if (J .EQ. 90) then
                if (IFF .GT. e(4)) then;  e(4) = IFF;  THETA_B = THETA;  end if
            else
                if (IFF .GT. e(5)) then;  e(5) = IFF;  THETA_C = THETA;  end if
            end if
        end if
        if (IFF .GT. e(2)) then;  e(2) = IFF;  THETAMAX = THETA;  J_COARSE = J;  end if
    end do

    ! === Pass 2: Fine (1-degree) in +-5 degree window – ~11 evaluations ===
    J_LOW  = max(0,   J_COARSE - 5)
    J_HIGH = min(180, J_COARSE + 5)
    do J = J_LOW, J_HIGH
        THETA = -Pi/2.0d0 + J*(Pi/180.0d0)
        SFP = STRESS(2)*cos(THETA)**2 + STRESS(3)*sin(THETA)**2 &
            + 2.0d0*STRESS(6)*sin(THETA)*cos(THETA)
        TNT = -STRESS(2)*sin(THETA)*cos(THETA) + STRESS(3)*sin(THETA)*cos(THETA) &
            + STRESS(6)*(cos(THETA)**2 - sin(THETA)**2)
        TN1 = STRESS(4)*cos(THETA) + STRESS(5)*sin(THETA)

        if (SFP .GE. 0.0d0) then
            IFF = SQRT((SFP/Rn)**2 + (TN1/(Rn1-Pn1*SFP))**2 + (TNT/(Rnt-Pnt*SFP))**2)
            if (IFF .GT. e(3)) then;  e(3) = IFF;  THETA_A = THETA;  end if
        else
            IFF = SQRT((TN1/(Rn1-Pn1*SFP))**2 + (TNT/(Rnt-Pnt*SFP))**2)
            if (J .EQ. 90) then
                if (IFF .GT. e(4)) then;  e(4) = IFF;  THETA_B = THETA;  end if
            else
                if (IFF .GT. e(5)) then;  e(5) = IFF;  THETA_C = THETA;  end if
            end if
        end if
        if (IFF .GT. e(2)) then;  e(2) = IFF;  THETAMAX = THETA;  end if
    end do

    ! Output fracture angles in degrees
    angles_out(2) = THETAMAX * 180.0d0 / Pi   ! overall max IFF angle
    angles_out(3) = THETA_A  * 180.0d0 / Pi   ! Mode A angle
    angles_out(4) = merge(THETA_B*180.0d0/Pi, -999.0d0, e(4).GT.0.0d0)
    angles_out(5) = merge(THETA_C*180.0d0/Pi, -999.0d0, e(5).GT.0.0d0)

    case default
        write(*,*) "ERROR: unsupported failure_id =", failure_id
    end select

end subroutine failure_calc


    subroutine damage_calc (failure_id,damage_id,STRESS,UPSTRAN,dmg,e,fflags,KINC,beta_ft,beta_fc,beta_mt,beta_mc,beta_s,&
    &EpsXt,EpsXc,EpsYt,EpsYc,EpsZt,EpsZc,GamS12,GamS13,GamS23,a_ft,a_fc,a_mt,a_mc,a_s,n_ft,n_fc,n_mt,n_mc,n_s,&
    &G_ft,G_fc,G_mt,G_mc,G_IIC,le,let,alpha,beta,E11,E22,E33,G12,G13,G23,Xt,Xc,Yt,Yc,Zt,Zc,S12,S13,S23)    
    
    implicit none
    real*8, dimension  (6) :: e,STRESS,UPSTRAN,dmg,d_index
    real*8 :: beta_ft,beta_fc,beta_mc,beta_mt,beta_s
    real*8 :: a_ft,a_fc,a_mt,a_mc,a_s,n_ft,n_fc,n_mt,n_mc,n_s
    real*8 :: EpsXt,EpsXc,EpsYt,EpsYc,EpsZt,EpsZc,GamS12,GamS13,GamS23
    real*8 :: E11,E22,E33,G12,G13,G23,Xt,Xc,Yt,Yc,Zt,Zc,S12,S13,S23    
    real*8 :: G_ft,G_fc,G_mt,G_mc,G_IIC,le,let,alpha,beta,E11D,E22D,E33D,G12D,G13D,G23D    
    integer:: damage_id,failure_id,I,KINC,fflags(6)
    d_index = 0.0d0
    select case (damage_id)

    case(1) 
!#----------------Instantaneous Degredation-----------------#

! Fiber Tensile/Compressive damage

    if((e(1).GT.1.0) .AND. (STRESS(1).GE.0.0)) then
     
        fflags(1)=KINC   
        dmg(1)=(0.99999999-beta_ft)

        if(failure_id.EQ.5) then ! Additional induced shear damage for Hashin

            dmg(4)=(0.99999999-beta_s)

        end if
                  
    else if((e(1).GT.1.0) .AND. (STRESS(1).LT.0.0)) then  

        fflags(1)=-KINC  
        dmg(1)=(0.99999999-beta_fc)  

    end if

! Matrix Tensile/Compressive damage

    if((e(2).GT.1.0) .AND. (STRESS(2).GE.0.0)) then  

        fflags(2)=KINC    
        dmg(2)= (0.99999999-beta_mt)  

        if(failure_id.EQ.5) then ! Additional induced shear damage for Hashin

            dmg(4)=(0.99999999-beta_s)
            dmg(6)=(0.99999999-beta_s)

        end if 
        
        if(failure_id.EQ.7) then ! Additional induced shear damage for Puck

            dmg(3)= (0.99999999-beta_mt)
            dmg(4)= (0.99999999-beta_s)
            dmg(6)= (0.99999999-beta_s)
            
        end if           
        
    else if((e(2).GT.1.0) .AND. (STRESS(2).LT.0.0)) then  

        fflags(2)=-KINC  
        dmg(2)= (0.99999999-beta_mc)  

        if(failure_id.EQ.5) then ! Additional induced shear damage for Hashin

            dmg(4)=(0.99999999-beta_s)
            dmg(6)=(0.99999999-beta_s)

        end if 

        if(failure_id.EQ.7) then ! Additional induced shear damage for Puck

            dmg(3)= (0.99999999-beta_mc)
            dmg(4)= (0.99999999-beta_s)
            dmg(6)= (0.99999999-beta_s)
            
        end if    

    end if

! Interlaminar Tensile/Compressive damage

    if((e(3).GT.1.0) .AND. (STRESS(3).GE.0.0)) then  

        fflags(3)=KINC    
        dmg(3)=(0.99999999-beta_mt) 

    else if((e(3).GT.1.0) .AND. (STRESS(2).LT.0.0)) then  

        fflags(3)=-KINC  
        dmg(3)=(0.99999999-beta_mc)  

    end if     

! Shear damage

    do I=4,6 

    if(e(I).GT.1.0) then 

        fflags(I)=KINC  
        dmg(I)= (0.99999999-beta_s)  

    end if
                      
    end do                  


    case(2)

!#----------------Recursive Degredation-----------------#

! Fiber Tensile/Compressive damage

    if(((e(1).GT.1.0).OR.(dmg(1).LT.1.0)) .AND. (STRESS(1).GE.0.0)) then 

        fflags(1)=KINC    
        dmg(1)= dmg(1)*(0.99999999d+0-beta_ft) 

        
        if(failure_id.EQ.5) then ! Additional induced shear damage for Hashin

            dmg(4)= dmg(4)*(0.99999999-beta_s)

        end if   

    else if(((e(1).GT.1.0).OR.(dmg(1).LT.1.0)) .AND. (STRESS(1).LT.0.0)) then  

        fflags(1)=-KINC 
        dmg(1)= dmg(1)*(0.99999999d+0-beta_fc)  

    end if

! Matrix Tensile/Compressive damage

    if(((e(2).GT.1.0).OR.(dmg(2).LT.1.0)) .AND. (STRESS(2).GE.0.0)) then  
    
        fflags(2)=KINC   
        dmg(2)=dmg(2)*(0.99999999d+0-beta_mt)

        if(failure_id.EQ.5) then ! Additional induced shear damage for Hashin

            dmg(4)=dmg(4)*(0.99999999-beta_s)
            dmg(6)=dmg(6)*(0.99999999-beta_s)

        end if 

        if(failure_id.EQ.7) then ! Additional induced shear damage for Puck

            dmg(3)= dmg(3)*(0.99999999-beta_mt)
            dmg(4)= dmg(4)*(0.99999999-beta_s)
            dmg(6)= dmg(6)*(0.99999999-beta_s)
            
        end if         

    else if(((e(2).GT.1.0).OR.(dmg(2).LT.1.0)) .AND. (STRESS(2).LT.0.0)) then 

        fflags(2)=-KINC  
        dmg(2)=dmg(2)*(0.99999999d+0-beta_mc)  
 
         if(failure_id.EQ.5) then ! Additional induced shear damage for Hashin

            dmg(4)=dmg(4)*(0.99999999-beta_s)
            dmg(6)=dmg(6)*(0.99999999-beta_s)

        end if  

        if(failure_id.EQ.7) then ! Additional induced shear damage for Puck

            dmg(3)= dmg(3)*(0.99999999-beta_mc)
            dmg(4)= dmg(4)*(0.99999999-beta_s)
            dmg(6)= dmg(6)*(0.99999999-beta_s)
            
        end if
 
    end if

! Interlaminar Tensile/Compressive damage

    if(((e(3).GT.1.0).OR.(dmg(3).LT.1.0)) .AND. (STRESS(3).GE.0.0)) then  

        fflags(3)=KINC    
        dmg(3)= dmg(3)*(0.99999999d+0-beta_mt)  

    else if(((e(3).GT.1.0).OR.(dmg(3).LT.1.0)) .AND. (STRESS(3).LT.0.0)) then  

        fflags(3)=-KINC  
        dmg(3)= dmg(3)*(0.99999999d+0-beta_mc)  

    end if     

! Shear damage

    do I=4,6 

    if((e(I).GT.1.0).OR.(dmg(I).LT.1.0)) THEN

        fflags(I)=KINC  
        dmg(I)=dmg(I)*(0.99999999d+0-beta_s)  

    end if
                    
    end do

    case(3) 
    
!#----------------Exponential Degredation-----------------#

! Fiber Tensile/Compressive damage

    if(((e(1).GT.1.0).OR.(dmg(1).LT.1.0)) .AND. (UPSTRAN(1).GE.0.0)) then  

        fflags(1)=KINC    
        dmg(1)= EXP(-a_ft*(UPSTRAN(1)- EpsXt)/(n_ft*EpsXt))

        if (failure_id.EQ.5) then ! Additional induced shear damage for Hashin

            dmg(4)= EXP(-a_ft*(UPSTRAN(1)- EpsXt)/(n_ft*EpsXt))

        end if

   
    else if(((e(1).GT.1.0).OR.(dmg(1).LT.1.0)) .AND. (UPSTRAN(1).LT.0.0)) then  

        fflags(1)=-KINC  
        dmg(1)= EXP(-a_fc*(abs(UPSTRAN(1))- EpsXc)/(n_fc*EpsXc))  


    end if

! Matrix Tensile/Compressive damage

    if(((e(2).GT.1.0).OR.(dmg(2).LT.1.0)) .AND. (UPSTRAN(2).GE.0.0)) then  

        fflags(2)=KINC    
        dmg(2)= EXP(-a_mt*(UPSTRAN(2)- EpsYt)/(n_mt*EpsYt))   

        if (failure_id.EQ.5) then ! Additional induced shear damage for Hashin

            dmg(4)= EXP(-a_mt*(UPSTRAN(2)- EpsYt)/(n_mt*EpsYt))
            dmg(6)= EXP(-a_mt*(UPSTRAN(2)- EpsYt)/(n_mt*EpsYt))       

        end if 
 
        if(failure_id.EQ.7) then ! Additional induced shear damage for Puck

            dmg(3)= EXP(-a_mt*(UPSTRAN(2)- EpsYt)/(n_mt*EpsYt))
            dmg(4)= EXP(-a_mt*(UPSTRAN(2)- EpsYt)/(n_mt*EpsYt))
            dmg(6)= EXP(-a_mt*(UPSTRAN(2)- EpsYt)/(n_mt*EpsYt))
            
        end if  
 
    else if(((e(2).GT.1.0).OR.(dmg(2).LT.1.0)) .AND. (UPSTRAN(2).LT.0.0)) then  

        fflags(2)=-KINC 
        dmg(2)= EXP(-a_mc*(abs(UPSTRAN(2))- EpsYc)/(n_mc*EpsYc))  
          
        if (failure_id.EQ.5) then ! Additional induced shear damage for Hashin

            dmg(4)= EXP(-a_mc*(abs(UPSTRAN(2))- EpsYc)/(n_mc*EpsYc))
            dmg(6)= EXP(-a_mc*(abs(UPSTRAN(2))- EpsYc)/(n_mc*EpsYc))       

        end if 

        if(failure_id.EQ.7) then ! Additional induced shear damage for Puck

            dmg(3)= EXP(-a_mc*(abs(UPSTRAN(2))- EpsYc)/(n_mc*EpsYc))
            dmg(4)= EXP(-a_mc*(abs(UPSTRAN(2))- EpsYc)/(n_mc*EpsYc))
            dmg(6)= EXP(-a_mc*(abs(UPSTRAN(2))- EpsYc)/(n_mc*EpsYc))
            
        end if  

    end if

! Interlaminar Tensile/Compressive damage

    if(((e(3).GT.1.0).OR.(dmg(3).LT.1.0)) .AND. (UPSTRAN(3).GE.0.0)) then  

        fflags(3)=KINC    
        dmg(3)= EXP(-a_mt*(UPSTRAN(3)- EpsZt)/(n_mt*EpsZt)) 

    else if(((e(3).GT.1.0).OR.(dmg(3).LT.1.0)) .AND. (UPSTRAN(3).LT.0.0)) then  

        fflags(3)=-KINC  
        dmg(3)= EXP(-a_mc*abs((UPSTRAN(3))- EpsZc)/(n_mc*EpsZc))  

    end if     

! Shear damage

    if((e(4).GT.1.0).OR.(dmg(4).LT.1.0)) then  

        fflags(4)=KINC  
        dmg(4)= EXP(-a_s*abs((UPSTRAN(4))- GamS12)/(n_s*GamS12))  

    end if                  

    if((e(5).GT.1.0).OR.(dmg(5).LT.1.0)) then  

        fflags(5)=KINC  
        dmg(5)= EXP(-a_s*abs((UPSTRAN(5))- GamS13)/(n_s*GamS13))  

    end if

    if((e(6).GT.1.0).OR.(dmg(6).LT.1.0)) then 

        fflags(6)=KINC  
        dmg(6)= EXP(-a_s*abs((UPSTRAN(6))- GamS23)/(n_s*GamS23))  

    end if                      
    
    case(4) 
    
!#--------------Constant Stress Degredation---------------#

! Fiber Tensile/Compressive damage

    if((e(1).GT.1.0) .AND. (STRESS(1).GE.0.0)) then
     
        fflags(1)=KINC   
        dmg(1)=dmg(1)*(1.0d+0/e(1)) 

        if (failure_id.EQ.5) then ! Additional induced shear damage for Hashin
        
            dmg(4)=dmg(4)*(1.0d+0/e(1))
            
        end if         


    else if((e(1).GT.1.0) .AND. (STRESS(1).LT.0.0)) then  

        fflags(1)=-KINC  
        dmg(1)=dmg(1)*(1.0d+0/e(1))   

    end if

! Matrix Tensile/Compressive damage

    if((e(2).GT.1.0) .AND. (STRESS(2).GE.0.0)) then  

        fflags(2)=KINC    
        dmg(2)=dmg(2)*(1.0d+0/e(2))   
 
        if(failure_id.EQ.5) then ! Additional induced shear damage for Hashin

            dmg(4)=dmg(4)*(1.0d+0/e(2))
            dmg(6)=dmg(4)*(1.0d+0/e(2))

        end if

        if(failure_id.EQ.7) then ! Additional induced shear damage for Puck

            dmg(3)= dmg(3)*(1.0d+0/e(2))
            dmg(4)= dmg(4)*(1.0d+0/e(2))
            dmg(6)= dmg(6)*(1.0d+0/e(2))
            
        end if 
 
 
    else if((e(2).GT.1.0) .AND. (STRESS(2).LT.0.0)) then  

        fflags(2)=-KINC  
        dmg(2)=dmg(2)*(1.0d+0/e(2)) 

        if(failure_id.EQ.5) then ! Additional induced shear damage for Hashin

            dmg(4)=dmg(4)*(1.0d+0/e(2))
            dmg(6)=dmg(4)*(1.0d+0/e(2))

        end if

        if(failure_id.EQ.7) then ! Additional induced shear damage for Puck

            dmg(3)= dmg(3)*(1.0d+0/e(2))
            dmg(4)= dmg(4)*(1.0d+0/e(2))
            dmg(6)= dmg(6)*(1.0d+0/e(2))
            
        end if 

    end if

! Interlaminar Tensile/Compressive damage

    if((e(3).GT.1.0) .AND. (STRESS(3).GE.0.0)) then  

        fflags(3)=KINC    
        dmg(3)=dmg(3)*(1.0d+0/e(3))  

    else if((e(3).GT.1.0) .AND. (STRESS(2).LT.0.0)) then  

        fflags(3)=-KINC  
        dmg(3)=dmg(3)*(1.0d+0/e(3))   

    end if     

! Shear damage

    do I=4,6 

    if(e(I).GT.1.0) then 

        fflags(I)=KINC  
        dmg(I)=dmg(I)*(1.0d+0/e(I)) 

    end if
                      
    end do      
            
    case(5) 
!#------Continium damage mechanics: Crack Band Theory-----#

! Fiber Tensile/Compressive damage

    if((e(1).GT.1.0) .AND. (UPSTRAN(1).GE.0.0)) then  

        fflags(1)=KINC   

        ! Calculate Degraded E11 modulus 

        E11D =((1.0/E11)+(UPSTRAN(1)- EpsXt)/(Xt*(1.0-(le*Xt*(UPSTRAN(1)- EpsXt))/(2*G_ft))))**(-1)
        d_index(1) = 1.0d+0 -(E11D/E11) 

    else if((e(1).GT.1.0) .AND. (UPSTRAN(1).LT.0.0)) then  

        fflags(1)=-KINC                   

        ! Calculate Degraded E11 modulus as per Eq (12)

        E11D =((1.0/E11)+(abs(UPSTRAN(1))-EpsXc))/(Xc*(1.0-(le*Xc*((abs(UPSTRAN(1)-EpsXc)))/(2*G_fc))))**(-1)
        d_index(1) = 1.0d+0 - (E11D/E11) 

    end if          

! Matrix Tensile/Compressive damage

    if((e(2).GT.1.0) .AND. (UPSTRAN(2).GE.0.0)) then  

        fflags(2)=KINC    

        ! Calculate Degraded E22,G12,G23 modulI                   

        E22D =((1.0/E22)+(UPSTRAN(2)-EpsYt)/(Yt*(1.0-(le*Yt*(UPSTRAN(2)-EpsYt))/(2*G_mt))))**(-1)

        G12D =((1.0/G12)+(abs(UPSTRAN(4))-GamS12))/&
        &(2*S12*(1.0-(le*S12*((abs(UPSTRAN(4)-GamS12)))/(4*G_IIC))))**(-1)

        G23D =((1.0/G23)+(abs(UPSTRAN(6))-GamS23))/&
        &(2*S23*(1.0-(let*S23*((abs(UPSTRAN(6)-GamS23)))/(4*G_IIC))))**(-1)                   

        d_index(2) = 1.0d+0 -(E22D/E22)
        d_index(4) = 1.0d+0 -(G12D/G12)                                      
        d_index(6) = 1.0d+0 -(G23D/G23)  


    else if((e(2).GT.1.0) .AND. (UPSTRAN(2).LT.0.0)) then

        fflags(2)=-KINC                      

        ! Calculate Degraded E22,G12,G23 modulI                    

        E22D =((1.0/E22)+(abs(UPSTRAN(2))-EpsYc)/(Yc*(1.0-(le*Yc*((abs(UPSTRAN(2))-EpsYc)))/(2*G_mc))))**(-1)

        G12D =((1.0/G12)+(abs(UPSTRAN(4))-GamS12)/&
        &(2*S12*(1.0-(le*S12*((abs(UPSTRAN(4))-GamS12)))/(4*G_IIC))))**(-1)

        G23D =((1.0/G23)+(abs(UPSTRAN(6))-GamS23)/&
        &(2*S23*(1.0-(let*S23*((abs(UPSTRAN(6))-GamS23)))/(4*G_IIC))))**(-1)                   

        d_index(2) = 1.0d+0 -(E22D/E22)
        d_index(4) = 1.0d+0 -(G12D/G12)                                        
        d_index(6) = 1.0d+0 -(G23D/G23)

    end if 

! Interlaminar Tensile/Compressive damage

    if((e(3).GT.1.0) .AND. (UPSTRAN(3).GE.0.0)) then 

        fflags(3)=KINC                

        ! Calculate Degraded (E33,G23,G13) moduli

        E33D =((1.0/E33)+(UPSTRAN(3)-EpsZt)/(Zt*(1.0-(let*Zt*(UPSTRAN(3)-EpsZt))/(2*G_mt))))**(-1)

        G23D =((1.0/G23)+((abs(UPSTRAN(6))-GamS23))/&
        &(2*S23*(1.0-(let*S23*((abs(UPSTRAN(6))-GamS23)))/(4*G_IIC))))**(-1)

        G13D =((1.0/G13)+((abs(UPSTRAN(5))-GamS13))/&
        &(2*S13*(1.0-(let*S13*((abs(UPSTRAN(5))-GamS13)))/(4*G_IIC))))**(-1)                    

        d_index(3) = 1.0d+0 -(E33D/E33)
        d_index(5) = 1.0d+0 -(G13D/G13)                                        
        d_index(6) = 1.0d+0 -(G23D/G23)                 

    else if((e(3).GT.1.0) .AND. (UPSTRAN(3).LT.0.0)) then  

        fflags(3)=-KINC    

        ! Calculate Degraded (E33,G23,G13) moduli 

        E33D =((1.0/E33)+(abs(UPSTRAN(3))-EpsZc)/(Zc*(1.0-(let*Zc*(UPSTRAN(3)-EpsZc))/(2*G_mc))))**(-1)

        G23D =((1.0/G23)+((abs(UPSTRAN(6))-GamS23))/&
        &(2*S23*(1.0-(let*S23*((abs(UPSTRAN(6))-GamS23)))/(4*G_IIC))))**(-1)

        G13D =((1.0/G13)+((abs(UPSTRAN(5))-GamS13))/&
        &(2*S13*(1.0-(let*S13*((abs(UPSTRAN(5))-GamS13)))/(4*G_IIC))))**(-1)                    

        d_index(3) = 1.0d+0 -(E33D/E33)
        d_index(5) = 1.0d+0 -(G13D/G13)                                        
        d_index(6) = 1.0d+0 -(G23D/G23) 

    end if 

! Assign Constiutive Matrix Damage Variables 

    dmg(1) = abs(0.999999-(d_index(1)))
    dmg(2) = abs(0.999999-(d_index(2)))
    dmg(3) = abs(0.999999-(d_index(3)))               
    dmg(4) = abs(0.999999-((1.0d+0-dmg(1))*(1.0d+0-alpha*dmg(2))*(1.0d+0-beta*d_index(4))))
    dmg(5) = abs(0.999999-((1.0d+0-dmg(1))*(1.0d+0-alpha*dmg(3))*(1.0d+0-beta*d_index(5))))
    dmg(6) = abs(0.999999-((1.0d+0-dmg(1))*(1.0d+0-alpha*dmg(2))*(1.0d+0-alpha*dmg(3))*(1.0d+0-beta*d_index(5))))     


    case default 

! Print error message and exit the program for invalid selection

    write(*,*) "ERROR: UNRECOGNISED DEGRADATION METHOD"
    stop
    
    end select 
    

    end subroutine damage_calc
