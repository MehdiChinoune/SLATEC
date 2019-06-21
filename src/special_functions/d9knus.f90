!** D9KNUS
SUBROUTINE D9KNUS(Xnu,X,Bknu,Bknu1,Iswtch)
  !> Compute Bessel functions EXP(X)*K-SUB-XNU(X) and EXP(X)*
  !            K-SUB-XNU+1(X) for 0.0 <= XNU < 1.0.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C10B3
  !***
  ! **Type:**      DOUBLE PRECISION (R9KNUS-S, D9KNUS-D)
  !***
  ! **Keywords:**  BESSEL FUNCTION, FNLIB, SPECIAL FUNCTIONS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! Compute Bessel functions EXP(X) * K-sub-XNU (X)  and
  ! EXP(X) * K-sub-XNU+1 (X) for 0.0 <= XNU < 1.0 .
  !
  ! Series for C0K        on the interval  0.          to  2.50000E-01
  !                                        with weighted error   2.16E-32
  !                                         log weighted error  31.67
  !                               significant figures required  30.86
  !                                    decimal places required  32.40
  !
  ! Series for ZNU1       on the interval -7.00000E-01 to  0.
  !                                        with weighted error   2.45E-33
  !                                         log weighted error  32.61
  !                               significant figures required  31.85
  !                                    decimal places required  33.26
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  D1MACH, DCSEVL, INITDS, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890911  Removed unnecessary intrinsics.  (WRB)
  !   890911  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900720  Routine changed from user-callable to subsidiary.  (WRB)
  !   900727  Added EXTERNAL statement.  (WRB)
  !   920618  Removed space from variable names.  (RWC, WRB)
  USE service, ONLY : XERMSG, D1MACH
  INTEGER :: Iswtch
  REAL(DP) :: Xnu, X, Bknu, Bknu1
  INTEGER :: i, ii, inu, n, nterms
  REAL(DP) :: alpha(32), beta(32), a(32), alnz, a0, bknud, bknu0, b0, c0, expx, &
    p1, p2, p3, qq, result, sqrtx, v, vlnz, xi, xmu, x2n, x2tov, z, ztov, an, bn
  INTEGER, SAVE :: ntc0k, ntznu1
  REAL(DP), PARAMETER :: eta = 0.1_DP*D1MACH(3), xnusml = SQRT(D1MACH(3)/8._DP), &
    xsml = 0.1_DP*D1MACH(3), alnsml = LOG(D1MACH(1)), alnbig = LOG(D1MACH(2)), &
    alneps = LOG(0.1_DP*D1MACH(3))
  REAL(DP), PARAMETER :: c0kcs(29) = [ +.60183057242626108387577445180329E-1_DP, &
    -.15364871433017286092959755943124E+0_DP, -.11751176008210492040068229226213E-1_DP, &
    -.85248788891979509827048401550987E-3_DP, -.61329838767496791874098176922111E-4_DP, &
    -.44052281245510444562679889548505E-5_DP, -.31631246728384488192915445892199E-6_DP, &
    -.22710719382899588330673771793396E-7_DP, -.16305644608077609552274620515360E-8_DP, &
    -.11706939299414776568756044043130E-9_DP, -.84052063786464437174546593413792E-11_DP, &
    -.60346670118979991487096050737198E-12_DP, -.43326960335681371952045997366903E-13_DP, &
    -.31107358030203546214634697772237E-14_DP, -.22334078226736982254486133409840E-15_DP, &
    -.16035146716864226300635791528610E-16_DP, -.11512717363666556196035697705305E-17_DP, &
    -.82657591746836959105169479089258E-19_DP, -.59345480806383948172333436695984E-20_DP, &
    -.42608138196467143926499613023976E-21_DP, -.30591266864812876299263698370542E-22_DP, &
    -.21963541426734575224975501815516E-23_DP, -.15769113261495836071105750684760E-24_DP, &
    -.11321713935950320948757731048056E-25_DP, -.81286248834598404082792349714433E-27_DP, &
    -.58360900893453226552829349315949E-28_DP, -.41901241623610922519452337780905E-29_DP, &
    -.30083737960206435069530504212862E-30_DP, -.21599152067808647728342168089832E-31_DP ]
  REAL(DP), PARAMETER :: znu1cs(20) = [ +.203306756994191729674444001216911E+0_DP, &
    +.140077933413219771062943670790563E+0_DP, +.791679696100161352840972241972320E-2_DP, &
    +.339801182532104045352930092205750E-3_DP, +.117419756889893366664507228352690E-4_DP, &
    +.339357570612261680333825865475121E-6_DP, +.842594176976219910194629891264803E-8_DP, &
    +.183336677024850089184748150900090E-9_DP, +.354969844704416310863007064469557E-11_DP, &
    +.619032496469887332205244342078407E-13_DP, +.981964535680439424960346115456527E-15_DP, &
    +.142851314396490474211473563005985E-16_DP, +.191894921887825298966162467488436E-18_DP, &
    +.239430979739498914162313140597128E-20_DP, +.278890246815347354835870465474995E-22_DP, &
    +.304606650633033442582845214092865E-24_DP, +.313173237042191815771564260932089E-26_DP, &
    +.304133098987854951645174908005034E-28_DP, +.279840384636833084343185097659733E-30_DP, &
    +.244637186274497596485238794922666E-32_DP ]
  REAL(DP), PARAMETER :: euler = 0.57721566490153286060651209008240_DP
  REAL(DP), PARAMETER :: sqpi2 = +1.2533141373155002512078826424055_DP
  REAL(DP), PARAMETER :: aln2 = 0.69314718055994530941723212145818_DP
  LOGICAL, SAVE :: first = .TRUE.
  !* FIRST EXECUTABLE STATEMENT  D9KNUS
  IF( first ) THEN
    ntc0k = INITDS(c0kcs,29,eta)
    ntznu1 = INITDS(znu1cs,20,eta)
    first = .FALSE.
  END IF
  !
  IF( Xnu<0._DP .OR. Xnu>=1._DP )&
    CALL XERMSG('D9KNUS','XNU MUST BE GE 0 AND LT 1',1,2)
  IF( X<=0. ) CALL XERMSG('D9KNUS','X MUST BE GT 0',2,2)
  !
  Iswtch = 0
  IF( X>2._DP ) THEN
    !
    ! X IS LARGE.  FIND K-SUB-XNU (X) AND K-SUB-XNU+1 (X) WITH Y. L. LUKE-S
    ! RATIONAL EXPANSION.
    !
    sqrtx = SQRT(X)
    IF( X>1._DP/xsml ) THEN
      !
      Bknu = sqpi2/sqrtx
      Bknu1 = Bknu
      RETURN
    END IF
  ELSE
    !
    ! X IS SMALL.  COMPUTE K-SUB-XNU (X) AND THE DERIVATIVE OF K-SUB-XNU (X)
    ! THEN FIND K-SUB-XNU+1 (X).  XNU IS REDUCED TO THE INTERVAL (-.5,+.5)
    ! THEN TO (0., .5), BECAUSE K OF NEGATIVE ORDER (-NU) = K OF POSITIVE
    ! ORDER (+NU).
    !
    v = Xnu
    IF( Xnu>0.5_DP ) v = 1._DP - Xnu
    !
    ! CAREFULLY FIND (X/2)**XNU AND Z**XNU WHERE Z = X*X/4.
    alnz = 2._DP*(LOG(X)-aln2)
    !
    IF( X<=Xnu ) THEN
      IF( -0.5_DP*Xnu*alnz-aln2-LOG(Xnu)>alnbig )&
        CALL XERMSG('D9KNUS',&
        'X SO SMALL BESSEL K-SUB-XNU OVERFLOWS',3,2)
    END IF
    !
    vlnz = v*alnz
    x2tov = EXP(0.5_DP*vlnz)
    ztov = 0._DP
    IF( vlnz>alnsml ) ztov = x2tov**2
    !
    a0 = 0.5_DP*GAMMA(1._DP+v)
    b0 = 0.5_DP*GAMMA(1._DP-v)
    c0 = -euler
    IF( ztov>0.5_DP .AND. v>xnusml ) c0 = -0.75_DP + DCSEVL((8._DP*v)*v-1._DP,c0kcs,ntc0k)
    !
    IF( ztov<=0.5_DP ) alpha(1) = (a0-ztov*b0)/v
    IF( ztov>0.5_DP ) alpha(1) = c0 - alnz*(0.75_DP+DCSEVL(vlnz/0.35_DP+1._DP,&
      znu1cs,ntznu1))*b0
    beta(1) = -0.5_DP*(a0+ztov*b0)
    !
    z = 0._DP
    IF( X>xsml ) z = 0.25_DP*X*X
    nterms = INT( MAX( 2._DP, 11._DP+(8._DP*alnz-25.19_DP-alneps)/(4.28_DP-alnz) ) )
    DO i = 2, nterms
      xi = i - 1
      a0 = a0/(xi*(xi-v))
      b0 = b0/(xi*(xi+v))
      alpha(i) = (alpha(i-1)+2._DP*xi*a0)/(xi*(xi+v))
      beta(i) = (xi-0.5_DP*v)*alpha(i) - ztov*b0
    END DO
    !
    Bknu = alpha(nterms)
    bknud = beta(nterms)
    DO ii = 2, nterms
      i = nterms + 1 - ii
      Bknu = alpha(i) + Bknu*z
      bknud = beta(i) + bknud*z
    END DO
    !
    expx = EXP(X)
    Bknu = expx*Bknu/x2tov
    !
    IF( -0.5_DP*(Xnu+1._DP)*alnz-2._DP*aln2>alnbig ) Iswtch = 1
    IF( Iswtch==1 ) RETURN
    bknud = expx*bknud*2._DP/(x2tov*X)
    !
    IF( Xnu<=0.5_DP ) Bknu1 = v*Bknu/X - bknud
    IF( Xnu<=0.5_DP ) RETURN
    !
    bknu0 = Bknu
    Bknu = -v*Bknu/X - bknud
    Bknu1 = 2._DP*Xnu*Bknu/X + bknu0
    RETURN
  END IF
  an = -0.60_DP - 1.02_DP/X
  bn = -0.27_DP - 0.53_DP/X
  nterms = MIN( 32, INT( MAX( 3._DP, an+bn*alneps ) ) )
  !
  DO inu = 1, 2
    xmu = 0._DP
    IF( inu==1 .AND. Xnu>xnusml ) xmu = (4._DP*Xnu)*Xnu
    IF( inu==2 ) xmu = 4._DP*(ABS(Xnu)+1._DP)**2
    !
    a(1) = 1._DP - xmu
    a(2) = 9._DP - xmu
    a(3) = 25._DP - xmu
    IF( a(2)==0._DP ) result = sqpi2*(16._DP*X+xmu+7._DP)/(16._DP*X*sqrtx)
    IF( a(2)/=0._DP ) THEN
      !
      alpha(1) = 1._DP
      alpha(2) = (16._DP*X+a(2))/a(2)
      alpha(3) = ((768._DP*X+48._DP*a(3))*X+a(2)*a(3))/(a(2)*a(3))
      !
      beta(1) = 1._DP
      beta(2) = (16._DP*X+(xmu+7._DP))/a(2)
      beta(3) = ((768._DP*X+48._DP*(xmu+23._DP))*X+((xmu+62._DP)*xmu+129._DP))&
        /(a(2)*a(3))
      !
      IF( nterms>=4 ) THEN
        DO i = 4, nterms
          n = i - 1
          x2n = 2*n - 1
          !
          a(i) = (x2n+2._DP)**2 - xmu
          qq = 16._DP*x2n/a(i)
          p1 = -x2n*((12*n*n-20*n)-a(1))/((x2n-2._DP)*a(i)) - qq*X
          p2 = ((12*n*n-28*n+8)-a(1))/a(i) - qq*X
          p3 = -x2n*a(i-3)/((x2n-2._DP)*a(i))
          !
          alpha(i) = -p1*alpha(i-1) - p2*alpha(i-2) - p3*alpha(i-3)
          beta(i) = -p1*beta(i-1) - p2*beta(i-2) - p3*beta(i-3)
        END DO
      END IF
      !
      result = sqpi2*beta(nterms)/(sqrtx*alpha(nterms))
    END IF
    !
    IF( inu==1 ) Bknu = result
    IF( inu==2 ) Bknu1 = result
  END DO
  RETURN
END SUBROUTINE D9KNUS
