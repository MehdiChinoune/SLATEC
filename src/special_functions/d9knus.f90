!** D9KNUS
SUBROUTINE D9KNUS(Xnu,X,Bknu,Bknu1,Iswtch)
  IMPLICIT NONE
  !>
  !***
  !  Compute Bessel functions EXP(X)*K-SUB-XNU(X) and EXP(X)*
  !            K-SUB-XNU+1(X) for 0.0 .LE. XNU .LT. 1.0.
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
  ! EXP(X) * K-sub-XNU+1 (X) for 0.0 .LE. XNU .LT. 1.0 .
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
  ! **Routines called:**  D1MACH, DCSEVL, DGAMMA, INITDS, XERMSG

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

  REAL an, bn, eta
  INTEGER i, ii, inu, Iswtch, n, nterms
  REAL(8) :: Xnu, X, Bknu, Bknu1, alpha(32), beta(32), a(32), alnz, a0, &
    bknud, bknu0, b0, c0, expx, p1, p2, p3, qq, result, sqrtx, v, &
    vlnz, xi, xmu, x2n, x2tov, z, ztov
  INTEGER, EXTERNAL :: INITDS
  REAL(8), EXTERNAL :: D1MACH, DCSEVL, DGAMMA
  INTEGER, SAVE :: ntc0k, ntznu1
  REAL(8), SAVE :: xnusml, xsml, alnsml, alnbig
  REAL, SAVE :: alneps
  REAL(8), PARAMETER :: c0kcs(29) = [ +.60183057242626108387577445180329D-1, &
    -.15364871433017286092959755943124D+0, -.11751176008210492040068229226213D-1, &
    -.85248788891979509827048401550987D-3, -.61329838767496791874098176922111D-4, &
    -.44052281245510444562679889548505D-5, -.31631246728384488192915445892199D-6, &
    -.22710719382899588330673771793396D-7, -.16305644608077609552274620515360D-8, &
    -.11706939299414776568756044043130D-9, -.84052063786464437174546593413792D-11, &
    -.60346670118979991487096050737198D-12, -.43326960335681371952045997366903D-13, &
    -.31107358030203546214634697772237D-14, -.22334078226736982254486133409840D-15, &
    -.16035146716864226300635791528610D-16, -.11512717363666556196035697705305D-17, &
    -.82657591746836959105169479089258D-19, -.59345480806383948172333436695984D-20, &
    -.42608138196467143926499613023976D-21, -.30591266864812876299263698370542D-22, &
    -.21963541426734575224975501815516D-23, -.15769113261495836071105750684760D-24, &
    -.11321713935950320948757731048056D-25, -.81286248834598404082792349714433D-27, &
    -.58360900893453226552829349315949D-28, -.41901241623610922519452337780905D-29, &
    -.30083737960206435069530504212862D-30, -.21599152067808647728342168089832D-31 ]
  REAL(8), PARAMETER :: znu1cs(20) = [ +.203306756994191729674444001216911D+0, &
    +.140077933413219771062943670790563D+0, +.791679696100161352840972241972320D-2, &
    +.339801182532104045352930092205750D-3, +.117419756889893366664507228352690D-4, &
    +.339357570612261680333825865475121D-6, +.842594176976219910194629891264803D-8, &
    +.183336677024850089184748150900090D-9, +.354969844704416310863007064469557D-11, &
    +.619032496469887332205244342078407D-13, +.981964535680439424960346115456527D-15, &
    +.142851314396490474211473563005985D-16, +.191894921887825298966162467488436D-18, &
    +.239430979739498914162313140597128D-20, +.278890246815347354835870465474995D-22, &
    +.304606650633033442582845214092865D-24, +.313173237042191815771564260932089D-26, &
    +.304133098987854951645174908005034D-28, +.279840384636833084343185097659733D-30, &
    +.244637186274497596485238794922666D-32 ]
  REAL(8), PARAMETER :: euler = 0.57721566490153286060651209008240D0
  REAL(8), PARAMETER :: sqpi2 = +1.2533141373155002512078826424055D0
  REAL(8), PARAMETER :: aln2 = 0.69314718055994530941723212145818D0
  LOGICAL :: first = .TRUE.
  !* FIRST EXECUTABLE STATEMENT  D9KNUS
  IF ( first ) THEN
    eta = REAL(  0.1D0*D1MACH(3), 4 )
    ntc0k = INITDS(c0kcs,29,eta)
    ntznu1 = INITDS(znu1cs,20,eta)
    !
    xnusml = SQRT(D1MACH(3)/8.D0)
    xsml = 0.1D0*D1MACH(3)
    alnsml = LOG(D1MACH(1))
    alnbig = LOG(D1MACH(2))
    alneps = REAL(  LOG(0.1D0*D1MACH(3)), 4 )
    first = .FALSE.
  ENDIF
  !
  IF ( Xnu<0.D0.OR.Xnu>=1.D0 )&
    CALL XERMSG('SLATEC','D9KNUS','XNU MUST BE GE 0 AND LT 1',1,2)
  IF ( X<=0. ) CALL XERMSG('SLATEC','D9KNUS','X MUST BE GT 0',2,2)
  !
  Iswtch = 0
  IF ( X>2.0D0 ) THEN
    !
    ! X IS LARGE.  FIND K-SUB-XNU (X) AND K-SUB-XNU+1 (X) WITH Y. L. LUKE-S
    ! RATIONAL EXPANSION.
    !
    sqrtx = SQRT(X)
    IF ( X>1.0D0/xsml ) THEN
      !
      Bknu = sqpi2/sqrtx
      Bknu1 = Bknu
      RETURN
    ENDIF
  ELSE
    !
    ! X IS SMALL.  COMPUTE K-SUB-XNU (X) AND THE DERIVATIVE OF K-SUB-XNU (X)
    ! THEN FIND K-SUB-XNU+1 (X).  XNU IS REDUCED TO THE INTERVAL (-.5,+.5)
    ! THEN TO (0., .5), BECAUSE K OF NEGATIVE ORDER (-NU) = K OF POSITIVE
    ! ORDER (+NU).
    !
    v = Xnu
    IF ( Xnu>0.5D0 ) v = 1.0D0 - Xnu
    !
    ! CAREFULLY FIND (X/2)**XNU AND Z**XNU WHERE Z = X*X/4.
    alnz = 2.D0*(LOG(X)-aln2)
    !
    IF ( X<=Xnu ) THEN
      IF ( -0.5D0*Xnu*alnz-aln2-LOG(Xnu)>alnbig )&
        CALL XERMSG('SLATEC','D9KNUS',&
        'X SO SMALL BESSEL K-SUB-XNU OVERFLOWS',3,2)
    ENDIF
    !
    vlnz = v*alnz
    x2tov = EXP(0.5D0*vlnz)
    ztov = 0.0D0
    IF ( vlnz>alnsml ) ztov = x2tov**2
    !
    a0 = 0.5D0*DGAMMA(1.0D0+v)
    b0 = 0.5D0*DGAMMA(1.0D0-v)
    c0 = -euler
    IF ( ztov>0.5D0.AND.v>xnusml ) c0 = -0.75D0 + DCSEVL((8.0D0*v)*v-1.0D0,c0kcs,ntc0k)
    !
    IF ( ztov<=0.5D0 ) alpha(1) = (a0-ztov*b0)/v
    IF ( ztov>0.5D0 ) alpha(1) = c0 - alnz*(0.75D0+DCSEVL(vlnz/0.35D0+1.0D0,&
      znu1cs,ntznu1))*b0
    beta(1) = -0.5D0*(a0+ztov*b0)
    !
    z = 0.0D0
    IF ( X>xsml ) z = 0.25D0*X*X
    nterms = INT( MAX(2.0,11.0+(8.*REAL(alnz)-25.19-alneps)/(4.28-REAL(alnz))) )
    DO i = 2, nterms
      xi = i - 1
      a0 = a0/(xi*(xi-v))
      b0 = b0/(xi*(xi+v))
      alpha(i) = (alpha(i-1)+2.0D0*xi*a0)/(xi*(xi+v))
      beta(i) = (xi-0.5D0*v)*alpha(i) - ztov*b0
    ENDDO
    !
    Bknu = alpha(nterms)
    bknud = beta(nterms)
    DO ii = 2, nterms
      i = nterms + 1 - ii
      Bknu = alpha(i) + Bknu*z
      bknud = beta(i) + bknud*z
    ENDDO
    !
    expx = EXP(X)
    Bknu = expx*Bknu/x2tov
    !
    IF ( -0.5D0*(Xnu+1.D0)*alnz-2.0D0*aln2>alnbig ) Iswtch = 1
    IF ( Iswtch==1 ) RETURN
    bknud = expx*bknud*2.0D0/(x2tov*X)
    !
    IF ( Xnu<=0.5D0 ) Bknu1 = v*Bknu/X - bknud
    IF ( Xnu<=0.5D0 ) RETURN
    !
    bknu0 = Bknu
    Bknu = -v*Bknu/X - bknud
    Bknu1 = 2.0D0*Xnu*Bknu/X + bknu0
    RETURN
  ENDIF
  an = -0.60 - 1.02/REAL(X)
  bn = -0.27 - 0.53/REAL(X)
  nterms = MIN(32,MAX1(3.0,an+bn*alneps))
  !
  DO inu = 1, 2
    xmu = 0.D0
    IF ( inu==1.AND.Xnu>xnusml ) xmu = (4.0D0*Xnu)*Xnu
    IF ( inu==2 ) xmu = 4.0D0*(ABS(Xnu)+1.D0)**2
    !
    a(1) = 1.0D0 - xmu
    a(2) = 9.0D0 - xmu
    a(3) = 25.0D0 - xmu
    IF ( a(2)==0.D0 ) result = sqpi2*(16.D0*X+xmu+7.D0)/(16.D0*X*sqrtx)
    IF ( a(2)/=0.D0 ) THEN
      !
      alpha(1) = 1.0D0
      alpha(2) = (16.D0*X+a(2))/a(2)
      alpha(3) = ((768.D0*X+48.D0*a(3))*X+a(2)*a(3))/(a(2)*a(3))
      !
      beta(1) = 1.0D0
      beta(2) = (16.D0*X+(xmu+7.D0))/a(2)
      beta(3) = ((768.D0*X+48.D0*(xmu+23.D0))*X+((xmu+62.D0)*xmu+129.D0))&
        /(a(2)*a(3))
      !
      IF ( nterms>=4 ) THEN
        DO i = 4, nterms
          n = i - 1
          x2n = 2*n - 1
          !
          a(i) = (x2n+2.D0)**2 - xmu
          qq = 16.D0*x2n/a(i)
          p1 = -x2n*((12*n*n-20*n)-a(1))/((x2n-2.D0)*a(i)) - qq*X
          p2 = ((12*n*n-28*n+8)-a(1))/a(i) - qq*X
          p3 = -x2n*a(i-3)/((x2n-2.D0)*a(i))
          !
          alpha(i) = -p1*alpha(i-1) - p2*alpha(i-2) - p3*alpha(i-3)
          beta(i) = -p1*beta(i-1) - p2*beta(i-2) - p3*beta(i-3)
        ENDDO
      ENDIF
      !
      result = sqpi2*beta(nterms)/(sqrtx*alpha(nterms))
    ENDIF
    !
    IF ( inu==1 ) Bknu = result
    IF ( inu==2 ) Bknu1 = result
  ENDDO
  RETURN
END SUBROUTINE D9KNUS
