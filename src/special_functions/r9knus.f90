!** R9KNUS
SUBROUTINE R9KNUS(Xnu,X,Bknu,Bknu1,Iswtch)
  !>
  !  Compute Bessel functions EXP(X)*K-SUB-XNU(X) and EXP(X)*
  !            K-SUB-XNU+1(X) for 0.0 .LE. XNU .LT. 1.0.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C10B3
  !***
  ! **Type:**      SINGLE PRECISION (R9KNUS-S, D9KNUS-D)
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
  ! Series for C0K        on the interval  0.          to  2.50000D-01
  !                                        with weighted error   1.60E-17
  !                                         log weighted error  16.79
  !                               significant figures required  15.99
  !                                    decimal places required  17.40
  !
  ! Series for ZNU1       on the interval -7.00000D-01 to  0.
  !                                        with weighted error   1.43E-17
  !                                         log weighted error  16.85
  !                               significant figures required  16.08
  !                                    decimal places required  17.38
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  CSEVL, INITS, R1MACH, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900720  Routine changed from user-callable to subsidiary.  (WRB)
  !   900727  Added EXTERNAL statement.  (WRB)
  !   920618  Removed space from variable names.  (RWC, WRB)
  USE service, ONLY : XERMSG, R1MACH
  INTEGER :: Iswtch
  REAL(SP) :: Bknu, Bknu1, X, Xnu
  INTEGER :: i, ii, inu, n, nterms
  REAL(SP) :: a(15), a0, alnz, alpha(15), an, b0, beta(15), bknu0, bknud, bn, c0, &
    expx, p1, p2, p3, qq, result, sqrtx, v, vlnz, x2n, x2tov, xi, xmu, z, ztov
  INTEGER, SAVE :: ntc0k, ntznu1
  REAL(SP), PARAMETER :: xnusml = SQRT(R1MACH(3)/8.0), xsml = 0.1*R1MACH(3), &
    alnsml = LOG(R1MACH(1)), alnbig = LOG(R1MACH(2)), alneps = LOG(0.1*R1MACH(3))
  REAL(SP), PARAMETER :: c0kcs(16) = [ .060183057242626108E0,-.15364871433017286E0, &
    -.011751176008210492E0,-.000852487888919795E0,-.000061329838767496E0, &
    -.000004405228124551E0,-.000000316312467283E0,-.000000022710719382E0, &
    -.000000001630564460E0,-.000000000117069392E0,-.000000000008405206E0, &
    -.000000000000603466E0,-.000000000000043326E0,-.000000000000003110E0, &
    -.000000000000000223E0,-.000000000000000016E0 ]
  REAL(SP), PARAMETER :: znu1cs(12) = [ .20330675699419173E0, .14007793341321977E0, &
    .007916796961001613E0, .000339801182532104E0, .000011741975688989E0, &
    .000000339357570612E0, .000000008425941769E0, .000000000183336677E0, &
    .000000000003549698E0, .000000000000061903E0, .000000000000000981E0, &
    .000000000000000014E0 ]
  REAL(SP), PARAMETER :: euler = 0.57721566490153286E0
  REAL(SP), PARAMETER :: sqpi2 = 1.2533141373155003E0
  REAL(SP), PARAMETER :: aln2 = 0.69314718055994531E0
  LOGICAL, SAVE :: first = .TRUE.
  !* FIRST EXECUTABLE STATEMENT  R9KNUS
  IF ( first ) THEN
    ntc0k = INITS(c0kcs,16,0.1*R1MACH(3))
    ntznu1 = INITS(znu1cs,12,0.1*R1MACH(3))
    first = .FALSE.
  END IF
  !
  IF ( Xnu<0..OR.Xnu>=1.0 )&
    CALL XERMSG('R9KNUS','XNU MUST BE GE 0 AND LT 1',1,2)
  IF ( X<=0. ) CALL XERMSG('R9KNUS','X MUST BE GT 0',2,2)
  !
  Iswtch = 0
  IF ( X>2.0 ) THEN
    !
    ! X IS LARGE.  FIND K-SUB-XNU (X) AND K-SUB-XNU+1 (X) WITH Y. L. LUKE-S
    ! RATIONAL EXPANSION.
    !
    sqrtx = SQRT(X)
    IF ( X>1.0/xsml ) THEN
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
    IF ( Xnu>0.5 ) v = 1.0 - Xnu
    !
    ! CAREFULLY FIND (X/2)**XNU AND Z**XNU WHERE Z = X*X/4.
    alnz = 2.0*(LOG(X)-aln2)
    !
    IF ( X<=Xnu ) THEN
      IF ( -0.5*Xnu*alnz-aln2-LOG(Xnu)>alnbig )&
        CALL XERMSG('R9KNUS',&
        'X SO SMALL BESSEL K-SUB-XNU OVERFLOWS',3,2)
    END IF
    !
    vlnz = v*alnz
    x2tov = EXP(0.5*vlnz)
    ztov = 0.0
    IF ( vlnz>alnsml ) ztov = x2tov**2
    !
    a0 = 0.5*GAMMA(1.0+v)
    b0 = 0.5*GAMMA(1.0-v)
    c0 = -euler
    IF ( ztov>0.5.AND.v>xnusml ) c0 = -0.75 + &
      CSEVL((8.0*v)*v-1.,c0kcs,ntc0k)
    !
    IF ( ztov<=0.5 ) alpha(1) = (a0-ztov*b0)/v
    IF ( ztov>0.5 ) alpha(1) = c0 - alnz*(0.75+CSEVL(vlnz/0.35+1.0,znu1cs,&
      ntznu1))*b0
    beta(1) = -0.5*(a0+ztov*b0)
    !
    z = 0.0
    IF ( X>xsml ) z = 0.25*X*X
    nterms = INT( MAX(2.0,11.0+(8.*alnz-25.19-alneps)/(4.28-alnz)) )
    DO i = 2, nterms
      xi = i - 1
      a0 = a0/(xi*(xi-v))
      b0 = b0/(xi*(xi+v))
      alpha(i) = (alpha(i-1)+2.0*xi*a0)/(xi*(xi+v))
      beta(i) = (xi-0.5*v)*alpha(i) - ztov*b0
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
    IF ( -0.5*(Xnu+1.)*alnz-2.0*aln2>alnbig ) Iswtch = 1
    IF ( Iswtch==1 ) RETURN
    bknud = expx*bknud*2.0/(x2tov*X)
    !
    IF ( Xnu<=0.5 ) Bknu1 = v*Bknu/X - bknud
    IF ( Xnu<=0.5 ) RETURN
    !
    bknu0 = Bknu
    Bknu = -v*Bknu/X - bknud
    Bknu1 = 2.0*Xnu*Bknu/X + bknu0
    RETURN
  END IF
  an = -1.56 + 4.0/X
  bn = -0.29 - 0.22/X
  nterms = MIN(15,MAX1(3.0,an+bn*alneps))
  !
  DO inu = 1, 2
    xmu = 0.
    IF ( inu==1.AND.Xnu>xnusml ) xmu = (4.0*Xnu)*Xnu
    IF ( inu==2 ) xmu = 4.0*(ABS(Xnu)+1.)**2
    !
    a(1) = 1.0 - xmu
    a(2) = 9.0 - xmu
    a(3) = 25.0 - xmu
    IF ( a(2)==0. ) result = sqpi2*(16.*X+xmu+7.)/(16.*X*sqrtx)
    IF ( a(2)/=0. ) THEN
      !
      alpha(1) = 1.0
      alpha(2) = (16.*X+a(2))/a(2)
      alpha(3) = ((768.*X+48.*a(3))*X+a(2)*a(3))/(a(2)*a(3))
      !
      beta(1) = 1.0
      beta(2) = (16.*X+(xmu+7.))/a(2)
      beta(3) = ((768.*X+48.*(xmu+23.))*X+((xmu+62.)*xmu+129.))/(a(2)*a(3))
      !
      IF ( nterms>=4 ) THEN
        DO i = 4, nterms
          n = i - 1
          x2n = 2*n - 1
          !
          a(i) = (x2n+2.)**2 - xmu
          qq = 16.*x2n/a(i)
          p1 = -x2n*(12*n*n-20*n-a(1))/((x2n-2.)*a(i)) - qq*X
          p2 = (12*n*n-28*n+8-a(1))/a(i) - qq*X
          p3 = -x2n*a(i-3)/((x2n-2.)*a(i))
          !
          alpha(i) = -p1*alpha(i-1) - p2*alpha(i-2) - p3*alpha(i-3)
          beta(i) = -p1*beta(i-1) - p2*beta(i-2) - p3*beta(i-3)
        END DO
      END IF
      !
      result = sqpi2*beta(nterms)/(sqrtx*alpha(nterms))
    END IF
    !
    IF ( inu==1 ) Bknu = result
    IF ( inu==2 ) Bknu1 = result
  END DO
  RETURN
END SUBROUTINE R9KNUS
