!** R9KNUS
ELEMENTAL SUBROUTINE R9KNUS(Xnu,X,Bknu,Bknu1,Iswtch)
  !> Compute Bessel functions EXP(X)*K_{XNU}(X) and EXP(X)*K_{XNU+1}(X)
  !  for 0.0 <= XNU < 1.0.
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
  ! Compute Bessel functions EXP(X) * K{XNU} (X)  and
  ! EXP(X) * K_{XNU+1} (X) for 0.0 <= XNU < 1.0 .
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
  USE service, ONLY : eps_2_sp, tiny_sp, huge_sp
  !
  INTEGER, INTENT(OUT) :: Iswtch
  REAL(SP), INTENT(IN) :: X, Xnu
  REAL(SP), INTENT(OUT) :: Bknu, Bknu1
  !
  INTEGER :: i, ii, inu, n, nterms
  REAL(SP) :: a(15), a0, alnz, alpha(15), an, b0, beta(15), bknu0, bknud, bn, c0, &
    expx, p1, p2, p3, qq, result, sqrtx, v, vlnz, x2n, x2tov, xi, xmu, z, ztov
  INTEGER, PARAMETER :: ntc0k = 8, ntznu1 = 7
  REAL(SP), PARAMETER :: xnusml = SQRT(eps_2_sp/8._SP), xsml = 0.1_SP*eps_2_sp, &
    alnsml = LOG(tiny_sp), alnbig = LOG(huge_sp), alneps = LOG(0.1_SP*eps_2_sp)
  REAL(SP), PARAMETER :: c0kcs(16) = [ .060183057242626108_SP,-.15364871433017286_SP, &
    -.011751176008210492_SP,-.000852487888919795_SP,-.000061329838767496_SP, &
    -.000004405228124551_SP,-.000000316312467283_SP,-.000000022710719382_SP, &
    -.000000001630564460_SP,-.000000000117069392_SP,-.000000000008405206_SP, &
    -.000000000000603466_SP,-.000000000000043326_SP,-.000000000000003110_SP, &
    -.000000000000000223_SP,-.000000000000000016_SP ]
  REAL(SP), PARAMETER :: znu1cs(12) = [ .20330675699419173_SP, .14007793341321977_SP, &
    .007916796961001613_SP, .000339801182532104_SP, .000011741975688989_SP, &
    .000000339357570612_SP, .000000008425941769_SP, .000000000183336677_SP, &
    .000000000003549698_SP, .000000000000061903_SP, .000000000000000981_SP, &
    .000000000000000014_SP ]
  REAL(SP), PARAMETER :: euler = 0.57721566490153286_SP
  REAL(SP), PARAMETER :: sqpi2 = 1.2533141373155003_SP
  REAL(SP), PARAMETER :: aln2 = 0.69314718055994531_SP
  !* FIRST EXECUTABLE STATEMENT  R9KNUS
  ! ntc0k = INITS(c0kcs,0.1_SP*eps_2_sp)
  ! ntznu1 = INITS(znu1cs,0.1_SP*eps_2_sp)
  !
  IF( Xnu<0._SP .OR. Xnu>=1._SP ) THEN
    ERROR STOP 'R9KNUS : XNU MUST 0 <= XNU < 1'
  ELSEIF( X<=0. ) THEN
    ERROR STOP 'R9KNUS : X MUST BE > 0'
  END IF
  !
  Iswtch = 0
  IF( X>2._SP ) THEN
    !
    ! X IS LARGE.  FIND K_{XNU} (X) AND K_{XNU+1} (X) WITH Y. L. LUKE-S
    ! RATIONAL EXPANSION.
    !
    sqrtx = SQRT(X)
    IF( X>1._SP/xsml ) THEN
      !
      Bknu = sqpi2/sqrtx
      Bknu1 = Bknu
      RETURN
    END IF
  ELSE
    !
    ! X IS SMALL.  COMPUTE K_{XNU} (X) AND THE DERIVATIVE OF K_{XNU} (X)
    ! THEN FIND K_{XNU+1} (X).  XNU IS REDUCED TO THE INTERVAL (-.5,+.5)
    ! THEN TO (0., .5), BECAUSE K OF NEGATIVE ORDER (-NU) = K OF POSITIVE
    ! ORDER (+NU).
    !
    v = Xnu
    IF( Xnu>0.5_SP ) v = 1._SP - Xnu
    !
    ! CAREFULLY FIND (X/2)**XNU AND Z**XNU WHERE Z = X*X/4.
    alnz = 2._SP*(LOG(X)-aln2)
    !
    IF( X<=Xnu .AND. -0.5_SP*Xnu*alnz-aln2-LOG(Xnu)>alnbig ) THEN
      ERROR STOP 'R9KNUS : X SO SMALL BESSEL K_{XNU} OVERFLOWS'
    END IF
    !
    vlnz = v*alnz
    x2tov = EXP(0.5_SP*vlnz)
    ztov = 0._SP
    IF( vlnz>alnsml ) ztov = x2tov**2
    !
    a0 = 0.5_SP*GAMMA(1._SP+v)
    b0 = 0.5_SP*GAMMA(1._SP-v)
    c0 = -euler
    IF( ztov>0.5_SP .AND. v>xnusml ) c0 = -0.75_SP + &
      CSEVL((8._SP*v)*v-1._SP,c0kcs(1:ntc0k))
    !
    IF( ztov<=0.5_SP ) alpha(1) = (a0-ztov*b0)/v
    IF( ztov>0.5_SP ) alpha(1) = c0 - alnz*(0.75_SP+CSEVL(vlnz/0.35_SP+1._SP,&
      znu1cs(1:ntznu1)))*b0
    beta(1) = -0.5_SP*(a0+ztov*b0)
    !
    z = 0._SP
    IF( X>xsml ) z = 0.25_SP*X*X
    nterms = INT( MAX(2._SP,11._SP+(8._SP*alnz-25.19_SP-alneps)/(4.28_SP-alnz)) )
    DO i = 2, nterms
      xi = i - 1
      a0 = a0/(xi*(xi-v))
      b0 = b0/(xi*(xi+v))
      alpha(i) = (alpha(i-1)+2._SP*xi*a0)/(xi*(xi+v))
      beta(i) = (xi-0.5_SP*v)*alpha(i) - ztov*b0
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
    IF( -0.5_SP*(Xnu+1._SP)*alnz-2._SP*aln2>alnbig ) Iswtch = 1
    IF( Iswtch==1 ) RETURN
    bknud = expx*bknud*2._SP/(x2tov*X)
    !
    IF( Xnu<=0.5_SP ) Bknu1 = v*Bknu/X - bknud
    IF( Xnu<=0.5_SP ) RETURN
    !
    bknu0 = Bknu
    Bknu = -v*Bknu/X - bknud
    Bknu1 = 2._SP*Xnu*Bknu/X + bknu0
    RETURN
  END IF
  an = -1.56_SP + 4._SP/X
  bn = -0.29_SP - 0.22_SP/X
  nterms = MIN(15, INT( MAX(3._SP,an+bn*alneps) ) )
  !
  DO inu = 1, 2
    xmu = 0.
    IF( inu==1 .AND. Xnu>xnusml ) xmu = (4._SP*Xnu)*Xnu
    IF( inu==2 ) xmu = 4._SP*(ABS(Xnu)+1._SP)**2
    !
    a(1) = 1._SP - xmu
    a(2) = 9._SP - xmu
    a(3) = 25._SP - xmu
    IF( a(2)==0. ) result = sqpi2*(16._SP*X+xmu+7._SP)/(16._SP*X*sqrtx)
    IF( a(2)/=0. ) THEN
      !
      alpha(1) = 1._SP
      alpha(2) = (16._SP*X+a(2))/a(2)
      alpha(3) = ((768._SP*X+48._SP*a(3))*X+a(2)*a(3))/(a(2)*a(3))
      !
      beta(1) = 1._SP
      beta(2) = (16._SP*X+(xmu+7._SP))/a(2)
      beta(3) = ((768._SP*X+48._SP*(xmu+23._SP))*X+((xmu+62._SP)*xmu+129._SP))/(a(2)*a(3))
      !
      IF( nterms>=4 ) THEN
        DO i = 4, nterms
          n = i - 1
          x2n = 2*n - 1
          !
          a(i) = (x2n+2._SP)**2 - xmu
          qq = 16._SP*x2n/a(i)
          p1 = -x2n*(12*n*n-20*n-a(1))/((x2n-2._SP)*a(i)) - qq*X
          p2 = (12*n*n-28*n+8-a(1))/a(i) - qq*X
          p3 = -x2n*a(i-3)/((x2n-2._SP)*a(i))
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
END SUBROUTINE R9KNUS