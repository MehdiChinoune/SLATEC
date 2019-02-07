!*==SPENC.f90  processed by SPAG 6.72Dc at 10:56 on  6 Feb 2019
!DECK SPENC
FUNCTION SPENC(X)
  IMPLICIT NONE
  !*--SPENC5
  !*** Start of declarations inserted by SPAG
  REAL aln, CSEVL, pi26, R1MACH, SPENC, spencs, X, xbig
  INTEGER INITS, nspenc
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  SPENC
  !***PURPOSE  Compute a form of Spence's integral due to K. Mitchell.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C5
  !***TYPE      SINGLE PRECISION (SPENC-S, DSPENC-D)
  !***KEYWORDS  FNLIB, SPECIAL FUNCTIONS, SPENCE'S INTEGRAL
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  ! Evaluate a form of Spence's function defined by
  !        integral from 0 to X of  -LOG(1-Y)/Y  DY.
  ! For ABS(X) .LE. 1, the uniformly convergent expansion
  !        SPENC = sum K=1,infinity  X**K / K**2     is valid.
  !
  ! Spence's function can be used to evaluate much more general integral
  ! forms.  For example,
  !        integral from 0 to Z of  LOG(A*X+B)/(C*X+D)  DX  =
  !             LOG(ABS(B-A*D/C))*LOG(ABS(A*(C*X+D)/(A*D-B*C)))/C
  !             - SPENC (A*(C*Z+D)/(A*D-B*C)) / C.
  !
  ! Ref -- K. Mitchell, Philosophical Magazine, 40, p. 351 (1949).
  !        Stegun and Abromowitz, AMS 55, p. 1004.
  !
  !
  ! Series for SPEN       on the interval  0.          to  5.00000D-01
  !                                        with weighted error   6.82E-17
  !                                         log weighted error  16.17
  !                               significant figures required  15.22
  !                                    decimal places required  16.81
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  CSEVL, INITS, R1MACH
  !***REVISION HISTORY  (YYMMDD)
  !   780201  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  SPENC
  DIMENSION spencs(19)
  LOGICAL first
  SAVE spencs, pi26, nspenc, xbig, first
  DATA spencs(1)/.1527365598892406E0/
  DATA spencs(2)/.08169658058051014E0/
  DATA spencs(3)/.00581415714077873E0/
  DATA spencs(4)/.00053716198145415E0/
  DATA spencs(5)/.00005724704675185E0/
  DATA spencs(6)/.00000667454612164E0/
  DATA spencs(7)/.00000082764673397E0/
  DATA spencs(8)/.00000010733156730E0/
  DATA spencs(9)/.00000001440077294E0/
  DATA spencs(10)/.00000000198444202E0/
  DATA spencs(11)/.00000000027940058E0/
  DATA spencs(12)/.00000000004003991E0/
  DATA spencs(13)/.00000000000582346E0/
  DATA spencs(14)/.00000000000085767E0/
  DATA spencs(15)/.00000000000012768E0/
  DATA spencs(16)/.00000000000001918E0/
  DATA spencs(17)/.00000000000000290E0/
  DATA spencs(18)/.00000000000000044E0/
  DATA spencs(19)/.00000000000000006E0/
  DATA pi26/1.644934066848226E0/
  DATA first/.TRUE./
  !***FIRST EXECUTABLE STATEMENT  SPENC
  IF ( first ) THEN
    nspenc = INITS(spencs,19,0.1*R1MACH(3))
    xbig = 1.0/R1MACH(3)
  ENDIF
  first = .FALSE.
  !
  IF ( X>2.0 ) THEN
    !
    ! X .GT. 2.0
    !
    SPENC = 2.0*pi26 - 0.5*LOG(X)**2
    IF ( X<xbig ) SPENC = SPENC - (1.0+CSEVL(4.0/X-1.0,spencs,nspenc))/X
    GOTO 99999
  ELSEIF ( X<=1.0 ) THEN
    IF ( X>0.5 ) THEN
      !
      ! 0.5 .LT. X .LE. 1.0
      !
      SPENC = pi26
      IF ( X/=1.0 ) SPENC = pi26 - LOG(X)*LOG(1.0-X) - (1.0-X)&
        *(1.0+CSEVL(4.0*(1.0-X)-1.0,spencs,nspenc))
      RETURN
    ELSEIF ( X>=0.0 ) THEN
      !
      ! 0.0 .LE. X .LE. 0.5
      !
      SPENC = X*(1.0+CSEVL(4.0*X-1.0,spencs,nspenc))
      RETURN
    ELSEIF ( X>(-1.) ) THEN
      !
      ! -1.0 .LT. X .LT. 0.0
      !
      SPENC = -0.5*LOG(1.0-X)&
        **2 - X*(1.0+CSEVL(4.0*X/(X-1.0)-1.0,spencs,nspenc))/(X-1.0)
      RETURN
    ELSE
      !
      ! HERE IF X .LE. -1.0
      !
      aln = LOG(1.0-X)
      SPENC = -pi26 - 0.5*aln*(2.0*LOG(-X)-aln)
      IF ( X>(-xbig) ) SPENC = SPENC + &
        (1.0+CSEVL(4.0/(1.0-X)-1.0,spencs,nspenc))&
        /(1.0-X)
      RETURN
    ENDIF
  ENDIF
  !
  ! 1.0 .LT. X .LE. 2.0
  !
  SPENC = pi26 - 0.5*LOG(X)*LOG((X-1.0)**2/X) + (X-1.)&
    *(1.0+CSEVL(4.0*(X-1.)/X-1.0,spencs,nspenc))/X
  RETURN
  !
  99999 CONTINUE
  END FUNCTION SPENC
