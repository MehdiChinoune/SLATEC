!** SPENC
REAL(SP) ELEMENTAL FUNCTION SPENC(X)
  !> Compute a form of Spence's integral due to K. Mitchell.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C5
  !***
  ! **Type:**      SINGLE PRECISION (SPENC-S, DSPENC-D)
  !***
  ! **Keywords:**  FNLIB, SPECIAL FUNCTIONS, SPENCE'S INTEGRAL
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! Evaluate a form of Spence's function defined by
  !        integral from 0 to X of  -LOG(1-Y)/Y  DY.
  ! For ABS(X) <= 1, the uniformly convergent expansion
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
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  CSEVL, INITS, R1MACH

  !* REVISION HISTORY  (YYMMDD)
  !   780201  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  USE service, ONLY : eps_2_sp
  !
  REAL(SP), INTENT(IN) :: X
  !
  REAL(SP) :: aln
  INTEGER, PARAMETER :: nspenc = 9
  REAL(SP), PARAMETER :: xbig = 1._SP/eps_2_sp
  REAL(SP), PARAMETER :: spencs(19) = [ .1527365598892406_SP, .08169658058051014_SP, &
    .00581415714077873_SP, .00053716198145415_SP, .00005724704675185_SP, &
    .00000667454612164_SP, .00000082764673397_SP, .00000010733156730_SP, &
    .00000001440077294_SP, .00000000198444202_SP, .00000000027940058_SP, &
    .00000000004003991_SP, .00000000000582346_SP, .00000000000085767_SP, &
    .00000000000012768_SP, .00000000000001918_SP, .00000000000000290_SP, &
    .00000000000000044_SP, .00000000000000006_SP ]
  REAL(SP), PARAMETER ::  pi26 = 1.644934066848226_SP
  !* FIRST EXECUTABLE STATEMENT  SPENC
  ! nspenc = INITS(spencs,0.1_SP*eps_2_sp)
  !
  IF( X>2._SP ) THEN
    ! X > 2.0
    SPENC = 2._SP*pi26 - 0.5_SP*LOG(X)**2
    IF( X<xbig ) SPENC = SPENC - (1._SP+CSEVL(4._SP/X-1._SP,spencs(1:nspenc)))/X
  ELSEIF( X>1._SP ) THEN
    ! 1.0 < X <= 2.0
    SPENC = pi26 - 0.5_SP*LOG(X)*LOG((X-1._SP)**2/X) + (X-1._SP)&
      *(1._SP+CSEVL(4._SP*(X-1._SP)/X-1._SP,spencs(1:nspenc)))/X
  ELSEIF( X>0.5_SP ) THEN
    ! 0.5 < X <= 1.0
    SPENC = pi26
    IF( X/=1._SP ) SPENC = pi26 - LOG(X)*LOG(1._SP-X) - (1._SP-X)&
      *(1._SP+CSEVL(4._SP*(1._SP-X)-1._SP,spencs(1:nspenc)))
  ELSEIF( X>=0._SP ) THEN
    ! 0.0 <= X <= 0.5
    SPENC = X*(1._SP+CSEVL(4._SP*X-1._SP,spencs(1:nspenc)))
  ELSEIF( X>(-1._SP) ) THEN
    ! -1.0 < X < 0.0
    SPENC = -0.5_SP*LOG(1._SP-X)**2 &
      - X*(1._SP+CSEVL(4._SP*X/(X-1._SP)-1._SP,spencs(1:nspenc)))/(X-1._SP)
  ELSE
    ! HERE IF X <= -1.0
    aln = LOG(1._SP-X)
    SPENC = -pi26 - 0.5_SP*aln*(2._SP*LOG(-X)-aln)
    IF( X>(-xbig) ) THEN
      SPENC = SPENC + (1._SP+CSEVL(4._SP/(1._SP-X)-1._SP,spencs(1:nspenc)))/(1._SP-X)
    END IF
  END IF

  RETURN
END FUNCTION SPENC