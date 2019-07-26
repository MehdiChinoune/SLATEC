!** DCHFCM
INTEGER ELEMENTAL FUNCTION DCHFCM(D1,D2,Delta)
  !> Check a single cubic for monotonicity.
  !***
  ! **Library:**   SLATEC (PCHIP)
  !***
  ! **Type:**      DOUBLE PRECISION (CHFCM-S, DCHFCM-D)
  !***
  ! **Author:**  Fritsch, F. N., (LLNL)
  !***
  ! **Description:**
  !
  !- Usage:
  !
  !        DOUBLE PRECISION  D1, D2, DELTA
  !        INTEGER  ISMON, DCHFCM
  !
  !        ISMON = DCHFCM (D1, D2, DELTA)
  !
  !- Arguments:
  !
  !     D1,D2:IN  are the derivative values at the ends of an interval.
  !
  !     DELTA:IN  is the data slope over that interval.
  !
  !- Function Return Values:
  !     ISMON : indicates the monotonicity of the cubic segment:
  !             ISMON = -3  if function is probably decreasing;
  !             ISMON = -1  if function is strictly decreasing;
  !             ISMON =  0  if function is constant;
  !             ISMON =  1  if function is strictly increasing;
  !             ISMON =  2  if function is non-monotonic;
  !             ISMON =  3  if function is probably increasing.
  !           If ABS(ISMON)=3, the derivative values are too close to the
  !           boundary of the monotonicity region to declare monotonicity
  !           in the presence of roundoff error.
  !
  !- Description:
  !
  !          DCHFCM:  Cubic Hermite Function -- Check Monotonicity.
  !
  !    Called by  DPCHCM  to determine the monotonicity properties of the
  !    cubic with boundary derivative values D1,D2 and chord slope DELTA.
  !
  !- Cautions:
  !     This is essentially the same as old DCHFMC, except that a
  !     new output value, -3, was added February 1989.  (Formerly, -3
  !     and +3 were lumped together in the single value 3.)  Codes that
  !     flag nonmonotonicity by "IF(ISMON=2)" need not be changed.
  !     Codes that check via "IF(ISMON>=3)" should change the test to
  !     "IF(IABS(ISMON)>=3)".  Codes that declare monotonicity via
  !     "IF(ISMON<=1)" should change to "IF(IABS(ISMON)<=1)".
  !
  !   REFER TO  DPCHCM
  !
  !***
  ! **Routines called:**  D1MACH

  !* REVISION HISTORY  (YYMMDD)
  !   820518  DATE WRITTEN
  !   820805  Converted to SLATEC library version.
  !   831201  Changed from  ISIGN  to SIGN  to correct bug that
  !           produced wrong sign when -1 < DELTA < 0 .
  !   890206  Added SAVE statements.
  !   890209  Added sign to returned value ISMON=3 and corrected
  !           argument description accordingly.
  !   890306  Added caution about changed output.
  !   890407  Changed name from DCHFMC to DCHFCM, as requested at the
  !           March 1989 SLATEC CML meeting, and made a few other
  !           minor modifications necessitated by this change.
  !   890407  Converted to new SLATEC format.
  !   890407  Modified DESCRIPTION to LDOC format.
  !   891214  Moved SAVE statements.  (WRB)
  USE service, ONLY : eps_dp
  !
  !  Fortran intrinsics used:  DSIGN.
  !  Other routines used:  D1MACH.
  !
  ! ----------------------------------------------------------------------
  !
  !  Programming notes:
  !
  !     TEN is actually a tuning parameter, which determines the width of
  !     the fuzz around the elliptical boundary.
  !
  !     To produce a single precision version, simply:
  !        a. Change DCHFCM to CHFCM wherever it occurs,
  !        b. Change the double precision declarations to REAL(SP), and
  !        c. Change the constants ZERO, ONE, ... to single precision.
  !
  !  DECLARE ARGUMENTS.
  !
  REAL(DP), INTENT(IN) :: D1, D2, Delta
  !
  !  DECLARE LOCAL VARIABLES.
  !
  INTEGER :: ismon, itrue
  REAL(DP) :: a, b, eps, phi
  !
  !        MACHINE-DEPENDENT PARAMETER -- SHOULD BE ABOUT 10*UROUND.
  !* FIRST EXECUTABLE STATEMENT  DCHFCM
  eps = 10._DP*eps_dp
  !
  !  MAKE THE CHECK.
  !
  IF( Delta/=0._DP ) THEN
    !        DATA IS NOT CONSTANT -- PICK UP SIGN.
    itrue = INT( SIGN(1._DP,Delta) )
    a = D1/Delta
    b = D2/Delta
    IF( (a<0._DP) .OR. (b<0._DP) ) THEN
      ismon = 2
    ELSEIF( (a<=3._DP-eps) .AND. (b<=3._DP-eps) ) THEN
      !           INSIDE SQUARE (0,3)X(0,3)  IMPLIES   OK.
      ismon = itrue
    ELSEIF( (a>4._DP+eps) .AND. (b>4._DP+eps) ) THEN
      !           OUTSIDE SQUARE (0,4)X(0,4)  IMPLIES   NONMONOTONIC.
      ismon = 2
    ELSE
      !           MUST CHECK AGAINST BOUNDARY OF ELLIPSE.
      a = a - 2._DP
      b = b - 2._DP
      phi = ((a*a+b*b)+a*b) - 3._DP
      IF( phi<-eps ) THEN
        ismon = itrue
      ELSEIF( phi>eps ) THEN
        ismon = 2
      ELSE
        !              TO CLOSE TO BOUNDARY TO TELL,
        !                  IN THE PRESENCE OF ROUND-OFF ERRORS.
        ismon = 3*itrue
      END IF
    END IF
    !        CASE OF CONSTANT DATA.
  ELSEIF( (D1==0._DP) .AND. (D2==0._DP) ) THEN
    ismon = 0
  ELSE
    ismon = 2
  END IF
  !
  !  RETURN VALUE.
  !
  DCHFCM = ismon
  !------------- LAST LINE OF DCHFCM FOLLOWS -----------------------------
END FUNCTION DCHFCM