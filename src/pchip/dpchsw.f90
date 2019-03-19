!** DPCHSW
SUBROUTINE DPCHSW(Dfmax,Iextrm,D1,D2,H,Slope,Ierr)
  IMPLICIT NONE
  !>
  !***
  !  Limits excursion from data for DPCHCS
  !***
  ! **Library:**   SLATEC (PCHIP)
  !***
  ! **Type:**      DOUBLE PRECISION (PCHSW-S, DPCHSW-D)
  !***
  ! **Author:**  Fritsch, F. N., (LLNL)
  !***
  ! **Description:**
  !
  !         DPCHSW:  DPCHCS Switch Excursion Limiter.
  !
  !     Called by  DPCHCS  to adjust D1 and D2 if necessary to insure that
  !     the extremum on this interval is not further than DFMAX from the
  !     extreme data value.
  !
  ! ----------------------------------------------------------------------
  !
  !  Calling sequence:
  !
  !        INTEGER  IEXTRM, IERR
  !        DOUBLE PRECISION  DFMAX, D1, D2, H, SLOPE
  !
  !        CALL  DPCHSW (DFMAX, IEXTRM, D1, D2, H, SLOPE, IERR)
  !
  !   Parameters:
  !
  !     DFMAX -- (input) maximum allowed difference between F(IEXTRM) and
  !           the cubic determined by derivative values D1,D2.  (assumes
  !           DFMAX.GT.0.)
  !
  !     IEXTRM -- (input) index of the extreme data value.  (assumes
  !           IEXTRM = 1 or 2 .  Any value .NE.1 is treated as 2.)
  !
  !     D1,D2 -- (input) derivative values at the ends of the interval.
  !           (Assumes D1*D2 .LE. 0.)
  !          (output) may be modified if necessary to meet the restriction
  !           imposed by DFMAX.
  !
  !     H -- (input) interval length.  (Assumes  H.GT.0.)
  !
  !     SLOPE -- (input) data slope on the interval.
  !
  !     IERR -- (output) error flag.  should be zero.
  !           If IERR=-1, assumption on D1 and D2 is not satisfied.
  !           If IERR=-2, quadratic equation locating extremum has
  !                       negative discriminant (should never occur).
  !
  !    -------
  !    WARNING:  This routine does no validity-checking of arguments.
  !    -------
  !
  !  Fortran intrinsics used:  ABS, SIGN, SQRT.
  !
  !***
  ! **See also:**  DPCHCS
  !***
  ! **Routines called:**  D1MACH, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   820218  DATE WRITTEN
  !   820805  Converted to SLATEC library version.
  !   870707  Corrected XERROR calls for d.p. name(s).
  !   870707  Replaced DATA statement for SMALL with a use of D1MACH.
  !   870813  Minor cosmetic changes.
  !   890206  Corrected XERROR calls.
  !   890411  1. Added SAVE statements (Vers. 3.2).
  !           2. Added DOUBLE PRECISION declaration for D1MACH.
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900328  Added TYPE section.  (WRB)
  !   910408  Updated AUTHOR and DATE WRITTEN sections in prologue.  (WRB)
  !   920526  Eliminated possible divide by zero problem.  (FNF)
  !   930503  Improved purpose.  (FNF)
  
  !
  !**End
  !
  !  DECLARE ARGUMENTS.
  !
  INTEGER Iextrm, Ierr
  REAL(8) :: Dfmax, D1, D2, H, Slope
  !
  !  DECLARE LOCAL VARIABLES.
  !
  REAL(8) :: cp, fact, hphi, lambda, nu, one, phi, radcal, &
    rho, sigma, small, that, third, three, two, zero
  SAVE zero, one, two, three, fact
  SAVE third
  REAL(8) :: D1MACH
  !
  DATA zero/0.D0/, one/1.D0/, two/2.D0/, three/3.D0/, fact/100.D0/
  !        THIRD SHOULD BE SLIGHTLY LESS THAN 1/3.
  DATA third/0.33333D0/
  !
  !  NOTATION AND GENERAL REMARKS.
  !
  !     RHO IS THE RATIO OF THE DATA SLOPE TO THE DERIVATIVE BEING TESTED.
  !     LAMBDA IS THE RATIO OF D2 TO D1.
  !     THAT = T-HAT(RHO) IS THE NORMALIZED LOCATION OF THE EXTREMUM.
  !     PHI IS THE NORMALIZED VALUE OF P(X)-F1 AT X = XHAT = X-HAT(RHO),
  !           WHERE  THAT = (XHAT - X1)/H .
  !        THAT IS, P(XHAT)-F1 = D*H*PHI,  WHERE D=D1 OR D2.
  !     SIMILARLY,  P(XHAT)-F2 = D*H*(PHI-RHO) .
  !
  !      SMALL SHOULD BE A FEW ORDERS OF MAGNITUDE GREATER THAN MACHEPS.
  !* FIRST EXECUTABLE STATEMENT  DPCHSW
  small = fact*D1MACH(4)
  !
  !  DO MAIN CALCULATION.
  !
  IF ( D1==zero ) THEN
    !
    !        SPECIAL CASE -- D1.EQ.ZERO .
    !
    !          IF D2 IS ALSO ZERO, THIS ROUTINE SHOULD NOT HAVE BEEN CALLED.
    IF ( D2==zero ) GOTO 200
    !
    rho = Slope/D2
    !          EXTREMUM IS OUTSIDE INTERVAL WHEN RHO .GE. 1/3 .
    IF ( rho<third ) THEN
      that = (two*(three*rho-one))/(three*(two*rho-one))
      phi = that**2*((three*rho-one)/three)
      !
      !          CONVERT TO DISTANCE FROM F2 IF IEXTRM.NE.1 .
      IF ( Iextrm/=1 ) phi = phi - rho
      !
      !          TEST FOR EXCEEDING LIMIT, AND ADJUST ACCORDINGLY.
      hphi = H*ABS(phi)
      !           AT THIS POINT, HPHI.GT.0, SO DIVIDE IS OK.
      IF ( hphi*ABS(D2)>Dfmax ) D2 = SIGN(Dfmax/hphi,D2)
    ENDIF
  ELSE
    !
    rho = Slope/D1
    lambda = -D2/D1
    IF ( D2==zero ) THEN
      !
      !           SPECIAL CASE -- D2.EQ.ZERO .
      !
      !             EXTREMUM IS OUTSIDE INTERVAL WHEN RHO .GE. 1/3 .
      IF ( rho>=third ) GOTO 100
      cp = two - three*rho
      nu = one - two*rho
      that = one/(three*nu)
    ELSE
      IF ( lambda<=zero ) GOTO 200
      !
      !           NORMAL CASE -- D1 AND D2 BOTH NONZERO, OPPOSITE SIGNS.
      !
      nu = one - lambda - two*rho
      sigma = one - rho
      cp = nu + sigma
      IF ( ABS(nu)>small ) THEN
        radcal = (nu-(two*rho+one))*nu + sigma**2
        IF ( radcal<zero ) THEN
          !
          !     NEGATIVE VALUE OF RADICAL (SHOULD NEVER OCCUR).
          Ierr = -2
          CALL XERMSG('SLATEC','DPCHSW','NEGATIVE RADICAL',Ierr,1)
          RETURN
        ELSE
          that = (cp-SQRT(radcal))/(three*nu)
        ENDIF
      ELSE
        that = one/(two*sigma)
      ENDIF
    ENDIF
    phi = that*((nu*that-cp)*that+one)
    !
    !          CONVERT TO DISTANCE FROM F2 IF IEXTRM.NE.1 .
    IF ( Iextrm/=1 ) phi = phi - rho
    !
    !          TEST FOR EXCEEDING LIMIT, AND ADJUST ACCORDINGLY.
    hphi = H*ABS(phi)
    IF ( hphi*ABS(D1)>Dfmax ) THEN
      !           AT THIS POINT, HPHI.GT.0, SO DIVIDE IS OK.
      D1 = SIGN(Dfmax/hphi,D1)
      D2 = -lambda*D1
    ENDIF
  ENDIF
  !
  !  NORMAL RETURN.
  !
  100  Ierr = 0
  RETURN
  !
  !  ERROR RETURNS.
  !
  !     D1 AND D2 BOTH ZERO, OR BOTH NONZERO AND SAME SIGN.
  200  Ierr = -1
  CALL XERMSG('SLATEC','DPCHSW','D1 AND/OR D2 INVALID',Ierr,1)
  RETURN
  !------------- LAST LINE OF DPCHSW FOLLOWS -----------------------------
  RETURN
END SUBROUTINE DPCHSW
