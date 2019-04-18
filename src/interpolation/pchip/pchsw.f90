!** PCHSW
SUBROUTINE PCHSW(Dfmax,Iextrm,D1,D2,H,Slope,Ierr)
  !>
  !***
  !  Limits excursion from data for PCHCS
  !***
  ! **Library:**   SLATEC (PCHIP)
  !***
  ! **Type:**      SINGLE PRECISION (PCHSW-S, DPCHSW-D)
  !***
  ! **Author:**  Fritsch, F. N., (LLNL)
  !***
  ! **Description:**
  !
  !         PCHSW:  PCHCS Switch Excursion Limiter.
  !
  !     Called by  PCHCS  to adjust D1 and D2 if necessary to insure that
  !     the extremum on this interval is not further than DFMAX from the
  !     extreme data value.
  !
  ! ----------------------------------------------------------------------
  !
  !  Calling sequence:
  !
  !        INTEGER  IEXTRM, IERR
  !        REAL  DFMAX, D1, D2, H, SLOPE
  !
  !        CALL  PCHSW (DFMAX, IEXTRM, D1, D2, H, SLOPE, IERR)
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
  ! **See also:**  PCHCS
  !***
  ! **Routines called:**  R1MACH, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   820218  DATE WRITTEN
  !   820805  Converted to SLATEC library version.
  !   870707  Replaced DATA statement for SMALL with a use of R1MACH.
  !   890411  1. Added SAVE statements (Vers. 3.2).
  !           2. Added REAL R1MACH for consistency with D.P. version.
  !   890411  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900328  Added TYPE section.  (WRB)
  !   910408  Updated AUTHOR and DATE WRITTEN sections in prologue.  (WRB)
  !   920526  Eliminated possible divide by zero problem.  (FNF)
  !   930503  Improved purpose.  (FNF)
  USE service, ONLY : XERMSG, R1MACH
  !
  !**End
  !
  !  DECLARE ARGUMENTS.
  !
  INTEGER Iextrm, Ierr
  REAL Dfmax, D1, D2, H, Slope
  !
  !  DECLARE LOCAL VARIABLES.
  !
  REAL cp, hphi, lambda, nu, phi, radcal, rho, sigma, small, that
  !
  REAL, PARAMETER :: zero = 0., one = 1., two = 2., three = 3., fact = 100.
  !        THIRD SHOULD BE SLIGHTLY LESS THAN 1/3.
  REAL, PARAMETER :: third = 0.33333
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
  !* FIRST EXECUTABLE STATEMENT  PCHSW
  small = fact*R1MACH(4)
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
    END IF
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
          CALL XERMSG('SLATEC','PCHSW','NEGATIVE RADICAL',Ierr,1)
          RETURN
        ELSE
          that = (cp-SQRT(radcal))/(three*nu)
        END IF
      ELSE
        that = one/(two*sigma)
      END IF
    END IF
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
    END IF
  END IF
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
  CALL XERMSG('SLATEC','PCHSW','D1 AND/OR D2 INVALID',Ierr,1)
  RETURN
  !------------- LAST LINE OF PCHSW FOLLOWS ------------------------------
  RETURN
END SUBROUTINE PCHSW
