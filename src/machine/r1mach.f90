!** R1MACH
REAL FUNCTION R1MACH(I)
  IMPLICIT NONE
  !>
  !***
  !  Return floating point machine dependent constants.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  R1
  !***
  ! **Type:**      SINGLE PRECISION (R1MACH-S, D1MACH-D)
  !***
  ! **Keywords:**  MACHINE CONSTANTS
  !***
  ! **Author:**  Fox, P. A., (Bell Labs)
  !           Hall, A. D., (Bell Labs)
  !           Schryer, N. L., (Bell Labs)
  !***
  ! **Description:**
  !
  !   R1MACH can be used to obtain machine-dependent parameters for the
  !   local machine environment.  It is a function subprogram with one
  !   (input) argument, and can be referenced as follows:
  !
  !        A = R1MACH(I)
  !
  !   where I=1,...,5.  The (output) value of A above is determined by
  !   the (input) value of I.  The results for various values of I are
  !   discussed below.
  !
  !   R1MACH(1) = B**(EMIN-1), the smallest positive magnitude.
  !   R1MACH(2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
  !   R1MACH(3) = B**(-T), the smallest relative spacing.
  !   R1MACH(4) = B**(1-T), the largest relative spacing.
  !   R1MACH(5) = LOG10(B)
  !
  !   Assume single precision numbers are represented in the T-digit,
  !   base-B form
  !
  !              sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
  !
  !   where 0 .LE. X(I) .LT. B for I=1,...,T, 0 .LT. X(1), and
  !   EMIN .LE. E .LE. EMAX.
  !
  !   The values of B, T, EMIN and EMAX are provided in I1MACH as
  !   follows:
  !   I1MACH(10) = B, the base.
  !   I1MACH(11) = T, the number of base-B digits.
  !   I1MACH(12) = EMIN, the smallest exponent E.
  !   I1MACH(13) = EMAX, the largest exponent E.
  !
  !   To alter this function for a particular environment, the desired
  !   set of DATA statements should be activated by removing the C from
  !   column 1.  Also, the values of R1MACH(1) - R1MACH(4) should be
  !   checked for consistency with the local operating system.
  !
  !***
  ! **References:**  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for
  !                 a portable library, ACM Transactions on Mathematical
  !                 Software 4, 2 (June 1978), pp. 177-188.
  !***
  ! **Routines called:**  XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   790101  DATE WRITTEN
  !   890213  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900618  Added DEC RISC constants.  (WRB)
  !   900723  Added IBM RS 6000 constants.  (WRB)
  !   910710  Added HP 730 constants.  (SMR)
  !   911114  Added Convex IEEE constants.  (WRB)
  !   920121  Added SUN -r8 compiler option constants.  (WRB)
  !   920229  Added Touchstone Delta i860 constants.  (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !   920625  Added CONVEX -p8 and -pd8 compiler option constants.
  !           (BKS, WRB)
  !   930201  Added DEC Alpha and SGI constants.  (RWC and WRB)
  
  INTEGER I
  REAL, SAVE :: rmach(5)
  LOGICAL :: FIRST = .TRUE.
  !
  !* FIRST EXECUTABLE STATEMENT  D1MACH
  !IF ( I<1 .OR. I>5 ) CALL XERMSG('SLATEC','R1MACH','I OUT OF BOUNDS',1,2)

  IF(FIRST) THEN
    rmach = (/ TINY(1.), HUGE(1.), SPACING(0.5), SPACING(1.), &
      LOG10( REAL(RADIX(1.),4) ) /)
    FIRST = .FALSE.
  END IF

  R1MACH = rmach(I)

END FUNCTION R1MACH
