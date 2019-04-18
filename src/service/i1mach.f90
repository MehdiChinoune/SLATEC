!** I1MACH
INTEGER FUNCTION I1MACH(I)
  !>
  !  Return integer machine dependent constants.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  R1
  !***
  ! **Type:**      INTEGER (I1MACH-I)
  !***
  ! **Keywords:**  MACHINE CONSTANTS
  !***
  ! **Author:**  Fox, P. A., (Bell Labs)
  !           Hall, A. D., (Bell Labs)
  !           Schryer, N. L., (Bell Labs)
  !***
  ! **Description:**
  !
  !   I1MACH can be used to obtain machine-dependent parameters for the
  !   local machine environment.  It is a function subprogram with one
  !   (input) argument and can be referenced as follows:
  !
  !        K = I1MACH(I)
  !
  !   where I=1,...,16.  The (output) value of K above is determined by
  !   the (input) value of I.  The results for various values of I are
  !   discussed below.
  !
  !   I/O unit numbers:
  !     I1MACH( 1) = the standard input unit.
  !     I1MACH( 2) = the standard output unit.
  !     I1MACH( 3) = the standard punch unit.
  !     I1MACH( 4) = the standard error message unit.
  !
  !   Words:
  !     I1MACH( 5) = the number of bits per integer storage unit.
  !     I1MACH( 6) = the number of characters per integer storage unit.
  !
  !   Integers:
  !     assume integers are represented in the S-digit, base-A form
  !
  !                sign ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
  !
  !                where 0 .LE. X(I) .LT. A for I=0,...,S-1.
  !     I1MACH( 7) = A, the base.
  !     I1MACH( 8) = S, the number of base-A digits.
  !     I1MACH( 9) = A**S - 1, the largest magnitude.
  !
  !   Floating-Point Numbers:
  !     Assume floating-point numbers are represented in the T-digit,
  !     base-B form
  !                sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
  !
  !                where 0 .LE. X(I) .LT. B for I=1,...,T,
  !                0 .LT. X(1), and EMIN .LE. E .LE. EMAX.
  !     I1MACH(10) = B, the base.
  !
  !   Single-Precision:
  !     I1MACH(11) = T, the number of base-B digits.
  !     I1MACH(12) = EMIN, the smallest exponent E.
  !     I1MACH(13) = EMAX, the largest exponent E.
  !
  !   Double-Precision:
  !     I1MACH(14) = T, the number of base-B digits.
  !     I1MACH(15) = EMIN, the smallest exponent E.
  !     I1MACH(16) = EMAX, the largest exponent E.
  !
  !   To alter this function for a particular environment, the desired
  !   set of DATA statements should be activated by removing the C from
  !   column 1.  Also, the values of I1MACH(1) - I1MACH(4) should be
  !   checked for consistency with the local operating system.
  !
  !***
  ! **References:**  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for
  !                 a portable library, ACM Transactions on Mathematical
  !                 Software 4, 2 (June 1978), pp. 177-188.
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   750101  DATE WRITTEN
  !   891012  Added VAX G-floating constants.  (WRB)
  !   891012  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900618  Added DEC RISC constants.  (WRB)
  !   900723  Added IBM RS 6000 constants.  (WRB)
  !   901009  Correct I1MACH(7) for IBM Mainframes. Should be 2 not 16.
  !           (RWC)
  !   910710  Added HP 730 constants.  (SMR)
  !   911114  Added Convex IEEE constants.  (WRB)
  !   920121  Added SUN -r8 compiler option constants.  (WRB)
  !   920229  Added Touchstone Delta i860 constants.  (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !   920625  Added Convex -p8 and -pd8 compiler option constants.
  !           (BKS, WRB)
  !   930201  Added DEC Alpha and SGI constants.  (RWC and WRB)
  !   930618  Corrected I1MACH(5) for Convex -p8 and -pd8 compiler
  !           options.  (DWL, RWC and WRB).
  use ISO_FORTRAN_ENV, ONLY : INPUT_UNIT, OUTPUT_UNIT, ERROR_UNIT, &
    NUMERIC_STORAGE_SIZE, CHARACTER_STORAGE_SIZE
  INTEGER I
  INTEGER, PARAMETER :: imach(16) = (/ INPUT_UNIT, OUTPUT_UNIT, OUTPUT_UNIT, &
    ERROR_UNIT, NUMERIC_STORAGE_SIZE, CHARACTER_STORAGE_SIZE, &
    RADIX(1), DIGITS(1), HUGE(1), RADIX(1.), &
    DIGITS(1.), MINEXPONENT(1.), MAXEXPONENT(1.), &
    DIGITS(1.D0), MINEXPONENT(1.D0), MAXEXPONENT(1.D0) /)

  !* FIRST EXECUTABLE STATEMENT  I1MACH
  !IF ( I<1 .OR. I>16 ) CALL XERMSG('SLATEC','I1MACH','I OUT OF BOUNDS',1,2)

  I1MACH = imach(I)

END FUNCTION I1MACH
