!DECK SNRM2
REAL FUNCTION SNRM2(N,Sx,Incx)
  IMPLICIT NONE
  INTEGER i, Incx, j, N, nn
  !***BEGIN PROLOGUE  SNRM2
  !***PURPOSE  Compute the Euclidean length (L2 norm) of a vector.
  !***LIBRARY   SLATEC (BLAS)
  !***CATEGORY  D1A3B
  !***TYPE      SINGLE PRECISION (SNRM2-S, DNRM2-D, SCNRM2-C)
  !***KEYWORDS  BLAS, EUCLIDEAN LENGTH, EUCLIDEAN NORM, L2,
  !             LINEAR ALGEBRA, UNITARY, VECTOR
  !***AUTHOR  Lawson, C. L., (JPL)
  !           Hanson, R. J., (SNLA)
  !           Kincaid, D. R., (U. of Texas)
  !           Krogh, F. T., (JPL)
  !***DESCRIPTION
  !
  !                B L A S  Subprogram
  !    Description of Parameters
  !
  !     --Input--
  !        N  number of elements in input vector(s)
  !       SX  single precision vector with N elements
  !     INCX  storage spacing between elements of SX
  !
  !     --Output--
  !    SNRM2  single precision result (zero if N .LE. 0)
  !
  !     Euclidean norm of the N-vector stored in SX with storage
  !     increment INCX .
  !     If N .LE. 0, return with result = 0.
  !     If N .GE. 1, then INCX must be .GE. 1
  !
  !     Four Phase Method using two built-in constants that are
  !     hopefully applicable to all machines.
  !         CUTLO = maximum of  SQRT(U/EPS)  over all known machines.
  !         CUTHI = minimum of  SQRT(V)      over all known machines.
  !     where
  !         EPS = smallest no. such that EPS + 1. .GT. 1.
  !         U   = smallest positive no.   (underflow limit)
  !         V   = largest  no.            (overflow  limit)
  !
  !     Brief Outline of Algorithm.
  !
  !     Phase 1 scans zero components.
  !     Move to phase 2 when a component is nonzero and .LE. CUTLO
  !     Move to phase 3 when a component is .GT. CUTLO
  !     Move to phase 4 when a component is .GE. CUTHI/M
  !     where M = N for X() real and M = 2*N for complex.
  !
  !     Values for CUTLO and CUTHI.
  !     From the environmental parameters listed in the IMSL converter
  !     document the limiting values are as follows:
  !     CUTLO, S.P.   U/EPS = 2**(-102) for  Honeywell.  Close seconds are
  !                   Univac and DEC at 2**(-103)
  !                   Thus CUTLO = 2**(-51) = 4.44089E-16
  !     CUTHI, S.P.   V = 2**127 for Univac, Honeywell, and DEC.
  !                   Thus CUTHI = 2**(63.5) = 1.30438E19
  !     CUTLO, D.P.   U/EPS = 2**(-67) for Honeywell and DEC.
  !                   Thus CUTLO = 2**(-33.5) = 8.23181D-11
  !     CUTHI, D.P.   same as S.P.  CUTHI = 1.30438D19
  !     DATA CUTLO, CUTHI /8.232D-11,  1.304D19/
  !     DATA CUTLO, CUTHI /4.441E-16,  1.304E19/
  !
  !***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
  !                 Krogh, Basic linear algebra subprograms for Fortran
  !                 usage, Algorithm No. 539, Transactions on Mathematical
  !                 Software 5, 3 (September 1979), pp. 308-323.
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   791001  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  SNRM2
  INTEGER next
  REAL Sx(*), cutlo, cuthi, hitest, sum, xmax, zero, one
  SAVE cutlo, cuthi, zero, one
  DATA zero, one/0.0E0, 1.0E0/
  !
  DATA cutlo, cuthi/4.441E-16, 1.304E19/
  !***FIRST EXECUTABLE STATEMENT  SNRM2
  IF ( N>0 ) THEN
    !
    next = 200
    sum = zero
    nn = N*Incx
    !
    !                                                 BEGIN MAIN LOOP
    !
    i = 1
  ELSE
    SNRM2 = zero
    RETURN
  ENDIF
  100 CONTINUE
  SELECT CASE(next)
    CASE(200)
      GOTO 200
    CASE(300)
      GOTO 300
    CASE(600)
      GOTO 600
    CASE(700)
      GOTO 700
  END SELECT
  200 CONTINUE
  IF ( ABS(Sx(i))>cutlo ) GOTO 800
  next = 300
  xmax = zero
  !
  !                        PHASE 1.  SUM IS ZERO
  !
  300 CONTINUE
  IF ( Sx(i)==zero ) GOTO 900
  IF ( ABS(Sx(i))>cutlo ) GOTO 800
  !
  !                                PREPARE FOR PHASE 2.
  !
  next = 600
  GOTO 500
  !
  !                                PREPARE FOR PHASE 4.
  !
  400  i = j
  next = 700
  sum = (sum/Sx(i))/Sx(i)
  500  xmax = ABS(Sx(i))
  !
  sum = sum + (Sx(i)/xmax)**2
  GOTO 900
  !
  !                   PHASE 2.  SUM IS SMALL.
  !                             SCALE TO AVOID DESTRUCTIVE UNDERFLOW.
  !
  600 CONTINUE
  IF ( ABS(Sx(i))>cutlo ) THEN
    !
    !                  PREPARE FOR PHASE 3.
    !
    sum = (sum*xmax)*xmax
    GOTO 800
  ENDIF
  !
  !                     COMMON CODE FOR PHASES 2 AND 4.
  !                     IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW.
  !
  700 CONTINUE
  IF ( ABS(Sx(i))<=xmax ) THEN
    sum = sum + (Sx(i)/xmax)**2
  ELSE
    sum = one + sum*(xmax/Sx(i))**2
    xmax = ABS(Sx(i))
  ENDIF
  GOTO 900
  !
  !     FOR REAL OR D.P. SET HITEST = CUTHI/N
  !     FOR COMPLEX      SET HITEST = CUTHI/(2*N)
  !
  800  hitest = cuthi/N
  !
  !                   PHASE 3.  SUM IS MID-RANGE.  NO SCALING.
  !
  DO j = i, nn, Incx
    IF ( ABS(Sx(j))>=hitest ) GOTO 400
    sum = sum + Sx(j)**2
  ENDDO
  SNRM2 = SQRT(sum)
  RETURN
  !
  900  i = i + Incx
  IF ( i<=nn ) GOTO 100
  !
  !              END OF MAIN LOOP.
  !
  !              COMPUTE SQUARE ROOT AND ADJUST FOR SCALING.
  !
  SNRM2 = xmax*SQRT(sum)
  RETURN
END FUNCTION SNRM2