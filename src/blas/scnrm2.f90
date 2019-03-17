!DECK SCNRM2
REAL FUNCTION SCNRM2(N,Cx,Incx)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  SCNRM2
  !***PURPOSE  Compute the unitary norm of a complex vector.
  !***LIBRARY   SLATEC (BLAS)
  !***CATEGORY  D1A3B
  !***TYPE      COMPLEX (SNRM2-S, DNRM2-D, SCNRM2-C)
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
  !       CX  complex vector with N elements
  !     INCX  storage spacing between elements of CX
  !
  !     --Output--
  !   SCNRM2  single precision result (zero if N .LE. 0)
  !
  !     Unitary norm of the complex N-vector stored in CX with storage
  !     increment INCX.
  !     If N .LE. 0, return with result = 0.
  !     If N .GE. 1, then INCX must be .GE. 1
  !
  !     Four phase method using two built-in constants that are
  !     hopefully applicable to all machines.
  !         CUTLO = maximum of  SQRT(U/EPS)  over all known machines.
  !         CUTHI = minimum of  SQRT(V)      over all known machines.
  !     where
  !         EPS = smallest no. such that EPS + 1. .GT. 1.
  !         U   = smallest positive no.   (underflow limit)
  !         V   = largest  no.            (overflow  limit)
  !
  !     Brief outline of algorithm.
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
  !***END PROLOGUE  SCNRM2
  INTEGER i, Incx, N, nn
  LOGICAL imag, scale
  INTEGER next
  REAL cutlo, cuthi, hitest, sum, xmax, absx, zero, one
  COMPLEX Cx(*)
  SAVE cutlo, cuthi, zero, one
  DATA zero, one/0.0E0, 1.0E0/
  !
  DATA cutlo, cuthi/4.441E-16, 1.304E19/
  !***FIRST EXECUTABLE STATEMENT  SCNRM2
  IF ( N>0 ) THEN
    !
    next = 20
    sum = zero
    nn = N*Incx
    !
    !                                                 BEGIN MAIN LOOP
    !
    DO i = 1, nn, Incx
      absx = ABS(REAL(Cx(i)))
      imag = .FALSE.
      SELECT CASE(next)
        CASE(20)
          GOTO 20
        CASE(40)
          GOTO 40
        CASE(80)
          GOTO 80
        CASE(100)
          GOTO 100
        CASE(140)
          GOTO 140
      END SELECT
      20       IF ( absx>cutlo ) GOTO 120
      next = 40
      scale = .FALSE.
      !
      !                        PHASE 1.  SUM IS ZERO
      !
      40       IF ( absx==zero ) GOTO 160
      IF ( absx>cutlo ) GOTO 120
      !
      !                                PREPARE FOR PHASE 2.
      !
      next = 80
      60       scale = .TRUE.
      xmax = absx
      !
      sum = sum + (absx/xmax)**2
      GOTO 160
      !
      !                   PHASE 2.  SUM IS SMALL.
      !                             SCALE TO AVOID DESTRUCTIVE UNDERFLOW.
      !
      80       IF ( absx>cutlo ) THEN
      !
      !                  PREPARE FOR PHASE 3.
      !
      sum = (sum*xmax)*xmax
      GOTO 120
  ENDIF
  !
  !                     COMMON CODE FOR PHASES 2 AND 4.
  !                     IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW.
  !
  100      IF ( absx<=xmax ) THEN
  sum = sum + (absx/xmax)**2
ELSE
  sum = one + sum*(xmax/absx)**2
  xmax = absx
ENDIF
GOTO 160
!
120      next = 140
scale = .FALSE.
!
!     FOR REAL OR D.P. SET HITEST = CUTHI/N
!     FOR COMPLEX      SET HITEST = CUTHI/(2*N)
!
hitest = cuthi/N
!
!                   PHASE 3.  SUM IS MID-RANGE.  NO SCALING.
!
140      IF ( absx>=hitest ) THEN
!
!                                PREPARE FOR PHASE 4.
!
next =100
sum = (sum/absx)/absx
GOTO 60
ELSE
sum = sum + absx**2
ENDIF
!
!                  CONTROL SELECTION OF REAL AND IMAGINARY PARTS.
!
160      IF ( .NOT.(imag) ) THEN
absx = ABS(AIMAG(Cx(i)))
imag = .TRUE.
SELECT CASE(next)
CASE(20)
  GOTO 20
CASE(40)
  GOTO 40
CASE(80)
  GOTO 80
CASE(100)
  GOTO 100
CASE(140)
  GOTO 140
END SELECT
ENDIF
!
ENDDO
!
!              END OF MAIN LOOP.
!              COMPUTE SQUARE ROOT AND ADJUST FOR SCALING.
!
SCNRM2 = SQRT(sum)
IF ( scale ) SCNRM2 = SCNRM2*xmax
ELSE
SCNRM2 = zero
ENDIF
END FUNCTION SCNRM2
