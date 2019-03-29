!** DNRM2
REAL(8) FUNCTION DNRM2(N,Dx,Incx)
  IMPLICIT NONE
  !>
  !***
  !  Compute the Euclidean length (L2 norm) of a vector.
  !***
  ! **Library:**   SLATEC (BLAS)
  !***
  ! **Category:**  D1A3B
  !***
  ! **Type:**      DOUBLE PRECISION (SNRM2-S, DNRM2-D, SCNRM2-C)
  !***
  ! **Keywords:**  BLAS, EUCLIDEAN LENGTH, EUCLIDEAN NORM, L2,
  !             LINEAR ALGEBRA, UNITARY, VECTOR
  !***
  ! **Author:**  Lawson, C. L., (JPL)
  !           Hanson, R. J., (SNLA)
  !           Kincaid, D. R., (U. of Texas)
  !           Krogh, F. T., (JPL)
  !***
  ! **Description:**
  !
  !                B L A S  Subprogram
  !    Description of parameters
  !
  !     --Input--
  !        N  number of elements in input vector(s)
  !       DX  double precision vector with N elements
  !     INCX  storage spacing between elements of DX
  !
  !     --Output--
  !    DNRM2  double precision result (zero if N .LE. 0)
  !
  !     Euclidean norm of the N-vector stored in DX with storage
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
  !     move to phase 2 when a component is nonzero and .LE. CUTLO
  !     move to phase 3 when a component is .GT. CUTLO
  !     move to phase 4 when a component is .GE. CUTHI/M
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
  !***
  ! **References:**  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
  !                 Krogh, Basic linear algebra subprograms for Fortran
  !                 usage, Algorithm No. 539, Transactions on Mathematical
  !                 Software 5, 3 (September 1979), pp. 308-323.
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   791001  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  INTEGER i, Incx, j, N, nn, next
  REAL(8) :: Dx(*), hitest, sum, xmax
  REAL(8), PARAMETER :: zero = 0.0D0, one = 1.0D0
  REAL(8), PARAMETER :: cutlo = 8.232D-11, cuthi = 1.304D19
  !* FIRST EXECUTABLE STATEMENT  DNRM2
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
    DNRM2 = zero
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
  IF ( ABS(Dx(i))>cutlo ) GOTO 800
  next = 300
  xmax = zero
  !
  !                        PHASE 1.  SUM IS ZERO
  !
  300 CONTINUE
  IF ( Dx(i)==zero ) GOTO 900
  IF ( ABS(Dx(i))>cutlo ) GOTO 800
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
  sum = (sum/Dx(i))/Dx(i)
  500  xmax = ABS(Dx(i))
  !
  sum = sum + (Dx(i)/xmax)**2
  GOTO 900
  !
  !                   PHASE 2.  SUM IS SMALL.
  !                             SCALE TO AVOID DESTRUCTIVE UNDERFLOW.
  !
  600 CONTINUE
  IF ( ABS(Dx(i))>cutlo ) THEN
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
  IF ( ABS(Dx(i))<=xmax ) THEN
    sum = sum + (Dx(i)/xmax)**2
  ELSE
    sum = one + sum*(xmax/Dx(i))**2
    xmax = ABS(Dx(i))
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
    IF ( ABS(Dx(j))>=hitest ) GOTO 400
    sum = sum + Dx(j)**2
  ENDDO
  DNRM2 = SQRT(sum)
  RETURN
  !
  900  i = i + Incx
  IF ( i<=nn ) GOTO 100
  !
  !              END OF MAIN LOOP.
  !
  !              COMPUTE SQUARE ROOT AND ADJUST FOR SCALING.
  !
  DNRM2 = xmax*SQRT(sum)
  RETURN
END FUNCTION DNRM2
