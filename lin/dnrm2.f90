!*==DNRM2.f90  processed by SPAG 6.72Dc at 10:58 on  6 Feb 2019
!DECK DNRM2
      DOUBLE PRECISION FUNCTION DNRM2(N,Dx,Incx)
      IMPLICIT NONE
!*--DNRM25
!*** Start of declarations inserted by SPAG
      INTEGER i , Incx , j , N , nn
!*** End of declarations inserted by SPAG
!***BEGIN PROLOGUE  DNRM2
!***PURPOSE  Compute the Euclidean length (L2 norm) of a vector.
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1A3B
!***TYPE      DOUBLE PRECISION (SNRM2-S, DNRM2-D, SCNRM2-C)
!***KEYWORDS  BLAS, EUCLIDEAN LENGTH, EUCLIDEAN NORM, L2,
!             LINEAR ALGEBRA, UNITARY, VECTOR
!***AUTHOR  Lawson, C. L., (JPL)
!           Hanson, R. J., (SNLA)
!           Kincaid, D. R., (U. of Texas)
!           Krogh, F. T., (JPL)
!***DESCRIPTION
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
!***END PROLOGUE  DNRM2
      INTEGER next
      DOUBLE PRECISION Dx(*) , cutlo , cuthi , hitest , sum , xmax , zero , one
      SAVE cutlo , cuthi , zero , one
      DATA zero , one/0.0D0 , 1.0D0/
!
      DATA cutlo , cuthi/8.232D-11 , 1.304D19/
!***FIRST EXECUTABLE STATEMENT  DNRM2
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
        GOTO 99999
      ENDIF
 100  SELECT CASE(next)
        CASE(200)
          GOTO 200
        CASE(300)
          GOTO 300
        CASE(600)
          GOTO 600
        CASE(700)
          GOTO 700
      END SELECT

 200  IF ( ABS(Dx(i))>cutlo ) GOTO 800
      next = 300
      xmax = zero
!
!                        PHASE 1.  SUM IS ZERO
!
 300  IF ( Dx(i)==zero ) GOTO 900
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
 600  IF ( ABS(Dx(i))>cutlo ) THEN
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
 700  IF ( ABS(Dx(i))<=xmax ) THEN
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
      DO j = i , nn , Incx
        IF ( ABS(Dx(j))>=hitest ) GOTO 400
        sum = sum + Dx(j)**2
      ENDDO
      DNRM2 = SQRT(sum)
      GOTO 99999
!
 900  i = i + Incx
      IF ( i<=nn ) GOTO 100
!
!              END OF MAIN LOOP.
!
!              COMPUTE SQUARE ROOT AND ADJUST FOR SCALING.
!
      DNRM2 = xmax*SQRT(sum)
99999 END FUNCTION DNRM2
