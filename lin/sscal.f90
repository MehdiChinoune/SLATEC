!*==SSCAL.f90  processed by SPAG 6.72Dc at 10:58 on  6 Feb 2019
!DECK SSCAL
      SUBROUTINE SSCAL(N,Sa,Sx,Incx)
      IMPLICIT NONE
!*--SSCAL5
!***BEGIN PROLOGUE  SSCAL
!***PURPOSE  Multiply a vector by a constant.
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1A6
!***TYPE      SINGLE PRECISION (SSCAL-S, DSCAL-D, CSCAL-C)
!***KEYWORDS  BLAS, LINEAR ALGEBRA, SCALE, VECTOR
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
!       SA  single precision scale factor
!       SX  single precision vector with N elements
!     INCX  storage spacing between elements of SX
!
!     --Output--
!       SX  single precision result (unchanged if N .LE. 0)
!
!     Replace single precision SX by single precision SA*SX.
!     For I = 0 to N-1, replace SX(IX+I*INCX) with  SA * SX(IX+I*INCX),
!     where IX = 1 if INCX .GE. 0, else IX = 1+(1-N)*INCX.
!
!***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
!                 Krogh, Basic linear algebra subprograms for Fortran
!                 usage, Algorithm No. 539, Transactions on Mathematical
!                 Software 5, 3 (September 1979), pp. 308-323.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   791001  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900821  Modified to correct problem with a negative increment.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  SSCAL
      REAL Sa , Sx(*)
      INTEGER i , Incx , ix , m , mp1 , N
!***FIRST EXECUTABLE STATEMENT  SSCAL
      IF ( N<=0 ) RETURN
      IF ( Incx==1 ) THEN
!
!     Code for increment equal to 1.
!
!     Clean-up loop so remaining vector length is a multiple of 5.
!
        m = MOD(N,5)
        IF ( m/=0 ) THEN
          DO i = 1 , m
            Sx(i) = Sa*Sx(i)
          ENDDO
          IF ( N<5 ) RETURN
        ENDIF
        mp1 = m + 1
        DO i = mp1 , N , 5
          Sx(i) = Sa*Sx(i)
          Sx(i+1) = Sa*Sx(i+1)
          Sx(i+2) = Sa*Sx(i+2)
          Sx(i+3) = Sa*Sx(i+3)
          Sx(i+4) = Sa*Sx(i+4)
        ENDDO
      ELSE
!
!     Code for increment not equal to 1.
!
        ix = 1
        IF ( Incx<0 ) ix = (-N+1)*Incx + 1
        DO i = 1 , N
          Sx(ix) = Sa*Sx(ix)
          ix = ix + Incx
        ENDDO
        RETURN
      ENDIF
      END SUBROUTINE SSCAL
