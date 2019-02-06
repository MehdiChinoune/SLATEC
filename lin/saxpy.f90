!*==SAXPY.f90  processed by SPAG 6.72Dc at 10:58 on  6 Feb 2019
!DECK SAXPY
      SUBROUTINE SAXPY(N,Sa,Sx,Incx,Sy,Incy)
      IMPLICIT NONE
!*--SAXPY5
!*** Start of declarations inserted by SPAG
      INTEGER i , Incx , Incy , ix , iy , m , mp1 , N , ns
!*** End of declarations inserted by SPAG
!***BEGIN PROLOGUE  SAXPY
!***PURPOSE  Compute a constant times a vector plus a vector.
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1A7
!***TYPE      SINGLE PRECISION (SAXPY-S, DAXPY-D, CAXPY-C)
!***KEYWORDS  BLAS, LINEAR ALGEBRA, TRIAD, VECTOR
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
!       SA  single precision scalar multiplier
!       SX  single precision vector with N elements
!     INCX  storage spacing between elements of SX
!       SY  single precision vector with N elements
!     INCY  storage spacing between elements of SY
!
!     --Output--
!       SY  single precision result (unchanged if N .LE. 0)
!
!     Overwrite single precision SY with single precision SA*SX +SY.
!     For I = 0 to N-1, replace  SY(LY+I*INCY) with SA*SX(LX+I*INCX) +
!       SY(LY+I*INCY),
!     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
!     defined in a similar way using INCY.
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
!   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  SAXPY
      REAL Sx(*) , Sy(*) , Sa
!***FIRST EXECUTABLE STATEMENT  SAXPY
      IF ( N<=0.OR.Sa==0.0E0 ) RETURN
      IF ( Incx==Incy ) THEN
        IF ( Incx<1 ) THEN
        ELSEIF ( Incx==1 ) THEN
!
!     Code for both increments equal to 1.
!
!     Clean-up loop so remaining vector length is a multiple of 4.
!
          m = MOD(N,4)
          IF ( m/=0 ) THEN
            DO i = 1 , m
              Sy(i) = Sy(i) + Sa*Sx(i)
            ENDDO
            IF ( N<4 ) RETURN
          ENDIF
          GOTO 100
        ELSE
!
!     Code for equal, positive, non-unit increments.
!
          ns = N*Incx
          DO i = 1 , ns , Incx
            Sy(i) = Sa*Sx(i) + Sy(i)
          ENDDO
          GOTO 99999
        ENDIF
      ENDIF
!
!     Code for unequal or nonpositive increments.
!
      ix = 1
      iy = 1
      IF ( Incx<0 ) ix = (-N+1)*Incx + 1
      IF ( Incy<0 ) iy = (-N+1)*Incy + 1
      DO i = 1 , N
        Sy(iy) = Sy(iy) + Sa*Sx(ix)
        ix = ix + Incx
        iy = iy + Incy
      ENDDO
      RETURN
 100  mp1 = m + 1
      DO i = mp1 , N , 4
        Sy(i) = Sy(i) + Sa*Sx(i)
        Sy(i+1) = Sy(i+1) + Sa*Sx(i+1)
        Sy(i+2) = Sy(i+2) + Sa*Sx(i+2)
        Sy(i+3) = Sy(i+3) + Sa*Sx(i+3)
      ENDDO
      RETURN
99999 END SUBROUTINE SAXPY
