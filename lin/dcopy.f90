!*==DCOPY.f90  processed by SPAG 6.72Dc at 10:58 on  6 Feb 2019
!DECK DCOPY
      SUBROUTINE DCOPY(N,Dx,Incx,Dy,Incy)
      IMPLICIT NONE
!*--DCOPY5
!*** Start of declarations inserted by SPAG
      INTEGER i , Incx , Incy , ix , iy , m , mp1 , N , ns
!*** End of declarations inserted by SPAG
!***BEGIN PROLOGUE  DCOPY
!***PURPOSE  Copy a vector.
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1A5
!***TYPE      DOUBLE PRECISION (SCOPY-S, DCOPY-D, CCOPY-C, ICOPY-I)
!***KEYWORDS  BLAS, COPY, LINEAR ALGEBRA, VECTOR
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
!       DX  double precision vector with N elements
!     INCX  storage spacing between elements of DX
!       DY  double precision vector with N elements
!     INCY  storage spacing between elements of DY
!
!     --Output--
!       DY  copy of vector DX (unchanged if N .LE. 0)
!
!     Copy double precision DX to double precision DY.
!     For I = 0 to N-1, copy DX(LX+I*INCX) to DY(LY+I*INCY),
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
!***END PROLOGUE  DCOPY
      DOUBLE PRECISION Dx(*) , Dy(*)
!***FIRST EXECUTABLE STATEMENT  DCOPY
      IF ( N<=0 ) RETURN
      IF ( Incx==Incy ) THEN
        IF ( Incx<1 ) THEN
        ELSEIF ( Incx==1 ) THEN
!
!     Code for both increments equal to 1.
!
!     Clean-up loop so remaining vector length is a multiple of 7.
!
          m = MOD(N,7)
          IF ( m/=0 ) THEN
            DO i = 1 , m
              Dy(i) = Dx(i)
            ENDDO
            IF ( N<7 ) RETURN
          ENDIF
          GOTO 100
        ELSE
!
!     Code for equal, positive, non-unit increments.
!
          ns = N*Incx
          DO i = 1 , ns , Incx
            Dy(i) = Dx(i)
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
        Dy(iy) = Dx(ix)
        ix = ix + Incx
        iy = iy + Incy
      ENDDO
      RETURN
 100  mp1 = m + 1
      DO i = mp1 , N , 7
        Dy(i) = Dx(i)
        Dy(i+1) = Dx(i+1)
        Dy(i+2) = Dx(i+2)
        Dy(i+3) = Dx(i+3)
        Dy(i+4) = Dx(i+4)
        Dy(i+5) = Dx(i+5)
        Dy(i+6) = Dx(i+6)
      ENDDO
      RETURN
99999 END SUBROUTINE DCOPY
