!*==DSDOT.f90  processed by SPAG 6.72Dc at 10:58 on  6 Feb 2019
!DECK DSDOT
      DOUBLE PRECISION FUNCTION DSDOT(N,Sx,Incx,Sy,Incy)
      IMPLICIT NONE
!*--DSDOT5
!*** Start of declarations inserted by SPAG
      INTEGER i , Incx , Incy , kx , ky , N , ns
!*** End of declarations inserted by SPAG
!***BEGIN PROLOGUE  DSDOT
!***PURPOSE  Compute the inner product of two vectors with extended
!            precision accumulation and result.
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1A4
!***TYPE      DOUBLE PRECISION (DSDOT-D, DCDOT-C)
!***KEYWORDS  BLAS, COMPLEX VECTORS, DOT PRODUCT, INNER PRODUCT,
!             LINEAR ALGEBRA, VECTOR
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
!       SY  single precision vector with N elements
!     INCY  storage spacing between elements of SY
!
!     --Output--
!    DSDOT  double precision dot product (zero if N.LE.0)
!
!     Returns D.P. dot product accumulated in D.P., for S.P. SX and SY
!     DSDOT = sum for I = 0 to N-1 of  SX(LX+I*INCX) * SY(LY+I*INCY),
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
!***END PROLOGUE  DSDOT
      REAL Sx(*) , Sy(*)
!***FIRST EXECUTABLE STATEMENT  DSDOT
      DSDOT = 0.0D0
      IF ( N<=0 ) RETURN
      IF ( Incx==Incy.AND.Incx>0 ) THEN
!
!     Code for equal, positive, non-unit increments.
!
        ns = N*Incx
        DO i = 1 , ns , Incx
          DSDOT = DSDOT + DBLE(Sx(i))*DBLE(Sy(i))
        ENDDO
        GOTO 99999
      ENDIF
!
!     Code for unequal or nonpositive increments.
!
      kx = 1
      ky = 1
      IF ( Incx<0 ) kx = 1 + (1-N)*Incx
      IF ( Incy<0 ) ky = 1 + (1-N)*Incy
      DO i = 1 , N
        DSDOT = DSDOT + DBLE(Sx(kx))*DBLE(Sy(ky))
        kx = kx + Incx
        ky = ky + Incy
      ENDDO
      RETURN
99999 END FUNCTION DSDOT
