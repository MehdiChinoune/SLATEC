!*==SDSDOT.f90  processed by SPAG 6.72Dc at 10:58 on  6 Feb 2019
!DECK SDSDOT
      REAL FUNCTION SDSDOT(N,Sb,Sx,Incx,Sy,Incy)
      IMPLICIT NONE
!*--SDSDOT5
!*** Start of declarations inserted by SPAG
      INTEGER i , Incx , Incy , kx , ky , N , ns
!*** End of declarations inserted by SPAG
!***BEGIN PROLOGUE  SDSDOT
!***PURPOSE  Compute the inner product of two vectors with extended
!            precision accumulation.
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1A4
!***TYPE      SINGLE PRECISION (SDSDOT-S, CDCDOT-C)
!***KEYWORDS  BLAS, DOT PRODUCT, INNER PRODUCT, LINEAR ALGEBRA, VECTOR
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
!       SB  single precision scalar to be added to inner product
!       SX  single precision vector with N elements
!     INCX  storage spacing between elements of SX
!       SY  single precision vector with N elements
!     INCY  storage spacing between elements of SY
!
!     --Output--
!   SDSDOT  single precision dot product (SB if N .LE. 0)
!
!     Returns S.P. result with dot product accumulated in D.P.
!     SDSDOT = SB + sum for I = 0 to N-1 of SX(LX+I*INCX)*SY(LY+I*INCY),
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
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  SDSDOT
      REAL Sx(*) , Sy(*) , Sb
      DOUBLE PRECISION dsdot
!***FIRST EXECUTABLE STATEMENT  SDSDOT
      dsdot = Sb
      IF ( N>0 ) THEN
        IF ( Incx==Incy.AND.Incx>0 ) THEN
!
!     Code for equal and positive increments.
!
          ns = N*Incx
          DO i = 1 , ns , Incx
            dsdot = dsdot + DBLE(Sx(i))*DBLE(Sy(i))
          ENDDO
          SDSDOT = dsdot
          GOTO 99999
        ELSE
!
!     Code for unequal or nonpositive increments.
!
          kx = 1
          ky = 1
          IF ( Incx<0 ) kx = 1 + (1-N)*Incx
          IF ( Incy<0 ) ky = 1 + (1-N)*Incy
          DO i = 1 , N
            dsdot = dsdot + DBLE(Sx(kx))*DBLE(Sy(ky))
            kx = kx + Incx
            ky = ky + Incy
          ENDDO
        ENDIF
      ENDIF
      SDSDOT = dsdot
      RETURN
99999 END FUNCTION SDSDOT
