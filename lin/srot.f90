!*==SROT.f90  processed by SPAG 6.72Dc at 10:58 on  6 Feb 2019
!DECK SROT
      SUBROUTINE SROT(N,Sx,Incx,Sy,Incy,Sc,Ss)
      IMPLICIT NONE
!*--SROT5
!*** Start of declarations inserted by SPAG
      INTEGER i , Incx , Incy , kx , ky , N , nsteps
!*** End of declarations inserted by SPAG
!***BEGIN PROLOGUE  SROT
!***PURPOSE  Apply a plane Givens rotation.
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1A8
!***TYPE      SINGLE PRECISION (SROT-S, DROT-D, CSROT-C)
!***KEYWORDS  BLAS, GIVENS ROTATION, GIVENS TRANSFORMATION,
!             LINEAR ALGEBRA, PLANE ROTATION, VECTOR
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
!       SC  element of rotation matrix
!       SS  element of rotation matrix
!
!     --Output--
!       SX  rotated vector SX (unchanged if N .LE. 0)
!       SY  rotated vector SY (unchanged if N .LE. 0)
!
!     Multiply the 2 x 2 matrix  ( SC SS) times the 2 x N matrix (SX**T)
!                                (-SS SC)                        (SY**T)
!     where **T indicates transpose.  The elements of SX are in
!     SX(LX+I*INCX), I = 0 to N-1, where LX = 1 if INCX .GE. 0, else
!     LX = 1+(1-N)*INCX, and similarly for SY using LY and INCY.
!
!***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
!                 Krogh, Basic linear algebra subprograms for Fortran
!                 usage, Algorithm No. 539, Transactions on Mathematical
!                 Software 5, 3 (September 1979), pp. 308-323.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   791001  DATE WRITTEN
!   861211  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  SROT
      REAL Sx , Sy , Sc , Ss , zero , one , w , z
      DIMENSION Sx(*) , Sy(*)
      SAVE zero , one
      DATA zero , one/0.0E0 , 1.0E0/
!***FIRST EXECUTABLE STATEMENT  SROT
      IF ( .NOT.(N<=0.OR.(Ss==zero.AND.Sc==one)) ) THEN
        IF ( Incx/=Incy.OR.Incx<=0 ) THEN
!
!     Code for unequal or nonpositive increments.
!
          kx = 1
          ky = 1
!
          IF ( Incx<0 ) kx = 1 - (N-1)*Incx
          IF ( Incy<0 ) ky = 1 - (N-1)*Incy
!
          DO i = 1 , N
            w = Sx(kx)
            z = Sy(ky)
            Sx(kx) = Sc*w + Ss*z
            Sy(ky) = -Ss*w + Sc*z
            kx = kx + Incx
            ky = ky + Incy
          ENDDO
        ELSE
!
!          Code for equal and positive increments.
!
          nsteps = Incx*N
          DO i = 1 , nsteps , Incx
            w = Sx(i)
            z = Sy(i)
            Sx(i) = Sc*w + Ss*z
            Sy(i) = -Ss*w + Sc*z
          ENDDO
        ENDIF
      ENDIF
!
      END SUBROUTINE SROT
