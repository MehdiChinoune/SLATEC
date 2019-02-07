!*==CDCDOT.f90  processed by SPAG 6.72Dc at 10:58 on  6 Feb 2019
!DECK CDCDOT
COMPLEX FUNCTION CDCDOT(N,Cb,Cx,Incx,Cy,Incy)
  IMPLICIT NONE
  !*--CDCDOT5
  !***BEGIN PROLOGUE  CDCDOT
  !***PURPOSE  Compute the inner product of two vectors with extended
  !            precision accumulation.
  !***LIBRARY   SLATEC (BLAS)
  !***CATEGORY  D1A4
  !***TYPE      COMPLEX (SDSDOT-S, CDCDOT-C)
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
  !       CB  complex scalar to be added to inner product
  !       CX  complex vector with N elements
  !     INCX  storage spacing between elements of CX
  !       CY  complex vector with N elements
  !     INCY  storage spacing between elements of CY
  !
  !     --Output--
  !   CDCDOT  complex dot product (CB if N .LE. 0)
  !
  !     Returns complex result with dot product accumulated in D.P.
  !     CDCDOT = CB + sum for I = 0 to N-1 of CX(LX+I*INCY)*CY(LY+I*INCY)
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
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  CDCDOT
  INTEGER N , Incx , Incy , i , kx , ky
  COMPLEX Cx(*) , Cy(*) , Cb
  REAL(8) :: dsdotr , dsdoti , dt1 , dt2 , dt3 , dt4
  !***FIRST EXECUTABLE STATEMENT  CDCDOT
  dsdotr = REAL(REAL(Cb), 8)
  dsdoti = REAL(AIMAG(Cb), 8)
  IF ( N>0 ) THEN
    kx = 1
    ky = 1
    IF ( Incx<0 ) kx = 1 + (1-N)*Incx
    IF ( Incy<0 ) ky = 1 + (1-N)*Incy
    DO i = 1 , N
      dt1 = REAL(REAL(Cx(kx)), 8)
      dt2 = REAL(REAL(Cy(ky)), 8)
      dt3 = REAL(AIMAG(Cx(kx)), 8)
      dt4 = REAL(AIMAG(Cy(ky)), 8)
      dsdotr = dsdotr + (dt1*dt2) - (dt3*dt4)
      dsdoti = dsdoti + (dt1*dt4) + (dt3*dt2)
      kx = kx + Incx
      ky = ky + Incy
    ENDDO
  ENDIF
  CDCDOT = CMPLX(REAL(dsdotr),REAL(dsdoti))
END FUNCTION CDCDOT
