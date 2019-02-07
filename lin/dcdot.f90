!*==DCDOT.f90  processed by SPAG 6.72Dc at 10:58 on  6 Feb 2019
!DECK DCDOT
SUBROUTINE DCDOT(N,Fm,Cx,Incx,Cy,Incy,Dcr,Dci)
  IMPLICIT NONE
  !*--DCDOT5
  !***BEGIN PROLOGUE  DCDOT
  !***PURPOSE  Compute the inner product of two vectors with extended
  !            precision accumulation and result.
  !***LIBRARY   SLATEC (BLAS)
  !***CATEGORY  D1A4
  !***TYPE      COMPLEX (DSDOT-D, DCDOT-C)
  !***KEYWORDS  BLAS, COMPLEX VECTORS, DOT PRODUCT, INNER PRODUCT,
  !             LINEAR ALGEBRA, VECTOR
  !***AUTHOR  (UNKNOWN)
  !***DESCRIPTION
  !
  !    Compute the dot product of 2 complex vectors, CX and CY, e.g.
  !    CX DOT CY, or, CXconjugate DOT CY.  The real and imaginary
  !    parts of CX and CY are converted to double precision, the dot
  !    product accumulation is done in double precision and the output
  !    is given as 2 double precision numbers, corresponding to the real
  !    and imaginary part of the result.
  !     Input
  !      N:  Number of complex components of CX and CY.
  !      FM: =+1.0   compute CX DOT CY.
  !          =-1.0   compute CXconjugate DOT CY.
  !      CX(N):
  !      CY(N):  Complex arrays of length N.
  !      INCX:(Integer)   Spacing of elements of CX to use
  !      INCY:(Integer)   Spacing of elements of CY to use.
  !     Output
  !      DCR:(Double Precision) Real part of dot product.
  !      DCI:(Double Precision) Imaginary part of dot product.
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   790101  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  DCDOT
  INTEGER i, Incx, Incy, kx, ky, N
  COMPLEX Cx(*), Cy(*)
  REAL(8) :: Dcr, Dci, dt1, dt2, dt3, dt4, Fm
  !***FIRST EXECUTABLE STATEMENT  DCDOT
  Dcr = 0.0D0
  Dci = 0.0D0
  IF ( N>0 ) THEN
    !
    kx = 1
    ky = 1
    IF ( Incx<0 ) kx = 1 + (1-N)*Incx
    IF ( Incy<0 ) ky = 1 + (1-N)*Incy
    DO i = 1, N
      dt1 = REAL(REAL(Cx(kx)), 8)
      dt2 = REAL(REAL(Cy(ky)), 8)
      dt3 = REAL(AIMAG(Cx(kx)), 8)
      dt4 = REAL(AIMAG(Cy(ky)), 8)
      Dcr = Dcr + (dt1*dt2) - Fm*(dt3*dt4)
      Dci = Dci + (dt1*dt4) + Fm*(dt3*dt2)
      kx = kx + Incx
      ky = ky + Incy
    ENDDO
  ENDIF
END SUBROUTINE DCDOT
