!DECK CDOTU
COMPLEX FUNCTION CDOTU(N,Cx,Incx,Cy,Incy)
  IMPLICIT NONE
  INTEGER i, Incx, Incy, kx, ky, N, ns
  !***BEGIN PROLOGUE  CDOTU
  !***PURPOSE  Compute the inner product of two vectors.
  !***LIBRARY   SLATEC (BLAS)
  !***CATEGORY  D1A4
  !***TYPE      COMPLEX (SDOT-S, DDOT-D, CDOTU-C)
  !***KEYWORDS  BLAS, INNER PRODUCT, LINEAR ALGEBRA, VECTOR
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
  !       CX  complex vector with N elements
  !     INCX  storage spacing between elements of CX
  !       CY  complex vector with N elements
  !     INCY  storage spacing between elements of CY
  !
  !     --Output--
  !    CDOTU  complex result (zero if N .LE. 0)
  !
  !     Returns the dot product of complex CX and CY, no conjugation
  !     CDOTU = SUM for I = 0 to N-1 of  CX(LX+I*INCX) * CY(LY+I*INCY),
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
  !***END PROLOGUE  CDOTU
  COMPLEX Cx(*), Cy(*)
  !***FIRST EXECUTABLE STATEMENT  CDOTU
  CDOTU = (0.0,0.0)
  IF ( N<=0 ) RETURN
  IF ( Incx==Incy.AND.Incx>0 ) THEN
    !
    !     Code for equal, positive increments.
    !
    ns = N*Incx
    DO i = 1, ns, Incx
      CDOTU = CDOTU + Cx(i)*Cy(i)
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
  DO i = 1, N
    CDOTU = CDOTU + Cx(kx)*Cy(ky)
    kx = kx + Incx
    ky = ky + Incy
  ENDDO
  RETURN
  99999 CONTINUE
  END FUNCTION CDOTU
