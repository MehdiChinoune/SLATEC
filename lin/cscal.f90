!DECK CSCAL
SUBROUTINE CSCAL(N,Ca,Cx,Incx)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  CSCAL
  !***PURPOSE  Multiply a vector by a constant.
  !***LIBRARY   SLATEC (BLAS)
  !***CATEGORY  D1A6
  !***TYPE      COMPLEX (SSCAL-S, DSCAL-D, CSCAL-C)
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
  !       CA  complex scale factor
  !       CX  complex vector with N elements
  !     INCX  storage spacing between elements of CX
  !
  !     --Output--
  !       CX  complex result (unchanged if N .LE. 0)
  !
  !     Replace complex CX by complex CA*CX.
  !     For I = 0 to N-1, replace CX(IX+I*INCX) with CA*CX(IX+I*INCX),
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
  !***END PROLOGUE  CSCAL
  COMPLEX Ca, Cx(*)
  INTEGER i, Incx, ix, N
  !***FIRST EXECUTABLE STATEMENT  CSCAL
  IF ( N<=0 ) RETURN
  !
  IF ( Incx==1 ) THEN
    !
    !     Code for increment equal to 1.
    !
    DO i = 1, N
      Cx(i) = Ca*Cx(i)
    ENDDO
    RETURN
  ENDIF
  !
  !     Code for increment not equal to 1.
  !
  ix = 1
  IF ( Incx<0 ) ix = (-N+1)*Incx + 1
  DO i = 1, N
    Cx(ix) = Ca*Cx(ix)
    ix = ix + Incx
  ENDDO
  RETURN
END SUBROUTINE CSCAL
