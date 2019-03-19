!** ICAMAX
INTEGER FUNCTION ICAMAX(N,Cx,Incx)
  IMPLICIT NONE
  !>
  !***
  !  Find the smallest index of the component of a complex
  !            vector having the maximum sum of magnitudes of real
  !            and imaginary parts.
  !***
  ! **Library:**   SLATEC (BLAS)
  !***
  ! **Category:**  D1A2
  !***
  ! **Type:**      COMPLEX (ISAMAX-S, IDAMAX-D, ICAMAX-C)
  !***
  ! **Keywords:**  BLAS, LINEAR ALGEBRA, MAXIMUM COMPONENT, VECTOR
  !***
  ! **Author:**  Lawson, C. L., (JPL)
  !           Hanson, R. J., (SNLA)
  !           Kincaid, D. R., (U. of Texas)
  !           Krogh, F. T., (JPL)
  !***
  ! **Description:**
  !
  !                B L A S  Subprogram
  !    Description of Parameters
  !
  !     --Input--
  !        N  number of elements in input vector(s)
  !       CX  complex vector with N elements
  !     INCX  storage spacing between elements of CX
  !
  !     --Output--
  !   ICAMAX  smallest index (zero if N .LE. 0)
  !
  !     Returns the smallest index of the component of CX having the
  !     largest sum of magnitudes of real and imaginary parts.
  !     ICAMAX = first I, I = 1 to N, to maximize
  !     ABS(REAL(CX(IX+(I-1)*INCX))) + ABS(IMAG(CX(IX+(I-1)*INCX))),
  !     where IX = 1 if INCX .GE. 0, else IX = 1+(1-N)*INCX.
  !
  !***
  ! **References:**  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
  !                 Krogh, Basic linear algebra subprograms for Fortran
  !                 usage, Algorithm No. 539, Transactions on Mathematical
  !                 Software 5, 3 (September 1979), pp. 308-323.
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   791001  DATE WRITTEN
  !   861211  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900821  Modified to correct problem with a negative increment.
  !           (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  
  COMPLEX Cx(*)
  REAL smax, xmag
  INTEGER i, Incx, ix, N
  REAL, EXTERNAL :: CABS1
  !* FIRST EXECUTABLE STATEMENT  ICAMAX
  ICAMAX = 0
  IF ( N<=0 ) RETURN
  ICAMAX = 1
  IF ( N==1 ) RETURN
  !
  IF ( Incx==1 ) THEN
    !
    !     Code for increment equal to 1.
    !
    smax = CABS1(Cx(1))
    DO i = 2, N
      xmag = CABS1(Cx(i))
      IF ( xmag>smax ) THEN
        ICAMAX = i
        smax = xmag
      ENDIF
    ENDDO
    RETURN
  ENDIF
  !
  !     Code for increment not equal to 1.
  !
  ix = 1
  IF ( Incx<0 ) ix = (-N+1)*Incx + 1
  smax = CABS1(Cx(ix))
  ix = ix + Incx
  DO i = 2, N
    xmag = CABS1(Cx(ix))
    IF ( xmag>smax ) THEN
      ICAMAX = i
      smax = xmag
    ENDIF
    ix = ix + Incx
  ENDDO
  RETURN
END FUNCTION ICAMAX
