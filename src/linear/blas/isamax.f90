!** ISAMAX
INTEGER FUNCTION ISAMAX(N,Sx,Incx)
  !>
  !***
  !  Find the smallest index of that component of a vector
  !            having the maximum magnitude.
  !***
  ! **Library:**   SLATEC (BLAS)
  !***
  ! **Category:**  D1A2
  !***
  ! **Type:**      SINGLE PRECISION (ISAMAX-S, IDAMAX-D, ICAMAX-C)
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
  !       SX  single precision vector with N elements
  !     INCX  storage spacing between elements of SX
  !
  !     --Output--
  !   ISAMAX  smallest index (zero if N .LE. 0)
  !
  !     Find smallest index of maximum magnitude of single precision SX.
  !     ISAMAX = first I, I = 1 to N, to maximize  ABS(SX(IX+(I-1)*INCX)),
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
  !   920618  Slight restructuring of code.  (RWC, WRB)

  REAL Sx(*), smax, xmag
  INTEGER i, Incx, ix, N
  !* FIRST EXECUTABLE STATEMENT  ISAMAX
  ISAMAX = 0
  IF ( N<=0 ) RETURN
  ISAMAX = 1
  IF ( N==1 ) RETURN
  !
  IF ( Incx==1 ) THEN
    !
    !     Code for increments equal to 1.
    !
    smax = ABS(Sx(1))
    DO i = 2, N
      xmag = ABS(Sx(i))
      IF ( xmag>smax ) THEN
        ISAMAX = i
        smax = xmag
      END IF
    END DO
    RETURN
  END IF
  !
  !     Code for increment not equal to 1.
  !
  ix = 1
  IF ( Incx<0 ) ix = (-N+1)*Incx + 1
  smax = ABS(Sx(ix))
  ix = ix + Incx
  DO i = 2, N
    xmag = ABS(Sx(ix))
    IF ( xmag>smax ) THEN
      ISAMAX = i
      smax = xmag
    END IF
    ix = ix + Incx
  END DO
  RETURN
END FUNCTION ISAMAX
