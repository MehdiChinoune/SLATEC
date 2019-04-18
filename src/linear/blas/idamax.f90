!** IDAMAX
INTEGER FUNCTION IDAMAX(N,Dx,Incx)
  !>
  !  Find the smallest index of that component of a vector
  !            having the maximum magnitude.
  !***
  ! **Library:**   SLATEC (BLAS)
  !***
  ! **Category:**  D1A2
  !***
  ! **Type:**      DOUBLE PRECISION (ISAMAX-S, IDAMAX-D, ICAMAX-C)
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
  !       DX  double precision vector with N elements
  !     INCX  storage spacing between elements of DX
  !
  !     --Output--
  !   IDAMAX  smallest index (zero if N .LE. 0)
  !
  !     Find smallest index of maximum magnitude of double precision DX.
  !     IDAMAX = first I, I = 1 to N, to maximize ABS(DX(IX+(I-1)*INCX)),
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
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900821  Modified to correct problem with a negative increment.
  !           (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  REAL(8) :: Dx(*), dmax, xmag
  INTEGER i, Incx, ix, N
  !* FIRST EXECUTABLE STATEMENT  IDAMAX
  IDAMAX = 0
  IF ( N<=0 ) RETURN
  IDAMAX = 1
  IF ( N==1 ) RETURN
  !
  IF ( Incx==1 ) THEN
    !
    !     Code for increments equal to 1.
    !
    dmax = ABS(Dx(1))
    DO i = 2, N
      xmag = ABS(Dx(i))
      IF ( xmag>dmax ) THEN
        IDAMAX = i
        dmax = xmag
      END IF
    END DO
    RETURN
  END IF
  !
  !     Code for increments not equal to 1.
  !
  ix = 1
  IF ( Incx<0 ) ix = (-N+1)*Incx + 1
  dmax = ABS(Dx(ix))
  ix = ix + Incx
  DO i = 2, N
    xmag = ABS(Dx(ix))
    IF ( xmag>dmax ) THEN
      IDAMAX = i
      dmax = xmag
    END IF
    ix = ix + Incx
  END DO
  RETURN
END FUNCTION IDAMAX
