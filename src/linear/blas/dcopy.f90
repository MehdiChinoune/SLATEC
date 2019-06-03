!** DCOPY
SUBROUTINE DCOPY(N,Dx,Incx,Dy,Incy)
  !>
  !  Copy a vector.
  !***
  ! **Library:**   SLATEC (BLAS)
  !***
  ! **Category:**  D1A5
  !***
  ! **Type:**      DOUBLE PRECISION (SCOPY-S, DCOPY-D, CCOPY-C, ICOPY-I)
  !***
  ! **Keywords:**  BLAS, COPY, LINEAR ALGEBRA, VECTOR
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
  !       DY  double precision vector with N elements
  !     INCY  storage spacing between elements of DY
  !
  !     --Output--
  !       DY  copy of vector DX (unchanged if N .LE. 0)
  !
  !     Copy double precision DX to double precision DY.
  !     For I = 0 to N-1, copy DX(LX+I*INCX) to DY(LY+I*INCY),
  !     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
  !     defined in a similar way using INCY.
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
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  INTEGER i, Incx, Incy, ix, iy, m, mp1, N, ns
  REAL(DP) :: Dx(*), Dy(*)
  !* FIRST EXECUTABLE STATEMENT  DCOPY
  IF ( N<=0 ) RETURN
  IF ( Incx==Incy ) THEN
    IF ( Incx<1 ) THEN
    ELSEIF ( Incx==1 ) THEN
      !
      !     Code for both increments equal to 1.
      !
      !     Clean-up loop so remaining vector length is a multiple of 7.
      !
      m = MOD(N,7)
      IF ( m/=0 ) THEN
        DO i = 1, m
          Dy(i) = Dx(i)
        END DO
        IF ( N<7 ) RETURN
      END IF
      GOTO 100
    ELSE
      !
      !     Code for equal, positive, non-unit increments.
      !
      ns = N*Incx
      DO i = 1, ns, Incx
        Dy(i) = Dx(i)
      END DO
      RETURN
    END IF
  END IF
  !
  !     Code for unequal or nonpositive increments.
  !
  ix = 1
  iy = 1
  IF ( Incx<0 ) ix = (-N+1)*Incx + 1
  IF ( Incy<0 ) iy = (-N+1)*Incy + 1
  DO i = 1, N
    Dy(iy) = Dx(ix)
    ix = ix + Incx
    iy = iy + Incy
  END DO
  RETURN
  100  mp1 = m + 1
  DO i = mp1, N, 7
    Dy(i) = Dx(i)
    Dy(i+1) = Dx(i+1)
    Dy(i+2) = Dx(i+2)
    Dy(i+3) = Dx(i+3)
    Dy(i+4) = Dx(i+4)
    Dy(i+5) = Dx(i+5)
    Dy(i+6) = Dx(i+6)
  END DO
  RETURN
END SUBROUTINE DCOPY
