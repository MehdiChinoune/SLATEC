!** DAXPY
SUBROUTINE DAXPY(N,Da,Dx,Incx,Dy,Incy)
  IMPLICIT NONE
  !>
  !***
  !  Compute a constant times a vector plus a vector.
  !***
  ! **Library:**   SLATEC (BLAS)
  !***
  ! **Category:**  D1A7
  !***
  ! **Type:**      DOUBLE PRECISION (SAXPY-S, DAXPY-D, CAXPY-C)
  !***
  ! **Keywords:**  BLAS, LINEAR ALGEBRA, TRIAD, VECTOR
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
  !       DA  double precision scalar multiplier
  !       DX  double precision vector with N elements
  !     INCX  storage spacing between elements of DX
  !       DY  double precision vector with N elements
  !     INCY  storage spacing between elements of DY
  !
  !     --Output--
  !       DY  double precision result (unchanged if N .LE. 0)
  !
  !     Overwrite double precision DY with double precision DA*DX + DY.
  !     For I = 0 to N-1, replace  DY(LY+I*INCY) with DA*DX(LX+I*INCX) +
  !       DY(LY+I*INCY),
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
  REAL(8) :: Dx(*), Dy(*), Da
  !* FIRST EXECUTABLE STATEMENT  DAXPY
  IF ( N<=0.OR.Da==0.0D0 ) RETURN
  IF ( Incx==Incy ) THEN
    IF ( Incx<1 ) THEN
    ELSEIF ( Incx==1 ) THEN
      !
      !     Code for both increments equal to 1.
      !
      !     Clean-up loop so remaining vector length is a multiple of 4.
      !
      m = MOD(N,4)
      IF ( m/=0 ) THEN
        DO i = 1, m
          Dy(i) = Dy(i) + Da*Dx(i)
        END DO
        IF ( N<4 ) RETURN
      END IF
      GOTO 100
    ELSE
      !
      !     Code for equal, positive, non-unit increments.
      !
      ns = N*Incx
      DO i = 1, ns, Incx
        Dy(i) = Da*Dx(i) + Dy(i)
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
    Dy(iy) = Dy(iy) + Da*Dx(ix)
    ix = ix + Incx
    iy = iy + Incy
  END DO
  RETURN
  100  mp1 = m + 1
  DO i = mp1, N, 4
    Dy(i) = Dy(i) + Da*Dx(i)
    Dy(i+1) = Dy(i+1) + Da*Dx(i+1)
    Dy(i+2) = Dy(i+2) + Da*Dx(i+2)
    Dy(i+3) = Dy(i+3) + Da*Dx(i+3)
  END DO
  RETURN
END SUBROUTINE DAXPY
