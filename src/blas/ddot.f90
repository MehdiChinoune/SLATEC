!** DDOT
REAL(8) FUNCTION DDOT(N,Dx,Incx,Dy,Incy)
  IMPLICIT NONE
  !>
  !***
  !  Compute the inner product of two vectors.
  !***
  ! **Library:**   SLATEC (BLAS)
  !***
  ! **Category:**  D1A4
  !***
  ! **Type:**      DOUBLE PRECISION (SDOT-S, DDOT-D, CDOTU-C)
  !***
  ! **Keywords:**  BLAS, INNER PRODUCT, LINEAR ALGEBRA, VECTOR
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
  !     DDOT  double precision dot product (zero if N .LE. 0)
  !
  !     Returns the dot product of double precision DX and DY.
  !     DDOT = sum for I = 0 to N-1 of  DX(LX+I*INCX) * DY(LY+I*INCY),
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
  REAL(8) :: Dx(*), Dy(*)
  !* FIRST EXECUTABLE STATEMENT  DDOT
  DDOT = 0.0D0
  IF ( N<=0 ) RETURN
  IF ( Incx==Incy ) THEN
    IF ( Incx<1 ) THEN
    ELSEIF ( Incx==1 ) THEN
      !
      !     Code for both increments equal to 1.
      !
      !     Clean-up loop so remaining vector length is a multiple of 5.
      !
      m = MOD(N,5)
      IF ( m/=0 ) THEN
        DO i = 1, m
          DDOT = DDOT + Dx(i)*Dy(i)
        ENDDO
        IF ( N<5 ) RETURN
      ENDIF
      GOTO 100
    ELSE
      !
      !     Code for equal, positive, non-unit increments.
      !
      ns = N*Incx
      DO i = 1, ns, Incx
        DDOT = DDOT + Dx(i)*Dy(i)
      ENDDO
      RETURN
    ENDIF
  ENDIF
  !
  !     Code for unequal or nonpositive increments.
  !
  ix = 1
  iy = 1
  IF ( Incx<0 ) ix = (-N+1)*Incx + 1
  IF ( Incy<0 ) iy = (-N+1)*Incy + 1
  DO i = 1, N
    DDOT = DDOT + Dx(ix)*Dy(iy)
    ix = ix + Incx
    iy = iy + Incy
  ENDDO
  RETURN
  100  mp1 = m + 1
  DO i = mp1, N, 5
    DDOT = DDOT + Dx(i)*Dy(i) + Dx(i+1)*Dy(i+1) + Dx(i+2)*Dy(i+2) + Dx(i+3)&
      *Dy(i+3) + Dx(i+4)*Dy(i+4)
  ENDDO
  RETURN
END FUNCTION DDOT
