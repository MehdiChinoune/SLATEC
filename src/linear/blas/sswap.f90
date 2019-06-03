!** SSWAP
SUBROUTINE SSWAP(N,Sx,Incx,Sy,Incy)
  !>
  !  Interchange two vectors.
  !***
  ! **Library:**   SLATEC (BLAS)
  !***
  ! **Category:**  D1A5
  !***
  ! **Type:**      SINGLE PRECISION (SSWAP-S, DSWAP-D, CSWAP-C, ISWAP-I)
  !***
  ! **Keywords:**  BLAS, INTERCHANGE, LINEAR ALGEBRA, VECTOR
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
  !       SY  single precision vector with N elements
  !     INCY  storage spacing between elements of SY
  !
  !     --Output--
  !       SX  input vector SY (unchanged if N .LE. 0)
  !       SY  input vector SX (unchanged if N .LE. 0)
  !
  !     Interchange single precision SX and single precision SY.
  !     For I = 0 to N-1, interchange  SX(LX+I*INCX) and SY(LY+I*INCY),
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
  REAL(SP) Sx(*), Sy(*), stemp1, stemp2, stemp3
  !* FIRST EXECUTABLE STATEMENT  SSWAP
  IF ( N<=0 ) RETURN
  IF ( Incx==Incy ) THEN
    IF ( Incx<1 ) THEN
    ELSEIF ( Incx==1 ) THEN
      !
      !     Code for both increments equal to 1.
      !
      !     Clean-up loop so remaining vector length is a multiple of 3.
      !
      m = MOD(N,3)
      IF ( m/=0 ) THEN
        DO i = 1, m
          stemp1 = Sx(i)
          Sx(i) = Sy(i)
          Sy(i) = stemp1
        END DO
        IF ( N<3 ) RETURN
      END IF
      GOTO 100
    ELSE
      !
      !     Code for equal, positive, non-unit increments.
      !
      ns = N*Incx
      DO i = 1, ns, Incx
        stemp1 = Sx(i)
        Sx(i) = Sy(i)
        Sy(i) = stemp1
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
    stemp1 = Sx(ix)
    Sx(ix) = Sy(iy)
    Sy(iy) = stemp1
    ix = ix + Incx
    iy = iy + Incy
  END DO
  RETURN
  100  mp1 = m + 1
  DO i = mp1, N, 3
    stemp1 = Sx(i)
    stemp2 = Sx(i+1)
    stemp3 = Sx(i+2)
    Sx(i) = Sy(i)
    Sx(i+1) = Sy(i+1)
    Sx(i+2) = Sy(i+2)
    Sy(i) = stemp1
    Sy(i+1) = stemp2
    Sy(i+2) = stemp3
  END DO
  RETURN
END SUBROUTINE SSWAP
