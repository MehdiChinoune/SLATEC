!DECK SDOT
REAL FUNCTION SDOT(N,Sx,Incx,Sy,Incy)
  IMPLICIT NONE
  INTEGER i, Incx, Incy, ix, iy, m, mp1, N, ns
  !***BEGIN PROLOGUE  SDOT
  !***PURPOSE  Compute the inner product of two vectors.
  !***LIBRARY   SLATEC (BLAS)
  !***CATEGORY  D1A4
  !***TYPE      SINGLE PRECISION (SDOT-S, DDOT-D, CDOTU-C)
  !***KEYWORDS  BLAS, INNER PRODUCT, LINEAR ALGEBRA, VECTOR
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
  !       SX  single precision vector with N elements
  !     INCX  storage spacing between elements of SX
  !       SY  single precision vector with N elements
  !     INCY  storage spacing between elements of SY
  !
  !     --Output--
  !     SDOT  single precision dot product (zero if N .LE. 0)
  !
  !     Returns the dot product of single precision SX and SY.
  !     SDOT = sum for I = 0 to N-1 of  SX(LX+I*INCX) * SY(LY+I*INCY),
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
  !***END PROLOGUE  SDOT
  REAL Sx(*), Sy(*)
  !***FIRST EXECUTABLE STATEMENT  SDOT
  SDOT = 0.0E0
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
          SDOT = SDOT + Sx(i)*Sy(i)
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
        SDOT = SDOT + Sx(i)*Sy(i)
      ENDDO
      GOTO 99999
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
    SDOT = SDOT + Sx(ix)*Sy(iy)
    ix = ix + Incx
    iy = iy + Incy
  ENDDO
  RETURN
  100  mp1 = m + 1
  DO i = mp1, N, 5
    SDOT = SDOT + Sx(i)*Sy(i) + Sx(i+1)*Sy(i+1) + Sx(i+2)*Sy(i+2) + Sx(i+3)&
      *Sy(i+3) + Sx(i+4)*Sy(i+4)
  ENDDO
  RETURN
  99999 CONTINUE
  END FUNCTION SDOT
