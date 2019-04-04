!** DSDOT
REAL(8) FUNCTION DSDOT(N,Sx,Incx,Sy,Incy)
  IMPLICIT NONE
  !>
  !***
  !  Compute the inner product of two vectors with extended
  !            precision accumulation and result.
  !***
  ! **Library:**   SLATEC (BLAS)
  !***
  ! **Category:**  D1A4
  !***
  ! **Type:**      DOUBLE PRECISION (DSDOT-D, DCDOT-C)
  !***
  ! **Keywords:**  BLAS, COMPLEX VECTORS, DOT PRODUCT, INNER PRODUCT,
  !             LINEAR ALGEBRA, VECTOR
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
  !    DSDOT  double precision dot product (zero if N.LE.0)
  !
  !     Returns D.P. dot product accumulated in D.P., for S.P. SX and SY
  !     DSDOT = sum for I = 0 to N-1 of  SX(LX+I*INCX) * SY(LY+I*INCY),
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

  INTEGER i, Incx, Incy, kx, ky, N, ns
  REAL Sx(*), Sy(*)
  !* FIRST EXECUTABLE STATEMENT  DSDOT
  DSDOT = 0.0D0
  IF ( N<=0 ) RETURN
  IF ( Incx==Incy.AND.Incx>0 ) THEN
    !
    !     Code for equal, positive, non-unit increments.
    !
    ns = N*Incx
    DO i = 1, ns, Incx
      DSDOT = DSDOT + REAL(Sx(i), 8)*REAL(Sy(i), 8)
    END DO
    RETURN
  END IF
  !
  !     Code for unequal or nonpositive increments.
  !
  kx = 1
  ky = 1
  IF ( Incx<0 ) kx = 1 + (1-N)*Incx
  IF ( Incy<0 ) ky = 1 + (1-N)*Incy
  DO i = 1, N
    DSDOT = DSDOT + REAL(Sx(kx), 8)*REAL(Sy(ky), 8)
    kx = kx + Incx
    ky = ky + Incy
  END DO
  RETURN
END FUNCTION DSDOT
