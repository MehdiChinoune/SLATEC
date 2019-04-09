!** CSWAP
SUBROUTINE CSWAP(N,Cx,Incx,Cy,Incy)
  IMPLICIT NONE
  !>
  !***
  !  Interchange two vectors.
  !***
  ! **Library:**   SLATEC (BLAS)
  !***
  ! **Category:**  D1A5
  !***
  ! **Type:**      COMPLEX (SSWAP-S, DSWAP-D, CSWAP-C, ISWAP-I)
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
  !       CX  complex vector with N elements
  !     INCX  storage spacing between elements of CX
  !       CY  complex vector with N elements
  !     INCY  storage spacing between elements of CY
  !
  !     --Output--
  !       CX  input vector CY (unchanged if N .LE. 0)
  !       CY  input vector CX (unchanged if N .LE. 0)
  !
  !     Interchange complex CX and complex CY
  !     For I = 0 to N-1, interchange  CX(LX+I*INCX) and CY(LY+I*INCY),
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
  COMPLEX Cx(*), Cy(*), ctemp
  !* FIRST EXECUTABLE STATEMENT  CSWAP
  IF ( N<=0 ) RETURN
  IF ( Incx==Incy.AND.Incx>0 ) THEN
    !
    !     Code for equal, positive, non-unit increments.
    !
    ns = N*Incx
    DO i = 1, ns, Incx
      ctemp = Cx(i)
      Cx(i) = Cy(i)
      Cy(i) = ctemp
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
    ctemp = Cx(kx)
    Cx(kx) = Cy(ky)
    Cy(ky) = ctemp
    kx = kx + Incx
    ky = ky + Incy
  END DO
  RETURN
END SUBROUTINE CSWAP
