!** CAXPY
SUBROUTINE CAXPY(N,Ca,Cx,Incx,Cy,Incy)
  !>
  !  Compute a constant times a vector plus a vector.
  !***
  ! **Library:**   SLATEC (BLAS)
  !***
  ! **Category:**  D1A7
  !***
  ! **Type:**      COMPLEX (SAXPY-S, DAXPY-D, CAXPY-C)
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
  !       CA  complex scalar multiplier
  !       CX  complex vector with N elements
  !     INCX  storage spacing between elements of CX
  !       CY  complex vector with N elements
  !     INCY  storage spacing between elements of CY
  !
  !     --Output--
  !       CY  complex result (unchanged if N .LE. 0)
  !
  !     Overwrite complex CY with complex  CA*CX + CY.
  !     For I = 0 to N-1, replace  CY(LY+I*INCY) with CA*CX(LX+I*INCX) +
  !       CY(LY+I*INCY),
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
  !   861211  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !   920801  Removed variable CANORM.  (RWC, WRB)

  INTEGER i, Incx, Incy, kx, ky, N, ns
  COMPLEX(SP) Cx(*), Cy(*), Ca
  !* FIRST EXECUTABLE STATEMENT  CAXPY
  IF ( N<=0.OR.Ca==(0.0E0,0.0E0) ) RETURN
  IF ( Incx==Incy.AND.Incx>0 ) THEN
    !
    !     Code for equal, positive, non-unit increments.
    !
    ns = N*Incx
    DO i = 1, ns, Incx
      Cy(i) = Ca*Cx(i) + Cy(i)
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
    Cy(ky) = Cy(ky) + Ca*Cx(kx)
    kx = kx + Incx
    ky = ky + Incy
  END DO
  RETURN
END SUBROUTINE CAXPY
