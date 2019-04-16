!** DSCAL
SUBROUTINE DSCAL(N,Da,Dx,Incx)
  !>
  !***
  !  Multiply a vector by a constant.
  !***
  ! **Library:**   SLATEC (BLAS)
  !***
  ! **Category:**  D1A6
  !***
  ! **Type:**      DOUBLE PRECISION (SSCAL-S, DSCAL-D, CSCAL-C)
  !***
  ! **Keywords:**  BLAS, LINEAR ALGEBRA, SCALE, VECTOR
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
  !       DA  double precision scale factor
  !       DX  double precision vector with N elements
  !     INCX  storage spacing between elements of DX
  !
  !     --Output--
  !       DX  double precision result (unchanged if N.LE.0)
  !
  !     Replace double precision DX by double precision DA*DX.
  !     For I = 0 to N-1, replace DX(IX+I*INCX) with  DA * DX(IX+I*INCX),
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
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900821  Modified to correct problem with a negative increment.
  !           (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  REAL(8) :: Da, Dx(*)
  INTEGER i, Incx, ix, m, mp1, N
  !* FIRST EXECUTABLE STATEMENT  DSCAL
  IF ( N<=0 ) RETURN
  IF ( Incx==1 ) THEN
    !
    !     Code for increment equal to 1.
    !
    !     Clean-up loop so remaining vector length is a multiple of 5.
    !
    m = MOD(N,5)
    IF ( m/=0 ) THEN
      DO i = 1, m
        Dx(i) = Da*Dx(i)
      END DO
      IF ( N<5 ) RETURN
    END IF
    mp1 = m + 1
    DO i = mp1, N, 5
      Dx(i) = Da*Dx(i)
      Dx(i+1) = Da*Dx(i+1)
      Dx(i+2) = Da*Dx(i+2)
      Dx(i+3) = Da*Dx(i+3)
      Dx(i+4) = Da*Dx(i+4)
    END DO
  ELSE
    !
    !     Code for increment not equal to 1.
    !
    ix = 1
    IF ( Incx<0 ) ix = (-N+1)*Incx + 1
    DO i = 1, N
      Dx(ix) = Da*Dx(ix)
      ix = ix + Incx
    END DO
    RETURN
  END IF
END SUBROUTINE DSCAL
