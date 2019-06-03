!** SASUM
REAL(SP) FUNCTION SASUM(N,Sx,Incx)
  !>
  !  Compute the sum of the magnitudes of the elements of a
  !            vector.
  !***
  ! **Library:**   SLATEC (BLAS)
  !***
  ! **Category:**  D1A3A
  !***
  ! **Type:**      SINGLE PRECISION (SASUM-S, DASUM-D, SCASUM-C)
  !***
  ! **Keywords:**  BLAS, LINEAR ALGEBRA, SUM OF MAGNITUDES OF A VECTOR
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
  !        N  number of elements in input vector(S)
  !       SX  single precision vector with N elements
  !     INCX  storage spacing between elements of SX
  !
  !     --Output--
  !    SASUM  single precision result (zero if N .LE. 0)
  !
  !     Returns sum of magnitudes of single precision SX.
  !     SASUM = sum from 0 to N-1 of ABS(SX(IX+I*INCX)),
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

  REAL(SP) Sx(*)
  INTEGER i, Incx, ix, m, mp1, N
  !* FIRST EXECUTABLE STATEMENT  SASUM
  SASUM = 0.0E0
  IF ( N<=0 ) RETURN
  !
  IF ( Incx==1 ) THEN
    !
    !     Code for increment equal to 1.
    !
    !     Clean-up loop so remaining vector length is a multiple of 6.
    !
    m = MOD(N,6)
    IF ( m/=0 ) THEN
      DO i = 1, m
        SASUM = SASUM + ABS(Sx(i))
      END DO
      IF ( N<6 ) RETURN
    END IF
    mp1 = m + 1
    DO i = mp1, N, 6
      SASUM = SASUM + ABS(Sx(i)) + ABS(Sx(i+1)) + ABS(Sx(i+2))&
        + ABS(Sx(i+3)) + ABS(Sx(i+4)) + ABS(Sx(i+5))
    END DO
  ELSE
    !
    !     Code for increment not equal to 1.
    !
    ix = 1
    IF ( Incx<0 ) ix = (-N+1)*Incx + 1
    DO i = 1, N
      SASUM = SASUM + ABS(Sx(ix))
      ix = ix + Incx
    END DO
    RETURN
  END IF
END FUNCTION SASUM
