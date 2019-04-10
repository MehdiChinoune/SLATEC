!** SCASUM
REAL FUNCTION SCASUM(N,Cx,Incx)
  IMPLICIT NONE
  !>
  !***
  !  Compute the sum of the magnitudes of the real and
  !            imaginary elements of a complex vector.
  !***
  ! **Library:**   SLATEC (BLAS)
  !***
  ! **Category:**  D1A3A
  !***
  ! **Type:**      COMPLEX (SASUM-S, DASUM-D, SCASUM-C)
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
  !        N  number of elements in input vector(s)
  !       CX  complex vector with N elements
  !     INCX  storage spacing between elements of CX
  !
  !     --Output--
  !   SCASUM  single precision result (zero if N .LE. 0)
  !
  !     Returns sums of magnitudes of real and imaginary parts of
  !     components of CX.  Note that this is not the L1 norm of CX.
  !     CASUM = sum from 0 to N-1 of ABS(REAL(CX(IX+I*INCX))) +
  !             ABS(IMAG(CX(IX+I*INCX))),
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

  COMPLEX Cx(*)
  INTEGER i, Incx, ix, N
  !* FIRST EXECUTABLE STATEMENT  SCASUM
  SCASUM = 0.0E0
  IF ( N<=0 ) RETURN
  !
  IF ( Incx==1 ) THEN
    !
    !     Code for increment equal to 1.
    !
    DO i = 1, N
      SCASUM = SCASUM + ABS(REAL(Cx(i))) + ABS(AIMAG(Cx(i)))
    END DO
    RETURN
  END IF
  !
  !     Code for increment not equal to 1.
  !
  ix = 1
  IF ( Incx<0 ) ix = (-N+1)*Incx + 1
  DO i = 1, N
    SCASUM = SCASUM + ABS(REAL(Cx(ix))) + ABS(AIMAG(Cx(ix)))
    ix = ix + Incx
  END DO
  RETURN
END FUNCTION SCASUM