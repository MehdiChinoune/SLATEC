!** CSSCAL
SUBROUTINE CSSCAL(N,Sa,Cx,Incx)
  !>
  !***
  !  Scale a complex vector.
  !***
  ! **Library:**   SLATEC (BLAS)
  !***
  ! **Category:**  D1A6
  !***
  ! **Type:**      COMPLEX (CSSCAL-C)
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
  !       SA  single precision scale factor
  !       CX  complex vector with N elements
  !     INCX  storage spacing between elements of CX
  !
  !     --Output--
  !       CX  scaled result (unchanged if N .LE. 0)
  !
  !     Replace complex CX by (single precision SA) * (complex CX)
  !     For I = 0 to N-1, replace CX(IX+I*INCX) with  SA * CX(IX+I*INCX),
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
  REAL Sa
  INTEGER i, Incx, ix, N
  !* FIRST EXECUTABLE STATEMENT  CSSCAL
  IF ( N<=0 ) RETURN
  !
  IF ( Incx==1 ) THEN
    !
    !     Code for increment equal to 1.
    !
    DO i = 1, N
      Cx(i) = Sa*Cx(i)
    END DO
    RETURN
  END IF
  !
  !     Code for increment not equal to 1.
  !
  ix = 1
  IF ( Incx<0 ) ix = (-N+1)*Incx + 1
  DO i = 1, N
    Cx(ix) = Sa*Cx(ix)
    ix = ix + Incx
  END DO
  RETURN
END SUBROUTINE CSSCAL
