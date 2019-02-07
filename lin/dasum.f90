!*==DASUM.f90  processed by SPAG 6.72Dc at 10:58 on  6 Feb 2019
!DECK DASUM
REAL(8) FUNCTION DASUM(N,Dx,Incx)
  IMPLICIT NONE
  !*--DASUM5
  !***BEGIN PROLOGUE  DASUM
  !***PURPOSE  Compute the sum of the magnitudes of the elements of a
  !            vector.
  !***LIBRARY   SLATEC (BLAS)
  !***CATEGORY  D1A3A
  !***TYPE      DOUBLE PRECISION (SASUM-S, DASUM-D, SCASUM-C)
  !***KEYWORDS  BLAS, LINEAR ALGEBRA, SUM OF MAGNITUDES OF A VECTOR
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
  !       DX  double precision vector with N elements
  !     INCX  storage spacing between elements of DX
  !
  !     --Output--
  !    DASUM  double precision result (zero if N .LE. 0)
  !
  !     Returns sum of magnitudes of double precision DX.
  !     DASUM = sum from 0 to N-1 of ABS(DX(IX+I*INCX)),
  !     where IX = 1 if INCX .GE. 0, else IX = 1+(1-N)*INCX.
  !
  !***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
  !                 Krogh, Basic linear algebra subprograms for Fortran
  !                 usage, Algorithm No. 539, Transactions on Mathematical
  !                 Software 5, 3 (September 1979), pp. 308-323.
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   791001  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900821  Modified to correct problem with a negative increment.
  !           (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  DASUM
  REAL(8) :: Dx(*)
  INTEGER i , Incx , ix , m , mp1 , N
  !***FIRST EXECUTABLE STATEMENT  DASUM
  DASUM = 0.0D0
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
      DO i = 1 , m
        DASUM = DASUM + ABS(Dx(i))
      ENDDO
      IF ( N<6 ) RETURN
    ENDIF
    mp1 = m + 1
    DO i = mp1 , N , 6
      DASUM = DASUM + ABS(Dx(i)) + ABS(Dx(i+1)) + ABS(Dx(i+2))&
        + ABS(Dx(i+3)) + ABS(Dx(i+4)) + ABS(Dx(i+5))
    ENDDO
  ELSE
    !
    !     Code for increment not equal to 1.
    !
    ix = 1
    IF ( Incx<0 ) ix = (-N+1)*Incx + 1
    DO i = 1 , N
      DASUM = DASUM + ABS(Dx(ix))
      ix = ix + Incx
    ENDDO
    RETURN
  ENDIF
END FUNCTION DASUM
