!** ZABS
REAL(DP) ELEMENTAL FUNCTION ZABS(Zr,Zi)
  !> Subsidiary to ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZAIRY and
  !            ZBIRY
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      ALL (ZABS-A)
  !***
  ! **Author:**  Amos, D. E., (SNL)
  !***
  ! **Description:**
  !
  !     ZABS COMPUTES THE ABSOLUTE VALUE OR MAGNITUDE OF A DOUBLE
  !     PRECISION COMPLEX VARIABLE CMPLX(ZR,ZI)
  !
  !***
  ! **See also:**  ZAIRY, ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZBIRY
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)
  USE service, ONLY : D1MACH

  REAL(DP), INTENT(IN) :: Zr, Zi
  REAL(DP), PARAMETER :: sqrt_tiny = SQRT( D1MACH(1) ), sqrt_huge = SQRT( D1MACH(2) )
  !* FIRST EXECUTABLE STATEMENT  ZABS

  IF( ABS(Zr)<sqrt_tiny .AND. ABS(Zi)<sqrt_tiny ) THEN
    ZABS = ABS( CMPLX( Zr/sqrt_tiny, Zi/sqrt_tiny, DP ) ) * sqrt_tiny
  ELSEIF( ABS(Zr)>sqrt_huge .OR. ABS(Zi)>sqrt_huge ) THEN
    ZABS = ABS( CMPLX( Zr/sqrt_huge, Zi/sqrt_huge, DP ) ) * sqrt_huge
  ELSE
    ZABS = ABS( CMPLX( Zr, Zi, DP ) )
  END IF

END FUNCTION ZABS