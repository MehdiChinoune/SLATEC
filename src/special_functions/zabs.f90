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

  REAL(DP), INTENT(IN) :: Zr, Zi
  !* FIRST EXECUTABLE STATEMENT  ZABS

  ZABS = ABS( CMPLX( Zr, Zi, DP ) )

END FUNCTION ZABS