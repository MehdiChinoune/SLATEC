!** ZABS
REAL(8) FUNCTION ZABS(Zr,Zi)
  !>
  !  Subsidiary to ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZAIRY and
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

  REAL(8) :: Zr, Zi, u, v, q, s
  !* FIRST EXECUTABLE STATEMENT  ZABS
  u = ABS(Zr)
  v = ABS(Zi)
  s = u + v
  !-----------------------------------------------------------------------
  !     S*1.0D0 MAKES AN UNNORMALIZED UNDERFLOW ON CDC MACHINES INTO A
  !     TRUE FLOATING ZERO
  !-----------------------------------------------------------------------
  s = s*1.0D+0
  IF ( s==0.0D+0 ) THEN
    ZABS = 0.0D+0
    RETURN
  ELSEIF ( u<=v ) THEN
    q = u/v
    ZABS = v*SQRT(1.D+0+q*q)
    RETURN
  END IF
  q = v/u
  ZABS = u*SQRT(1.D+0+q*q)
  RETURN
END FUNCTION ZABS
