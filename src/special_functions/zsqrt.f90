!** ZSQRT
SUBROUTINE ZSQRT(Ar,Ai,Br,Bi)
  !>
  !  Subsidiary to ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZAIRY and
  !            ZBIRY
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      ALL (ZSQRT-A)
  !***
  ! **Author:**  Amos, D. E., (SNL)
  !***
  ! **Description:**
  !
  !     DOUBLE PRECISION COMPLEX SQUARE ROOT, B=SQRT(A)
  !
  !***
  ! **See also:**  ZAIRY, ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZBIRY
  !***
  ! **Routines called:**  ZABS

  !* REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)

  REAL(DP) :: Ar, Ai, Br, Bi, zm, dtheta
  REAL(DP), PARAMETER :: drt = 7.071067811865475244008443621D-1, &
    dpi = 3.141592653589793238462643383D+0
  !* FIRST EXECUTABLE STATEMENT  ZSQRT
  zm = ZABS(Ar,Ai)
  zm = SQRT(zm)
  IF ( Ar==0.0D+0 ) THEN
    IF ( Ai>0.0D+0 ) THEN
      Br = zm*drt
      Bi = zm*drt
      RETURN
    ELSEIF ( Ai<0.0D+0 ) THEN
      Br = zm*drt
      Bi = -zm*drt
      RETURN
    ELSE
      Br = 0.0D+0
      Bi = 0.0D+0
      RETURN
    END IF
  ELSEIF ( Ai==0.0D+0 ) THEN
    IF ( Ar>0.0D+0 ) THEN
      Br = SQRT(Ar)
      Bi = 0.0D+0
      RETURN
    ELSE
      Br = 0.0D+0
      Bi = SQRT(ABS(Ar))
      RETURN
    END IF
  ELSE
    dtheta = ATAN(Ai/Ar)
    IF ( dtheta<=0.0D+0 ) THEN
      IF ( Ar<0.0D+0 ) dtheta = dtheta + dpi
    ELSE
      IF ( Ar<0.0D+0 ) dtheta = dtheta - dpi
    END IF
  END IF
  dtheta = dtheta*0.5D+0
  Br = zm*COS(dtheta)
  Bi = zm*SIN(dtheta)
  RETURN
END SUBROUTINE ZSQRT
