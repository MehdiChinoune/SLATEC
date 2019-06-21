!** ZLOG
SUBROUTINE ZLOG(Ar,Ai,Br,Bi,Ierr)
  !> Subsidiary to ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZAIRY and
  !            ZBIRY
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      ALL (ZLOG-A)
  !***
  ! **Author:**  Amos, D. E., (SNL)
  !***
  ! **Description:**
  !
  !     DOUBLE PRECISION COMPLEX LOGARITHM B=LOG(A)
  !     IERR=0,NORMAL RETURN      IERR=1, Z=CMPLX(0.0,0.0)
  !***
  ! **See also:**  ZAIRY, ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZBIRY
  !***
  ! **Routines called:**  ZABS

  !* REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)

  REAL(DP) :: Ar, Ai, Br, Bi, zm, dtheta
  INTEGER :: Ierr
  REAL(DP), PARAMETER :: dpi = 3.141592653589793238462643383E+0_DP, &
    dhpi = 1.570796326794896619231321696_DP
  !* FIRST EXECUTABLE STATEMENT  ZLOG
  Ierr = 0
  IF( Ar==0._DP ) THEN
    IF( Ai==0._DP ) THEN
      Ierr = 1
      RETURN
    ELSE
      Bi = dhpi
      Br = LOG(ABS(Ai))
      IF( Ai<0._DP ) Bi = -Bi
      RETURN
    END IF
  ELSEIF( Ai==0._DP ) THEN
    IF( Ar>0._DP ) THEN
      Br = LOG(Ar)
      Bi = 0._DP
      RETURN
    ELSE
      Br = LOG(ABS(Ar))
      Bi = dpi
      RETURN
    END IF
  ELSE
    dtheta = ATAN(Ai/Ar)
    IF( dtheta<=0._DP ) THEN
      IF( Ar<0._DP ) dtheta = dtheta + dpi
    ELSE
      IF( Ar<0._DP ) dtheta = dtheta - dpi
    END IF
  END IF
  zm = ZABS(Ar,Ai)
  Br = LOG(zm)
  Bi = dtheta
  RETURN
END SUBROUTINE ZLOG
