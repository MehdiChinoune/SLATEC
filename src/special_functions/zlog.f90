!** ZLOG
SUBROUTINE ZLOG(Ar,Ai,Br,Bi,Ierr)
  IMPLICIT NONE
  !>
  !***
  !  Subsidiary to ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZAIRY and
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
  !     DOUBLE PRECISION COMPLEX LOGARITHM B=CLOG(A)
  !     IERR=0,NORMAL RETURN      IERR=1, Z=CMPLX(0.0,0.0)
  !***
  ! **See also:**  ZAIRY, ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZBIRY
  !***
  ! **Routines called:**  ZABS

  !* REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)
  
  REAL(8) :: Ar, Ai, Br, Bi, zm, dtheta, dpi, dhpi
  REAL(8) :: ZABS
  INTEGER Ierr
  EXTERNAL :: ZABS
  DATA dpi, dhpi/3.141592653589793238462643383D+0, &
    1.570796326794896619231321696D+0/
  !* FIRST EXECUTABLE STATEMENT  ZLOG
  Ierr = 0
  IF ( Ar==0.0D+0 ) THEN
    IF ( Ai==0.0D+0 ) THEN
      Ierr = 1
      RETURN
    ELSE
      Bi = dhpi
      Br = LOG(ABS(Ai))
      IF ( Ai<0.0D+0 ) Bi = -Bi
      RETURN
    ENDIF
  ELSEIF ( Ai==0.0D+0 ) THEN
    IF ( Ar>0.0D+0 ) THEN
      Br = LOG(Ar)
      Bi = 0.0D+0
      RETURN
    ELSE
      Br = LOG(ABS(Ar))
      Bi = dpi
      RETURN
    ENDIF
  ELSE
    dtheta = DATAN(Ai/Ar)
    IF ( dtheta<=0.0D+0 ) THEN
      IF ( Ar<0.0D+0 ) dtheta = dtheta + dpi
    ELSE
      IF ( Ar<0.0D+0 ) dtheta = dtheta - dpi
    ENDIF
  ENDIF
  zm = ZABS(Ar,Ai)
  Br = LOG(zm)
  Bi = dtheta
  RETURN
END SUBROUTINE ZLOG
