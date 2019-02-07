!*==ZSQRT.f90  processed by SPAG 6.72Dc at 11:02 on  6 Feb 2019
!DECK ZSQRT
SUBROUTINE ZSQRT(Ar,Ai,Br,Bi)
  IMPLICIT NONE
  !*--ZSQRT5
  !***BEGIN PROLOGUE  ZSQRT
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZAIRY and
  !            ZBIRY
  !***LIBRARY   SLATEC
  !***TYPE      ALL (ZSQRT-A)
  !***AUTHOR  Amos, D. E., (SNL)
  !***DESCRIPTION
  !
  !     DOUBLE PRECISION COMPLEX SQUARE ROOT, B=CSQRT(A)
  !
  !***SEE ALSO  ZAIRY, ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZBIRY
  !***ROUTINES CALLED  ZABS
  !***REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  ZSQRT
  DOUBLE PRECISION Ar , Ai , Br , Bi , zm , dtheta , dpi , drt
  DOUBLE PRECISION ZABS
  EXTERNAL ZABS
  DATA drt , dpi/7.071067811865475244008443621D-1 , &
    3.141592653589793238462643383D+0/
  !***FIRST EXECUTABLE STATEMENT  ZSQRT
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
      GOTO 99999
    ELSE
      Br = 0.0D+0
      Bi = 0.0D+0
      RETURN
    ENDIF
  ELSEIF ( Ai==0.0D+0 ) THEN
    IF ( Ar>0.0D+0 ) THEN
      Br = SQRT(Ar)
      Bi = 0.0D+0
      RETURN
    ELSE
      Br = 0.0D+0
      Bi = SQRT(ABS(Ar))
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
  dtheta = dtheta*0.5D+0
  Br = zm*COS(dtheta)
  Bi = zm*SIN(dtheta)
  RETURN
  99999 END SUBROUTINE ZSQRT
