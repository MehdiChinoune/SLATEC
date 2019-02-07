!*==ZABS.f90  processed by SPAG 6.72Dc at 11:02 on  6 Feb 2019
!DECK ZABS
DOUBLE PRECISION FUNCTION ZABS(Zr,Zi)
  IMPLICIT NONE
  !*--ZABS5
  !***BEGIN PROLOGUE  ZABS
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZAIRY and
  !            ZBIRY
  !***LIBRARY   SLATEC
  !***TYPE      ALL (ZABS-A)
  !***AUTHOR  Amos, D. E., (SNL)
  !***DESCRIPTION
  !
  !     ZABS COMPUTES THE ABSOLUTE VALUE OR MAGNITUDE OF A DOUBLE
  !     PRECISION COMPLEX VARIABLE CMPLX(ZR,ZI)
  !
  !***SEE ALSO  ZAIRY, ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZBIRY
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  ZABS
  DOUBLE PRECISION Zr , Zi , u , v , q , s
  !***FIRST EXECUTABLE STATEMENT  ZABS
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
    GOTO 99999
  ELSEIF ( u<=v ) THEN
    q = u/v
    ZABS = v*SQRT(1.D+0+q*q)
    RETURN
  ENDIF
  q = v/u
  ZABS = u*SQRT(1.D+0+q*q)
  RETURN
  99999 END FUNCTION ZABS
