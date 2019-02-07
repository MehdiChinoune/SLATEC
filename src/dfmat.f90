!*==DFMAT.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK DFMAT
SUBROUTINE DFMAT(X,Y,Yp)
  IMPLICIT NONE
  !*--DFMAT5
  !***BEGIN PROLOGUE  DFMAT
  !***PURPOSE  Subsidiary to
  !***LIBRARY   SLATEC
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  (NONE)
  !***COMMON BLOCKS    DSAVEX
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  DFMAT
  REAL(8) :: X, Y, Yp, XSAve, TERm, tanx
  DIMENSION Y(*), Yp(*)
  COMMON /DSAVEX/ XSAve, TERm
  !***FIRST EXECUTABLE STATEMENT  DFMAT
  Yp(1) = Y(2)
  IF ( X/=XSAve ) THEN
    XSAve = X
    tanx = TAN(X/57.2957795130823D0)
    TERm = 3.0D0/tanx + 2.0D0*tanx
  ENDIF
  Yp(2) = -TERm*Y(2) - 0.7D0*Y(1)
END SUBROUTINE DFMAT
