!** FMAT
SUBROUTINE FMAT(X,Y,Yp)
  !>
  !  Subsidiary to
  !***
  ! **Library:**   SLATEC
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Routines called:**  (NONE)
  !***
  ! COMMON BLOCKS    SAVEX

  !* REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  USE SAVEX, ONLY : XSAve, TERm
  REAL :: X, Y(*), Yp(*)
  REAL :: tanx
  !* FIRST EXECUTABLE STATEMENT  FMAT
  Yp(1) = Y(2)
  IF ( X/=XSAve ) THEN
    XSAve = X
    tanx = TAN(X/57.2957795130823)
    TERm = 3.0/tanx + 2.0*tanx
  END IF
  Yp(2) = -TERm*Y(2) - 0.7*Y(1)
END SUBROUTINE FMAT
