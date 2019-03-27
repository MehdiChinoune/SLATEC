!** DFMAT
SUBROUTINE DFMAT(X,Y,Yp)
  IMPLICIT NONE
  !>
  !***
  !  Subsidiary to
  !***
  ! **Library:**   SLATEC
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Routines called:**  (NONE)
  !***
  ! COMMON BLOCKS    DSAVEX

  !* REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  
  REAL(8) :: X, Y(*), Yp(*), XSAve, TERm, tanx
  COMMON /DSAVEX/ XSAve, TERm
  !* FIRST EXECUTABLE STATEMENT  DFMAT
  Yp(1) = Y(2)
  IF ( X/=XSAve ) THEN
    XSAve = X
    tanx = TAN(X/57.2957795130823D0)
    TERm = 3.0D0/tanx + 2.0D0*tanx
  ENDIF
  Yp(2) = -TERm*Y(2) - 0.7D0*Y(1)
END SUBROUTINE DFMAT
