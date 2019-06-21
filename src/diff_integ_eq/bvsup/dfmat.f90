!** DFMAT
SUBROUTINE DFMAT(X,Y,Yp)
  !> Subsidiary to
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
  USE DSAVEX, ONLY : xsave_com, term_com
  REAL(DP) :: X, Y(:), Yp(:)
  REAL(DP) :: tanx
  !* FIRST EXECUTABLE STATEMENT  DFMAT
  Yp(1) = Y(2)
  IF( X/=xsave_com ) THEN
    xsave_com = X
    tanx = TAN(X/57.2957795130823_DP)
    term_com = 3._DP/tanx + 2._DP*tanx
  END IF
  Yp(2) = -term_com*Y(2) - 0.7_DP*Y(1)
END SUBROUTINE DFMAT
