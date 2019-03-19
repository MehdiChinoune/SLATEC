!** CSHCH
SUBROUTINE CSHCH(Z,Csh,Cch)
  IMPLICIT NONE
  !>
  !***
  !  Subsidiary to CBESH and CBESK
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      ALL (CSHCH-A, ZSHCH-A)
  !***
  ! **Author:**  Amos, D. E., (SNL)
  !***
  ! **Description:**
  !
  !     CSHCH COMPUTES THE COMPLEX HYPERBOLIC FUNCTIONS CSH=SINH(X+I*Y)
  !     AND CCH=COSH(X+I*Y), WHERE I**2=-1.
  !
  !***
  ! **See also:**  CBESH, CBESK
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)
  
  COMPLEX Cch, Csh, Z
  REAL cchi, cchr, ch, cn, cshi, cshr, sh, sn, x, y
  !* FIRST EXECUTABLE STATEMENT  CSHCH
  x = REAL(Z)
  y = AIMAG(Z)
  sh = SINH(x)
  ch = COSH(x)
  sn = SIN(y)
  cn = COS(y)
  cshr = sh*cn
  cshi = ch*sn
  Csh = CMPLX(cshr,cshi)
  cchr = ch*cn
  cchi = sh*sn
  Cch = CMPLX(cchr,cchi)
END SUBROUTINE CSHCH
