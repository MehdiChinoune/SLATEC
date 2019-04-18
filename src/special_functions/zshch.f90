!** ZSHCH
SUBROUTINE ZSHCH(Zr,Zi,Cshr,Cshi,Cchr,Cchi)
  !>
  !  Subsidiary to ZBESH and ZBESK
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      ALL (CSHCH-A, ZSHCH-A)
  !***
  ! **Author:**  Amos, D. E., (SNL)
  !***
  ! **Description:**
  !
  !     ZSHCH COMPUTES THE COMPLEX HYPERBOLIC FUNCTIONS CSH=SINH(X+I*Y)
  !     AND CCH=COSH(X+I*Y), WHERE I**2=-1.
  !
  !***
  ! **See also:**  ZBESH, ZBESK
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)
  
  !
  REAL(8) :: Cchi, Cchr, ch, cn, Cshi, Cshr, sh, sn, Zi, Zr
  !* FIRST EXECUTABLE STATEMENT  ZSHCH
  sh = SINH(Zr)
  ch = COSH(Zr)
  sn = SIN(Zi)
  cn = COS(Zi)
  Cshr = sh*cn
  Cshi = ch*sn
  Cchr = ch*cn
  Cchi = sh*sn
END SUBROUTINE ZSHCH
