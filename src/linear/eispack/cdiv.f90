!** CDIV
SUBROUTINE CDIV(Ar,Ai,Br,Bi,Cr,Ci)
  !>
  !  Compute the complex quotient of two complex numbers.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      COMPLEX (CDIV-C)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !     Complex division, (CR,CI) = (AR,AI)/(BR,BI)
  !
  !***
  ! **See also:**  EISDOC
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   811101  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)

  REAL(SP) Ar, Ai, Br, Bi, Cr, Ci
  !
  REAL(SP) s, ars, ais, brs, bis
  !* FIRST EXECUTABLE STATEMENT  CDIV
  s = ABS(Br) + ABS(Bi)
  ars = Ar/s
  ais = Ai/s
  brs = Br/s
  bis = Bi/s
  s = brs**2 + bis**2
  Cr = (ars*brs+ais*bis)/s
  Ci = (ais*brs-ars*bis)/s
END SUBROUTINE CDIV
