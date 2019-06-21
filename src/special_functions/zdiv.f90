!** ZDIV
SUBROUTINE ZDIV(Ar,Ai,Br,Bi,Cr,Ci)
  !> Subsidiary to ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZAIRY and
  !            ZBIRY
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      ALL (ZDIV-A)
  !***
  ! **Author:**  Amos, D. E., (SNL)
  !***
  ! **Description:**
  !
  !     DOUBLE PRECISION COMPLEX DIVIDE C=A/B.
  !
  !***
  ! **See also:**  ZAIRY, ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZBIRY
  !***
  ! **Routines called:**  ZABS

  !* REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)

  REAL(DP) :: Ar, Ai, Br, Bi, Cr, Ci, bm, ca, cb, cc, cd
  !* FIRST EXECUTABLE STATEMENT  ZDIV
  bm = 1._DP/ZABS(Br,Bi)
  cc = Br*bm
  cd = Bi*bm
  ca = (Ar*cc+Ai*cd)*bm
  cb = (Ai*cc-Ar*cd)*bm
  Cr = ca
  Ci = cb
END SUBROUTINE ZDIV
