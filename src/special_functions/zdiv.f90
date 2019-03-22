!** ZDIV
SUBROUTINE ZDIV(Ar,Ai,Br,Bi,Cr,Ci)
  IMPLICIT NONE
  !>
  !***
  !  Subsidiary to ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZAIRY and
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
  
  REAL(8) :: Ar, Ai, Br, Bi, Cr, Ci, bm, ca, cb, cc, cd
  REAL(8) :: ZABS
  EXTERNAL :: ZABS
  !* FIRST EXECUTABLE STATEMENT  ZDIV
  bm = 1.0D0/ZABS(Br,Bi)
  cc = Br*bm
  cd = Bi*bm
  ca = (Ar*cc+Ai*cd)*bm
  cb = (Ai*cc-Ar*cd)*bm
  Cr = ca
  Ci = cb
END SUBROUTINE ZDIV
