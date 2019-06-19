!** ZMLT
SUBROUTINE ZMLT(Ar,Ai,Br,Bi,Cr,Ci)
  !> Subsidiary to ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZAIRY and
  !            ZBIRY
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      ALL (ZMLT-A)
  !***
  ! **Author:**  Amos, D. E., (SNL)
  !***
  ! **Description:**
  !
  !     DOUBLE PRECISION COMPLEX MULTIPLY, C=A*B.
  !
  !***
  ! **See also:**  ZAIRY, ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZBIRY
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)
  
  REAL(DP) :: Ar, Ai, Br, Bi, Cr, Ci, ca, cb
  !* FIRST EXECUTABLE STATEMENT  ZMLT
  ca = Ar*Br - Ai*Bi
  cb = Ar*Bi + Ai*Br
  Cr = ca
  Ci = cb
END SUBROUTINE ZMLT
