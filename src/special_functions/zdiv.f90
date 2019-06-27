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

  COMPLEX(DP) :: C
  REAL(DP) :: Ar, Ai, Br, Bi, Cr, Ci
  !* FIRST EXECUTABLE STATEMENT  ZDIV
  C = CMPLX(Ar,Ai,DP)/CMPLX(Br,Bi,DP)
  Cr = REAL(C,DP)
  Ci = AIMAG(C)

END SUBROUTINE ZDIV