!** ZSHCH
ELEMENTAL SUBROUTINE ZSHCH(Z,Csh,Cch)
  !> Subsidiary to ZBESH and ZBESK
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

  COMPLEX(DP), INTENT(IN) :: Z
  COMPLEX(DP), INTENT(OUT) :: Cch, Csh
  !* FIRST EXECUTABLE STATEMENT  ZSHCH
  Csh = SINH( Z )
  Cch = COSH( Z )
  !
END SUBROUTINE ZSHCH