!** CSHCH
ELEMENTAL SUBROUTINE CSHCH(Z,Csh,Cch)
  !> Subsidiary to CBESH and CBESK
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

  COMPLEX(SP), INTENT(IN) :: Z
  COMPLEX(SP), INTENT(OUT) :: Cch, Csh
  !* FIRST EXECUTABLE STATEMENT  CSHCH
  Csh = SINH( Z )
  Cch = COSH( Z )

END SUBROUTINE CSHCH