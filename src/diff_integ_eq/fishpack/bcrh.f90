!** BCRH
REAL FUNCTION BCRH(Xll,Xrr,Iz,C,A,Bh,F,Sgn)
  !>
  !  Subsidiary to CBLKTR
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (BCRH-S, BSRH-S)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **See also:**  CBLKTR
  !***
  ! **Routines called:**  (NONE)
  !***
  ! COMMON BLOCKS    CCBLK

  !* REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  USE CCBLK, ONLY : cnv_com
  INTERFACE
    REAL FUNCTION F(X,Iz,C,A,Bh)
      INTEGER :: Iz
      REAL :: X, A(Iz), Bh(Iz), C(Iz)
    END FUNCTION F
  END INTERFACE
  INTEGER :: Iz
  REAL :: A(:), Bh(:), C(:)
  REAL :: Sgn, Xll, Xrr
  REAL :: dx, x, xl, xr
  !* FIRST EXECUTABLE STATEMENT  BCRH
  xl = Xll
  xr = Xrr
  dx = .5*ABS(xr-xl)
  100  x = .5*(xl+xr)
  IF ( Sgn*F(x,Iz,C,A,Bh)<0 ) THEN
    xl = x
  ELSEIF ( Sgn*F(x,Iz,C,A,Bh)==0 ) THEN
    BCRH = .5*(xl+xr)
    RETURN
  ELSE
    xr = x
  END IF
  dx = .5*dx
  IF ( dx>cnv_com ) GOTO 100
  BCRH = .5*(xl+xr)
  RETURN
END FUNCTION BCRH
