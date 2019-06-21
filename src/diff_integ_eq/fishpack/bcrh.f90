!** BCRH
REAL(SP) FUNCTION BCRH(Xll,Xrr,Iz,C,A,Bh,F,Sgn)
  !> Subsidiary to CBLKTR
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
    REAL(SP) FUNCTION F(X,Iz,C,A,Bh)
      IMPORT SP
      INTEGER :: Iz
      REAL(SP) :: X, A(Iz), Bh(Iz), C(Iz)
    END FUNCTION F
  END INTERFACE
  INTEGER :: Iz
  REAL(SP) :: A(:), Bh(:), C(:)
  REAL(SP) :: Sgn, Xll, Xrr
  REAL(SP) :: dx, x, xl, xr
  !* FIRST EXECUTABLE STATEMENT  BCRH
  xl = Xll
  xr = Xrr
  dx = 0.5_SP*ABS(xr-xl)
  100  x = 0.5_SP*(xl+xr)
  IF( Sgn*F(x,Iz,C,A,Bh)<0 ) THEN
    xl = x
  ELSEIF( Sgn*F(x,Iz,C,A,Bh)==0 ) THEN
    BCRH = 0.5_SP*(xl+xr)
    RETURN
  ELSE
    xr = x
  END IF
  dx = 0.5_SP*dx
  IF( dx>cnv_com ) GOTO 100
  BCRH = 0.5_SP*(xl+xr)
  RETURN
END FUNCTION BCRH
