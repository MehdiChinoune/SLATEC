!*==DEFEHL.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK DEFEHL
SUBROUTINE DEFEHL(F,Neq,T,Y,H,Yp,F1,F2,F3,F4,F5,Ys,Rpar,Ipar)
  IMPLICIT NONE
  !*--DEFEHL5
  !*** Start of declarations inserted by SPAG
  REAL ch, F1, F2, F3, F4, F5, H, Rpar, T, Y, Yp, Ys
  INTEGER Ipar, k, Neq
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  DEFEHL
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to DERKF
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (DEFEHL-S, DFEHL-D)
  !***AUTHOR  Watts, H. A., (SNLA)
  !***DESCRIPTION
  !
  !     Fehlberg Fourth-Fifth order Runge-Kutta Method
  ! **********************************************************************
  !
  !    DEFEHL integrates a system of NEQ first order
  !    ordinary differential equations of the form
  !               dU/DX = F(X,U)
  !    over one step when the vector Y(*) of initial values for U(*) and
  !    the vector YP(*) of initial derivatives, satisfying  YP = F(T,Y),
  !    are given at the starting point X=T.
  !
  !    DEFEHL advances the solution over the fixed step H and returns
  !    the fifth order (sixth order accurate locally) solution
  !    approximation at T+H in the array YS(*).
  !    F1,---,F5 are arrays of dimension NEQ which are needed
  !    for internal storage.
  !    The formulas have been grouped to control loss of significance.
  !    DEFEHL should be called with an H not smaller than 13 units of
  !    roundoff in T so that the various independent arguments can be
  !    distinguished.
  !
  !    This subroutine has been written with all variables and statement
  !    numbers entirely compatible with DERKFS. For greater efficiency,
  !    the call to DEFEHL can be replaced by the module beginning with
  !    line 222 and extending to the last line just before the return
  !    statement.
  !
  ! **********************************************************************
  !
  !***SEE ALSO  DERKF
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   800501  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   891009  Removed unreferenced statement label.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910722  Updated AUTHOR section.  (ALS)
  !***END PROLOGUE  DEFEHL
  !
  !
  DIMENSION Y(*), Yp(*), F1(*), F2(*), F3(*), F4(*), F5(*), Ys(*), &
    Rpar(*), Ipar(*)
  !
  !***FIRST EXECUTABLE STATEMENT  DEFEHL
  ch = H/4.
  DO k = 1, Neq
    Ys(k) = Y(k) + ch*Yp(k)
  ENDDO
  CALL F(T+ch,Ys,F1,Rpar,Ipar)
  !
  ch = 3.*H/32.
  DO k = 1, Neq
    Ys(k) = Y(k) + ch*(Yp(k)+3.*F1(k))
  ENDDO
  CALL F(T+3.*H/8.,Ys,F2,Rpar,Ipar)
  !
  ch = H/2197.
  DO k = 1, Neq
    Ys(k) = Y(k) + ch*(1932.*Yp(k)+(7296.*F2(k)-7200.*F1(k)))
  ENDDO
  CALL F(T+12.*H/13.,Ys,F3,Rpar,Ipar)
  !
  ch = H/4104.
  DO k = 1, Neq
    Ys(k) = Y(k) + ch*((8341.*Yp(k)-845.*F3(k))+(29440.*F2(k)-32832.*F1(k)))
  ENDDO
  CALL F(T+H,Ys,F4,Rpar,Ipar)
  !
  ch = H/20520.
  DO k = 1, Neq
    Ys(k) = Y(k) + ch*((-6080.*Yp(k)+(9295.*F3(k)-5643.*F4(k)))&
      +(41040.*F1(k)-28352.*F2(k)))
  ENDDO
  CALL F(T+H/2.,Ys,F5,Rpar,Ipar)
  !
  !     COMPUTE APPROXIMATE SOLUTION AT T+H
  !
  ch = H/7618050.
  DO k = 1, Neq
    Ys(k) = Y(k) + ch*((902880.*Yp(k)+(3855735.*F3(k)-1371249.*F4(k)))&
      +(3953664.*F2(k)+277020.*F5(k)))
  ENDDO
  !
END SUBROUTINE DEFEHL
