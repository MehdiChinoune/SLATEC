!*==DFEHL.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK DFEHL
SUBROUTINE DFEHL(DF,Neq,T,Y,H,Yp,F1,F2,F3,F4,F5,Ys,Rpar,Ipar)
  IMPLICIT NONE
  !*--DFEHL5
  !***BEGIN PROLOGUE  DFEHL
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to DDERKF
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (DEFEHL-S, DFEHL-D)
  !***AUTHOR  Watts, H. A., (SNLA)
  !***DESCRIPTION
  !
  !     Fehlberg Fourth-Fifth Order Runge-Kutta Method
  ! **********************************************************************
  !
  !    DFEHL integrates a system of NEQ first order
  !    ordinary differential equations of the form
  !               DU/DX = DF(X,U)
  !    over one step when the vector Y(*) of initial values for U(*) and
  !    the vector YP(*) of initial derivatives, satisfying  YP = DF(T,Y),
  !    are given at the starting point X=T.
  !
  !    DFEHL advances the solution over the fixed step H and returns
  !    the fifth order (sixth order accurate locally) solution
  !    approximation at T+H in the array YS(*).
  !    F1,---,F5 are arrays of dimension NEQ which are needed
  !    for internal storage.
  !    The formulas have been grouped to control loss of significance.
  !    DFEHL should be called with an H not smaller than 13 units of
  !    roundoff in T so that the various independent arguments can be
  !    distinguished.
  !
  !    This subroutine has been written with all variables and statement
  !    numbers entirely compatible with DRKFS. For greater efficiency,
  !    the call to DFEHL can be replaced by the module beginning with
  !    line 222 and extending to the last line just before the return
  !    statement.
  !
  ! **********************************************************************
  !
  !***SEE ALSO  DDERKF
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   820301  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910722  Updated AUTHOR section.  (ALS)
  !***END PROLOGUE  DFEHL
  !
  INTEGER Ipar , k , Neq
  DOUBLE PRECISION ch , F1 , F2 , F3 , F4 , F5 , H , Rpar , T , Y , Yp , Ys
  DIMENSION Y(*) , Yp(*) , F1(*) , F2(*) , F3(*) , F4(*) , F5(*) , Ys(*) , &
    Rpar(*) , Ipar(*)
  !
  !***FIRST EXECUTABLE STATEMENT  DFEHL
  ch = H/4.0D0
  DO k = 1 , Neq
    Ys(k) = Y(k) + ch*Yp(k)
  ENDDO
  CALL DF(T+ch,Ys,F1,Rpar,Ipar)
  !
  ch = 3.0D0*H/32.0D0
  DO k = 1 , Neq
    Ys(k) = Y(k) + ch*(Yp(k)+3.0D0*F1(k))
  ENDDO
  CALL DF(T+3.0D0*H/8.0D0,Ys,F2,Rpar,Ipar)
  !
  ch = H/2197.0D0
  DO k = 1 , Neq
    Ys(k) = Y(k) + ch*(1932.0D0*Yp(k)+(7296.0D0*F2(k)-7200.0D0*F1(k)))
  ENDDO
  CALL DF(T+12.0D0*H/13.0D0,Ys,F3,Rpar,Ipar)
  !
  ch = H/4104.0D0
  DO k = 1 , Neq
    Ys(k) = Y(k)&
      + ch*((8341.0D0*Yp(k)-845.0D0*F3(k))+(29440.0D0*F2(k)-32832.0D0*&
      F1(k)))
  ENDDO
  CALL DF(T+H,Ys,F4,Rpar,Ipar)
  !
  ch = H/20520.0D0
  DO k = 1 , Neq
    Ys(k) = Y(k) + ch*((-6080.0D0*Yp(k)+(9295.0D0*F3(k)-5643.0D0*F4(k)))&
      +(41040.0D0*F1(k)-28352.0D0*F2(k)))
  ENDDO
  CALL DF(T+H/2.0D0,Ys,F5,Rpar,Ipar)
  !
  !     COMPUTE APPROXIMATE SOLUTION AT T+H
  !
  ch = H/7618050.0D0
  DO k = 1 , Neq
    Ys(k) = Y(k)&
      + ch*((902880.0D0*Yp(k)+(3855735.0D0*F3(k)-1371249.0D0*F4(k)))&
      +(3953664.0D0*F2(k)+277020.0D0*F5(k)))
  ENDDO
  !
END SUBROUTINE DFEHL
