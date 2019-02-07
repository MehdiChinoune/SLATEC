!*==JAC.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK JAC
SUBROUTINE JAC(T,U,Pd,Nrowpd,Rpar,Ipar)
  IMPLICIT NONE
  !*--JAC5
  !***BEGIN PROLOGUE  JAC
  !***SUBSIDIARY
  !***PURPOSE  Evaluate Jacobian for DEBDF quick check.
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (JAC-S, DJAC-D)
  !***AUTHOR  Chow, Jeff (LANL)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   810801  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900415  Minor clean-up of prologue and code.  (WRB)
  !***END PROLOGUE  JAC
  INTEGER Ipar , Nrowpd
  REAL Pd , r , r5 , Rpar , rsq , T , U , u1sq , u2sq , u1u2
  DIMENSION U(*) , Pd(Nrowpd,*) , Rpar(*) , Ipar(*)
  !***FIRST EXECUTABLE STATEMENT  JAC
  u1sq = U(1)*U(1)
  u2sq = U(2)*U(2)
  u1u2 = U(1)*U(2)
  rsq = u1sq + u2sq
  r = SQRT(rsq)
  r5 = rsq*rsq*r
  Pd(3,1) = (3.E0*u1sq-rsq)/r5
  Pd(4,1) = 3.E0*u1u2/r5
  Pd(3,2) = Pd(4,1)
  Pd(4,2) = (3.E0*u2sq-rsq)/r5
  Pd(1,3) = 1.E0
  Pd(2,4) = 1.E0
END SUBROUTINE JAC
