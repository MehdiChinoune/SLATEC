!*==DQWGTC.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK DQWGTC
REAL(8) FUNCTION DQWGTC(X,C,P2,P3,P4,Kp)
  IMPLICIT NONE
  !*--DQWGTC5
  !***BEGIN PROLOGUE  DQWGTC
  !***SUBSIDIARY
  !***PURPOSE  This function subprogram is used together with the
  !            routine DQAWC and defines the WEIGHT function.
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (QWGTC-S, DQWGTC-D)
  !***KEYWORDS  CAUCHY PRINCIPAL VALUE, WEIGHT FUNCTION
  !***AUTHOR  Piessens, Robert
  !             Applied Mathematics and Programming Division
  !             K. U. Leuven
  !           de Doncker, Elise
  !             Applied Mathematics and Programming Division
  !             K. U. Leuven
  !***SEE ALSO  DQK15W
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   810101  DATE WRITTEN
  !   830518  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !***END PROLOGUE  DQWGTC
  !
  REAL(8) :: C , P2 , P3 , P4 , X
  INTEGER Kp
  !***FIRST EXECUTABLE STATEMENT  DQWGTC
  DQWGTC = 0.1D+01/(X-C)
END FUNCTION DQWGTC
