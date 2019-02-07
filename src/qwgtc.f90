!*==QWGTC.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK QWGTC
REAL FUNCTION QWGTC(X,C,P2,P3,P4,Kp)
  IMPLICIT NONE
  !*--QWGTC5
  !***BEGIN PROLOGUE  QWGTC
  !***SUBSIDIARY
  !***PURPOSE  This function subprogram is used together with the
  !            routine QAWC and defines the WEIGHT function.
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (QWGTC-S, DQWGTC-D)
  !***KEYWORDS  CAUCHY PRINCIPAL VALUE, WEIGHT FUNCTION
  !***AUTHOR  Piessens, Robert
  !             Applied Mathematics and Programming Division
  !             K. U. Leuven
  !           de Doncker, Elise
  !             Applied Mathematics and Programming Division
  !             K. U. Leuven
  !***SEE ALSO  QK15W
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   810101  DATE WRITTEN
  !   830518  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !***END PROLOGUE  QWGTC
  !
  REAL C , P2 , P3 , P4 , X
  INTEGER Kp
  !***FIRST EXECUTABLE STATEMENT  QWGTC
  QWGTC = 0.1E+01/(X-C)
END FUNCTION QWGTC
