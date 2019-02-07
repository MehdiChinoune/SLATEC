!*==DQWGTF.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK DQWGTF
REAL(8) FUNCTION DQWGTF(X,Omega,P2,P3,P4,Integr)
  IMPLICIT NONE
  !*--DQWGTF5
  !***BEGIN PROLOGUE  DQWGTF
  !***SUBSIDIARY
  !***PURPOSE  This function subprogram is used together with the
  !            routine DQAWF and defines the WEIGHT function.
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (QWGTF-S, DQWGTF-D)
  !***KEYWORDS  COS OR SIN IN WEIGHT FUNCTION
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
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !***END PROLOGUE  DQWGTF
  !
  REAL(8) :: Omega, omx, P2, P3, P4, X
  INTEGER Integr
  !***FIRST EXECUTABLE STATEMENT  DQWGTF
  omx = Omega*X
  IF ( Integr==2 ) THEN
    DQWGTF = SIN(omx)
  ELSE
    DQWGTF = COS(omx)
  ENDIF
END FUNCTION DQWGTF
