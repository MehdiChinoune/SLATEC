!*==QWGTF.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK QWGTF
      REAL FUNCTION QWGTF(X,Omega,P2,P3,P4,Integr)
      IMPLICIT NONE
!*--QWGTF5
!***BEGIN PROLOGUE  QWGTF
!***SUBSIDIARY
!***PURPOSE  This function subprogram is used together with the
!            routine QAWF and defines the WEIGHT function.
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (QWGTF-S, DQWGTF-D)
!***KEYWORDS  COS OR SIN IN WEIGHT FUNCTION
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
!***END PROLOGUE  QWGTF
!
      REAL Omega , omx , P2 , P3 , P4 , X
      INTEGER Integr
!***FIRST EXECUTABLE STATEMENT  QWGTF
      omx = Omega*X
      IF ( Integr==2 ) THEN
        QWGTF = SIN(omx)
      ELSE
        QWGTF = COS(omx)
      ENDIF
      END FUNCTION QWGTF
