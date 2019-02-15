!DECK MPMLP
SUBROUTINE MPMLP(U,V,W,J)
  IMPLICIT NONE
  INTEGER i, J
  !***BEGIN PROLOGUE  MPMLP
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to DQDOTA and DQDOTI
  !***LIBRARY   SLATEC
  !***TYPE      ALL (MPMLP-A)
  !***AUTHOR  (UNKNOWN)
  !***DESCRIPTION
  !
  ! Performs inner multiplication loop for MPMUL. Carries are not pro-
  ! pagated in inner loop, which saves time at the expense of space.
  !
  !***SEE ALSO  DQDOTA, DQDOTI
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   791001  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !***END PROLOGUE  MPMLP
  INTEGER U(*), V(*), W
  !***FIRST EXECUTABLE STATEMENT  MPMLP
  DO i = 1, J
    U(i) = U(i) + W*V(i)
  ENDDO
END SUBROUTINE MPMLP
