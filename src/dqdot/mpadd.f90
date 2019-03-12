!DECK MPADD
SUBROUTINE MPADD(X,Y,Z)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  MPADD
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to DQDOTA and DQDOTI
  !***LIBRARY   SLATEC
  !***TYPE      ALL (MPADD-A)
  !***AUTHOR  (UNKNOWN)
  !***DESCRIPTION
  !
  ! Adds X and Y, forming result in Z, where X, Y and Z are 'mp'
  !  (multiple precision) numbers.  Four guard digits are used,
  !  and then R*-rounding.
  !
  !***SEE ALSO  DQDOTA, DQDOTI
  !***ROUTINES CALLED  MPADD2
  !***REVISION HISTORY  (YYMMDD)
  !   791001  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !***END PROLOGUE  MPADD
  INTEGER X(*), Y(*), Z(*)
  !***FIRST EXECUTABLE STATEMENT  MPADD
  CALL MPADD2(X,Y,Z,Y,0)
END SUBROUTINE MPADD
