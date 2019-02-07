!*==UVEC.f90  processed by SPAG 6.72Dc at 11:02 on  6 Feb 2019
!DECK UVEC
SUBROUTINE UVEC(X,Y,Yp)
  IMPLICIT NONE
  !*--UVEC5
  !***BEGIN PROLOGUE  UVEC
  !***PURPOSE  Dummy routine for BVSUP quick check.
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (UVEC-S, DUVEC-D)
  !***KEYWORDS  QUICK CHECK
  !***AUTHOR  Watts, H. A., (SNLA)
  !***DESCRIPTION
  !
  !   This routine is never called;  it is here to prevent loaders from
  !   complaining about undefined externals while testing BVSUP.
  !
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   750601  DATE WRITTEN
  !   890618  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920401  Variables declaration and TYPE sections added.  (WRB)
  !***END PROLOGUE  UVEC
  !     .. Scalar Arguments ..
  REAL X
  !     .. Array Arguments ..
  REAL Y(*), Yp(*)
  !***FIRST EXECUTABLE STATEMENT  UVEC
  STOP
END SUBROUTINE UVEC
