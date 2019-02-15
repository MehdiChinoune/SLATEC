!DECK CSROOT
SUBROUTINE CSROOT(Xr,Xi,Yr,Yi)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  CSROOT
  !***SUBSIDIARY
  !***PURPOSE  Compute the complex square root of a complex number.
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (CSROOT-S)
  !***AUTHOR  (UNKNOWN)
  !***DESCRIPTION
  !
  !     (YR,YI) = complex sqrt(XR,XI)
  !
  !***SEE ALSO  EISDOC
  !***ROUTINES CALLED  PYTHAG
  !***REVISION HISTORY  (YYMMDD)
  !   811101  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !***END PROLOGUE  CSROOT
  REAL Xr, Xi, Yr, Yi, s, tr, ti, PYTHAG
  !
  !     BRANCH CHOSEN SO THAT YR .GE. 0.0 AND SIGN(YI) .EQ. SIGN(XI)
  !***FIRST EXECUTABLE STATEMENT  CSROOT
  tr = Xr
  ti = Xi
  s = SQRT(0.5E0*(PYTHAG(tr,ti)+ABS(tr)))
  IF ( tr>=0.0E0 ) Yr = s
  IF ( ti<0.0E0 ) s = -s
  IF ( tr<=0.0E0 ) Yi = s
  IF ( tr<0.0E0 ) Yr = 0.5E0*(ti/Yi)
  IF ( tr>0.0E0 ) Yi = 0.5E0*(ti/Yr)
END SUBROUTINE CSROOT
