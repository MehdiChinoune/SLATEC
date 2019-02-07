!*==INDXC.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK INDXC
SUBROUTINE INDXC(I,Ir,Idxc,Nc)
  IMPLICIT NONE
  !*--INDXC5
  !*** Start of declarations inserted by SPAG
  REAL CNV, EPS
  INTEGER I, Idxc, IK, Ir, K, Nc, NCMplx, NM, NPP
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  INDXC
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to BLKTRI
  !***LIBRARY   SLATEC
  !***TYPE      INTEGER (INDXC-I)
  !***AUTHOR  (UNKNOWN)
  !***SEE ALSO  BLKTRI
  !***ROUTINES CALLED  (NONE)
  !***COMMON BLOCKS    CBLKT
  !***REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !***END PROLOGUE  INDXC
  COMMON /CBLKT / NPP, K, EPS, CNV, NM, NCMplx, IK
  !***FIRST EXECUTABLE STATEMENT  INDXC
  Nc = 2**Ir
  Idxc = I
  IF ( Idxc+Nc-1>NM ) Nc = 0
END SUBROUTINE INDXC
