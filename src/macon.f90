!*==MACON.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK MACON
SUBROUTINE MACON
  IMPLICIT NONE
  !*--MACON5
  !*** Start of declarations inserted by SPAG
  REAL dd, EPS, FOUru, R1MACH, SQOvfl, SRU, TWOu, URO
  INTEGER ke, LPAr
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  MACON
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to BVSUP
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (MACON-S, DMACON-D)
  !***AUTHOR  (UNKNOWN)
  !***DESCRIPTION
  !
  !    Sets up machine constants using R1MACH
  !
  !***SEE ALSO  BVSUP
  !***ROUTINES CALLED  R1MACH
  !***COMMON BLOCKS    ML5MCO
  !***REVISION HISTORY  (YYMMDD)
  !   750601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890921  Realigned order of variables in certain COMMON blocks.
  !           (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !***END PROLOGUE  MACON
  COMMON /ML5MCO/ URO, SRU, EPS, SQOvfl, TWOu, FOUru, LPAr
  !***FIRST EXECUTABLE STATEMENT  MACON
  URO = R1MACH(4)
  SRU = SQRT(URO)
  dd = -LOG10(URO)
  LPAr = 0.5*dd
  ke = 0.5 + 0.75*dd
  EPS = 10.**(-2*ke)
  SQOvfl = SQRT(R1MACH(2))
  TWOu = 2.0*URO
  FOUru = 4.0*URO
END SUBROUTINE MACON
