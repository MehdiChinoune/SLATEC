!*==DMACON.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK DMACON
      SUBROUTINE DMACON
      IMPLICIT NONE
!*--DMACON5
!***BEGIN PROLOGUE  DMACON
!***SUBSIDIARY
!***PURPOSE  Subsidiary to DBVSUP
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (MACON-S, DMACON-D)
!***AUTHOR  (UNKNOWN)
!***SEE ALSO  DBVSUP
!***ROUTINES CALLED  D1MACH
!***COMMON BLOCKS    DML5MC
!***REVISION HISTORY  (YYMMDD)
!   750601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890921  Realigned order of variables in certain COMMON blocks.
!           (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!***END PROLOGUE  DMACON
      DOUBLE PRECISION D1MACH
      INTEGER ke , LPAr
      DOUBLE PRECISION dd , EPS , FOUru , SQOvfl , SRU , TWOu , URO
      COMMON /DML5MC/ URO , SRU , EPS , SQOvfl , TWOu , FOUru , LPAr
!***FIRST EXECUTABLE STATEMENT  DMACON
      URO = D1MACH(4)
      SRU = SQRT(URO)
      dd = -LOG10(URO)
      LPAr = 0.5D0*dd
      ke = 0.5D0 + 0.75D0*dd
      EPS = 10.0D0**(-2*ke)
      SQOvfl = SQRT(D1MACH(2))
      TWOu = 2.0D0*URO
      FOUru = 4.0D0*URO
      END SUBROUTINE DMACON
