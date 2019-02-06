!*==PPGSF.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK PPGSF
      FUNCTION PPGSF(X,Iz,C,A,Bh)
      IMPLICIT NONE
!*--PPGSF5
!*** Start of declarations inserted by SPAG
      REAL A , Bh , C , PPGSF , sum , X
      INTEGER Iz , j
!*** End of declarations inserted by SPAG
!***BEGIN PROLOGUE  PPGSF
!***SUBSIDIARY
!***PURPOSE  Subsidiary to CBLKTR
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (PPGSF-S)
!***AUTHOR  (UNKNOWN)
!***SEE ALSO  CBLKTR
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  PPGSF
      DIMENSION A(*) , C(*) , Bh(*)
!***FIRST EXECUTABLE STATEMENT  PPGSF
      sum = 0.
      DO j = 1 , Iz
        sum = sum - 1./(X-Bh(j))**2
      ENDDO
      PPGSF = sum
      END FUNCTION PPGSF
