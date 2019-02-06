!*==INXCA.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK INXCA
      SUBROUTINE INXCA(I,Ir,Idxa,Na)
      IMPLICIT NONE
!*--INXCA5
!*** Start of declarations inserted by SPAG
      REAL CNV , EPS
      INTEGER I , Idxa , IK , Ir , K , Na , NCMplx , NM , NPP
!*** End of declarations inserted by SPAG
!***BEGIN PROLOGUE  INXCA
!***SUBSIDIARY
!***PURPOSE  Subsidiary to CBLKTR
!***LIBRARY   SLATEC
!***TYPE      INTEGER (INXCA-I)
!***AUTHOR  (UNKNOWN)
!***SEE ALSO  CBLKTR
!***ROUTINES CALLED  (NONE)
!***COMMON BLOCKS    CCBLK
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  INXCA
      COMMON /CCBLK / NPP , K , EPS , CNV , NM , NCMplx , IK
!***FIRST EXECUTABLE STATEMENT  INXCA
      Na = 2**Ir
      Idxa = I - Na + 1
      IF ( I>NM ) Na = 0
      END SUBROUTINE INXCA
