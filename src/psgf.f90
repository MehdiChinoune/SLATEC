!*==PSGF.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK PSGF
      FUNCTION PSGF(X,Iz,C,A,Bh)
      IMPLICIT NONE
!*--PSGF5
!*** Start of declarations inserted by SPAG
      REAL A , Bh , C , dd , fsg , hsg , PSGF , X
      INTEGER Iz , j
!*** End of declarations inserted by SPAG
!***BEGIN PROLOGUE  PSGF
!***SUBSIDIARY
!***PURPOSE  Subsidiary to BLKTRI
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (PSGF-S)
!***AUTHOR  (UNKNOWN)
!***SEE ALSO  BLKTRI
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  PSGF
      DIMENSION A(*) , C(*) , Bh(*)
!***FIRST EXECUTABLE STATEMENT  PSGF
      fsg = 1.
      hsg = 1.
      DO j = 1 , Iz
        dd = 1./(X-Bh(j))
        fsg = fsg*A(j)*dd
        hsg = hsg*C(j)*dd
      ENDDO
      IF ( MOD(Iz,2)/=0 ) THEN
        PSGF = 1. + fsg + hsg
        GOTO 99999
      ENDIF
      PSGF = 1. - fsg - hsg
      RETURN
99999 END FUNCTION PSGF
