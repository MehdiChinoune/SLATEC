!*==CHKPR4.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK CHKPR4
SUBROUTINE CHKPR4(Iorder,A,B,M,Mbdcnd,C,D,N,Nbdcnd,COFX,Idmn,Ierror)
  IMPLICIT NONE
  !*--CHKPR45
  !*** Start of declarations inserted by SPAG
  REAL A , ai , B , bi , C , ci , D , dlx , xi
  INTEGER i , Idmn , Ierror , Iorder , M , Mbdcnd , N , Nbdcnd
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  CHKPR4
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to SEPX4
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (CHKPR4-S)
  !***AUTHOR  (UNKNOWN)
  !***DESCRIPTION
  !
  !     This program checks the input parameters for errors.
  !
  !***SEE ALSO  SEPX4
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !***END PROLOGUE  CHKPR4
  EXTERNAL COFX
  !***FIRST EXECUTABLE STATEMENT  CHKPR4
  Ierror = 1
  IF ( A>=B.OR.C>=D ) RETURN
  !
  !     CHECK BOUNDARY SWITCHES
  !
  Ierror = 2
  IF ( Mbdcnd<0.OR.Mbdcnd>4 ) RETURN
  Ierror = 3
  IF ( Nbdcnd<0.OR.Nbdcnd>4 ) RETURN
  !
  !     CHECK FIRST DIMENSION IN CALLING ROUTINE
  !
  Ierror = 5
  IF ( Idmn<7 ) RETURN
  !
  !     CHECK M
  !
  Ierror = 6
  IF ( M>(Idmn-1).OR.M<6 ) RETURN
  !
  !     CHECK N
  !
  Ierror = 7
  IF ( N<5 ) RETURN
  !
  !     CHECK IORDER
  !
  Ierror = 8
  IF ( Iorder/=2.AND.Iorder/=4 ) RETURN
  !
  !     CHECK THAT EQUATION IS ELLIPTIC
  !
  dlx = (B-A)/M
  DO i = 2 , M
    xi = A + (i-1)*dlx
    CALL COFX(xi,ai,bi,ci)
    IF ( ai<=0.0 ) THEN
      Ierror = 10
      RETURN
    ENDIF
  ENDDO
  !
  !     NO ERROR FOUND
  !
  Ierror = 0
END SUBROUTINE CHKPR4
