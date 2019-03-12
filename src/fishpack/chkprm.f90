!DECK CHKPRM
SUBROUTINE CHKPRM(Intl,Iorder,A,B,M,Mbdcnd,C,D,N,Nbdcnd,COFX,COFY,Idmn,&
    Ierror)
  IMPLICIT NONE
  REAL A, ai, B, bi, C, ci, D, dj, dlx, dly, ej, fj, xi, yj
  INTEGER i, Idmn, Ierror, Intl, Iorder, j, M, Mbdcnd, N, Nbdcnd
  !***BEGIN PROLOGUE  CHKPRM
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to SEPELI
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (CHKPRM-S)
  !***AUTHOR  (UNKNOWN)
  !***DESCRIPTION
  !
  !     This program checks the input parameters for errors.
  !
  !***SEE ALSO  SEPELI
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !***END PROLOGUE  CHKPRM
  !
  EXTERNAL COFX, COFY
  !***FIRST EXECUTABLE STATEMENT  CHKPRM
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
  !     CHECK INTL
  !
  Ierror = 9
  IF ( Intl/=0.AND.Intl/=1 ) RETURN
  !
  !     CHECK THAT EQUATION IS ELLIPTIC
  !
  dlx = (B-A)/M
  dly = (D-C)/N
  DO i = 2, M
    xi = A + (i-1)*dlx
    CALL COFX(xi,ai,bi,ci)
    DO j = 2, N
      yj = C + (j-1)*dly
      CALL COFY(yj,dj,ej,fj)
      IF ( ai*dj<=0.0 ) THEN
        Ierror = 10
        RETURN
      ENDIF
    ENDDO
  ENDDO
  !
  !     NO ERROR FOUND
  !
  Ierror = 0
END SUBROUTINE CHKPRM
