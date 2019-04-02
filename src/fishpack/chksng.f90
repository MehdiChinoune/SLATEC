!** CHKSNG
SUBROUTINE CHKSNG(Mbdcnd,Nbdcnd,Alpha,Beta,Gama,Xnu,COFX,COFY,Singlr)
  USE SPLPCM
  IMPLICIT NONE
  !>
  !***
  !  Subsidiary to SEPELI
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (CHKSNG-S)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !     This subroutine checks if the PDE SEPELI
  !     must solve is a singular operator.
  !
  !***
  ! **See also:**  SEPELI
  !***
  ! **Routines called:**  (NONE)
  !***
  ! COMMON BLOCKS    SPLPCM

  !* REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)

  REAL ai, Alpha, Beta, bi, ci, dj, ej, fj, Gama, xi, Xnu, yj
  INTEGER i, j, Mbdcnd, Nbdcnd
  LOGICAL Singlr
  !* FIRST EXECUTABLE STATEMENT  CHKSNG
  Singlr = .FALSE.
  !
  !     CHECK IF THE BOUNDARY CONDITIONS ARE
  !     ENTIRELY PERIODIC AND/OR MIXED
  !
  IF ( (Mbdcnd/=0.AND.Mbdcnd/=3).OR.(Nbdcnd/=0.AND.Nbdcnd/=3) ) RETURN
  !
  !     CHECK THAT MIXED CONDITIONS ARE PURE NEUMAN
  !
  IF ( Mbdcnd==3 ) THEN
    IF ( Alpha/=0.0.OR.Beta/=0.0 ) RETURN
  ENDIF
  IF ( Nbdcnd==3 ) THEN
    IF ( Gama/=0.0.OR.Xnu/=0.0 ) RETURN
  ENDIF
  !
  !     CHECK THAT NON-DERIVATIVE COEFFICIENT FUNCTIONS
  !     ARE ZERO
  !
  DO i = IS, MS
    xi = AIT + (i-1)*DLX
    CALL COFX(xi,ai,bi,ci)
    IF ( ci/=0.0 ) RETURN
  ENDDO
  DO j = JS, NS
    yj = CIT + (j-1)*DLY
    CALL COFY(yj,dj,ej,fj)
    IF ( fj/=0.0 ) RETURN
  ENDDO
  !
  !     THE OPERATOR MUST BE SINGULAR IF THIS POINT IS REACHED
  !
  Singlr = .TRUE.
END SUBROUTINE CHKSNG
