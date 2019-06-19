!** CHKSNG
SUBROUTINE CHKSNG(Mbdcnd,Nbdcnd,Alpha,Beta,Gama,Xnu,COFX,COFY,Singlr)
  !> Subsidiary to SEPELI
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
  USE SPLPCM, ONLY : ait_com, cit_com, dlx_com, dly_com, is_com, js_com, ms_com, ns_com
  INTERFACE
    SUBROUTINE COFX(X,A,B,C)
      IMPORT SP
      REAL(SP) :: X, A, B, C
    END SUBROUTINE COFX
    SUBROUTINE COFY(Y,D,E,F)
      IMPORT SP
      REAL(SP) :: Y, D, E, F
    END SUBROUTINE COFY
  END INTERFACE
  INTEGER :: Mbdcnd, Nbdcnd
  REAL(SP) :: Alpha, Beta, Gama, Xnu
  LOGICAL :: Singlr
  INTEGER :: i, j
  REAL(SP) :: ai, bi, ci, dj, ej, fj, xi, yj
  !* FIRST EXECUTABLE STATEMENT  CHKSNG
  Singlr = .FALSE.
  !
  !     CHECK IF THE BOUNDARY CONDITIONS ARE
  !     ENTIRELY PERIODIC AND/OR MIXED
  !
  IF( (Mbdcnd/=0 .AND. Mbdcnd/=3) .OR. (Nbdcnd/=0 .AND. Nbdcnd/=3) ) RETURN
  !
  !     CHECK THAT MIXED CONDITIONS ARE PURE NEUMAN
  !
  IF( Mbdcnd==3 ) THEN
    IF( Alpha/=0.0 .OR. Beta/=0.0 ) RETURN
  END IF
  IF( Nbdcnd==3 ) THEN
    IF( Gama/=0.0 .OR. Xnu/=0.0 ) RETURN
  END IF
  !
  !     CHECK THAT NON-DERIVATIVE COEFFICIENT FUNCTIONS
  !     ARE ZERO
  !
  DO i = is_com, ms_com
    xi = ait_com + (i-1)*dlx_com
    CALL COFX(xi,ai,bi,ci)
    IF( ci/=0.0 ) RETURN
  END DO
  DO j = js_com, ns_com
    yj = cit_com + (j-1)*dly_com
    CALL COFY(yj,dj,ej,fj)
    IF( fj/=0.0 ) RETURN
  END DO
  !
  !     THE OPERATOR MUST BE SINGULAR IF THIS POINT IS REACHED
  !
  Singlr = .TRUE.
END SUBROUTINE CHKSNG
