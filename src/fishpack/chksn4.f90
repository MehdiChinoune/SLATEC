!** CHKSN4
SUBROUTINE CHKSN4(Mbdcnd,Nbdcnd,Alpha,Beta,COFX,Singlr)
  IMPLICIT NONE
  !>
  !***
  !  Subsidiary to SEPX4
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (CHKSN4-S)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !     This subroutine checks if the PDE SEPX4
  !     must solve is a singular operator.
  !
  !***
  ! **See also:**  SEPX4
  !***
  ! **Routines called:**  (NONE)
  !***
  ! COMMON BLOCKS    SPL4

  !* REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  
  REAL ai, AIT, Alpha, Beta, bi, BIT, ci, CIT, DIT, DLX, DLX4, &
    DLY, DLY4, TDLx3, TDLy3, xi
  INTEGER i, IS, JS, K, KSWx, KSWy, L, Mbdcnd, MIT, MS, Nbdcnd, NIT, NS
  COMMON /SPL4  / KSWx, KSWy, K, L, AIT, BIT, CIT, DIT, MIT, NIT, &
    IS, MS, JS, NS, DLX, DLY, TDLx3, TDLy3, DLX4, DLY4
  LOGICAL Singlr
  EXTERNAL :: COFX
  !* FIRST EXECUTABLE STATEMENT  CHKSN4
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
  !
  !     CHECK THAT NON-DERIVATIVE COEFFICIENT FUNCTIONS
  !     ARE ZERO
  !
  DO i = IS, MS
    xi = AIT + (i-1)*DLX
    CALL COFX(xi,ai,bi,ci)
    IF ( ci/=0.0 ) RETURN
  ENDDO
  !
  !     THE OPERATOR MUST BE SINGULAR IF THIS POINT IS REACHED
  !
  Singlr = .TRUE.
END SUBROUTINE CHKSN4
