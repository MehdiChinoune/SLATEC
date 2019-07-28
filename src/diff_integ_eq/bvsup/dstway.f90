!** DSTWAY
SUBROUTINE DSTWAY(U,V,Yhp,Inout,Stowa)
  !> Subsidiary to DBVSUP
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      DOUBLE PRECISION (STWAY-S, DSTWAY-D)
  !***
  ! **Author:**  Watts, H. A., (SNLA)
  !***
  ! **Description:**
  !
  !  This subroutine stores (recalls) integration data in the event
  !  that a restart is needed (the homogeneous solution vectors become
  !  too dependent to continue).
  !
  !***
  ! **See also:**  DBVSUP
  !***
  ! **Routines called:**  DSTOR1
  !***
  ! COMMON BLOCKS    DML15T, DML18J, DML8SZ

  !* REVISION HISTORY  (YYMMDD)
  !   750601  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890921  Realigned order of variables in certain COMMON blocks.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910722  Updated AUTHOR section.  (ALS)
  USE DML, ONLY : ncomp_com, nfc_com, x_com, xop_com, info_com, istkop_com, &
    kop_com, ndisk_com, ntape_com, neqivp_com
  !
  INTEGER, INTENT(IN) :: Inout
  REAL(DP), INTENT(IN) :: U(:), V(:)
  REAL(DP), INTENT(INOUT) :: Stowa(:), Yhp(:)
  !
  INTEGER :: j, k, ko, ks, ksj
  !* FIRST EXECUTABLE STATEMENT  DSTWAY
  IF( Inout==1 ) THEN
    !
    !        RECALL FROM STOWA ARRAY AND ISTKOP
    !
    ks = nfc_com*ncomp_com
    CALL DSTOR1(Yhp,Stowa,Yhp(ks+1:),Stowa(ks+1:),1,0,0)
    ks = ks + ncomp_com
    IF( neqivp_com>=1 ) THEN
      DO j = 1, neqivp_com
        ksj = ks + j
        Yhp(ksj) = Stowa(ksj)
      END DO
    END IF
    ks = ks + neqivp_com
    x_com = Stowa(ks+1)
    info_com(1) = 0
    ko = kop_com - istkop_com
    kop_com = istkop_com
    IF( ndisk_com/=0 .AND. ko/=0 ) THEN
      DO k = 1, ko
        BACKSPACE ntape_com
      END DO
    END IF
  ELSE
    !
    !        SAVE IN STOWA ARRAY AND ISTKOP
    !
    ks = nfc_com*ncomp_com
    CALL DSTOR1(Stowa,U,Stowa(ks+1:),V,1,0,0)
    ks = ks + ncomp_com
    IF( neqivp_com>=1 ) THEN
      DO j = 1, neqivp_com
        ksj = ks + j
        Stowa(ksj) = Yhp(ksj)
      END DO
    END IF
    ks = ks + neqivp_com
    Stowa(ks+1) = x_com
    istkop_com = kop_com
    IF( xop_com==x_com ) istkop_com = kop_com + 1
  END IF
  !
END SUBROUTINE DSTWAY