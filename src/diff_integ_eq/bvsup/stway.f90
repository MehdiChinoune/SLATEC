!** STWAY
SUBROUTINE STWAY(U,V,Yhp,Inout,Stowa)
  !>
  !  Subsidiary to BVSUP
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (STWAY-S, DSTWAY-D)
  !***
  ! **Author:**  Watts, H. A., (SNLA)
  !***
  ! **Description:**
  !
  !  This subroutine stores (recalls) integration data in the event
  !  that a restart is needed (the homogeneous solution vectors become
  !  too dependent to continue)
  !
  !***
  ! **See also:**  BVSUP
  !***
  ! **Routines called:**  STOR1
  !***
  ! COMMON BLOCKS    ML15TO, ML18JR, ML8SZ

  !* REVISION HISTORY  (YYMMDD)
  !   750601  DATE WRITTEN
  !   890921  Realigned order of variables in certain COMMON blocks.
  !           (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910722  Updated AUTHOR section.  (ALS)
  USE ML, ONLY: ncomp_com, nfc_com, x_com, xop_com, info_com, istkop_com, kop_com, &
    ndisk_com, ntape_com, neqivp_com
  INTEGER :: Inout
  REAL(SP) :: Stowa(:), U(:), V(:), Yhp(:)
  INTEGER :: j, k, ko, ks, ksj
  !
  !* FIRST EXECUTABLE STATEMENT  STWAY
  IF ( Inout==1 ) THEN
    !
    !     RECALL FROM STOWA ARRAY AND ISTKOP
    !
    ks = nfc_com*ncomp_com
    CALL STOR1(Yhp,Stowa,Yhp(ks+1:),Stowa(ks+1:),1,0,0)
    ks = ks + ncomp_com
    IF ( neqivp_com/=0 ) THEN
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
    IF ( ndisk_com==0.OR.ko==0 ) RETURN
    DO k = 1, ko
      BACKSPACE ntape_com
    END DO
  ELSE
    !
    !     SAVE IN STOWA ARRAY AND ISTKOP
    !
    ks = nfc_com*ncomp_com
    CALL STOR1(Stowa,U,Stowa(ks+1:),V,1,0,0)
    ks = ks + ncomp_com
    IF ( neqivp_com/=0 ) THEN
      DO j = 1, neqivp_com
        ksj = ks + j
        Stowa(ksj) = Yhp(ksj)
      END DO
    END IF
    ks = ks + neqivp_com
    Stowa(ks+1) = x_com
    istkop_com = kop_com
    IF ( xop_com==x_com ) istkop_com = kop_com + 1
    RETURN
  END IF
END SUBROUTINE STWAY
