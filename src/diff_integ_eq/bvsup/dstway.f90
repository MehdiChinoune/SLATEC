!** DSTWAY
SUBROUTINE DSTWAY(U,V,Yhp,Inout,Stowa)
  !>
  !  Subsidiary to DBVSUP
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
  !   890921  Realigned order of variables in certain COMMON blocks.
  !           (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910722  Updated AUTHOR section.  (ALS)
  USE DML, ONLY : NCOmp, NFC, X, XOP, INFo, ISTkop, KOP, NDIsk, NTApe, NEQivp
  !
  INTEGER :: Inout
  REAL(8) :: Stowa(:), U(:), V(:), Yhp(:)
  INTEGER :: j, k, ko, ks, ksj
  !* FIRST EXECUTABLE STATEMENT  DSTWAY
  IF ( Inout==1 ) THEN
    !
    !        RECALL FROM STOWA ARRAY AND ISTKOP
    !
    ks = NFC*NCOmp
    CALL DSTOR1(Yhp,Stowa,Yhp(ks+1:),Stowa(ks+1:),1,0,0)
    ks = ks + NCOmp
    IF ( NEQivp>=1 ) THEN
      DO j = 1, NEQivp
        ksj = ks + j
        Yhp(ksj) = Stowa(ksj)
      END DO
    END IF
    ks = ks + NEQivp
    X = Stowa(ks+1)
    INFo(1) = 0
    ko = KOP - ISTkop
    KOP = ISTkop
    IF ( NDIsk/=0.AND.ko/=0 ) THEN
      DO k = 1, ko
        BACKSPACE NTApe
      END DO
    END IF
  ELSE
    !
    !        SAVE IN STOWA ARRAY AND ISTKOP
    !
    ks = NFC*NCOmp
    CALL DSTOR1(Stowa,U,Stowa(ks+1:),V,1,0,0)
    ks = ks + NCOmp
    IF ( NEQivp>=1 ) THEN
      DO j = 1, NEQivp
        ksj = ks + j
        Stowa(ksj) = Yhp(ksj)
      END DO
    END IF
    ks = ks + NEQivp
    Stowa(ks+1) = X
    ISTkop = KOP
    IF ( XOP==X ) ISTkop = KOP + 1
  END IF
END SUBROUTINE DSTWAY
