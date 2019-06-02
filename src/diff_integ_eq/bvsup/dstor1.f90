!** DSTOR1
SUBROUTINE DSTOR1(U,Yh,V,Yp,Ntemp,Ndisk,Ntape)
  !>
  !  Subsidiary to DBVSUP
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      DOUBLE PRECISION (STOR1-S, DSTOR1-D)
  !***
  ! **Author:**  Watts, H. A., (SNLA)
  !***
  ! **Description:**
  !
  !- *********************************************************************
  !             0 -- storage at output points.
  !     NTEMP =
  !             1 -- temporary storage
  !- *********************************************************************
  !
  !***
  ! **See also:**  DBVSUP
  !***
  ! **Routines called:**  (NONE)
  !***
  ! COMMON BLOCKS    DML8SZ

  !* REVISION HISTORY  (YYMMDD)
  !   750601  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890921  Realigned order of variables in certain COMMON blocks.
  !           (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910722  Updated AUTHOR section.  (ALS)
  USE DML, ONLY : c_com, inhomo_com, ncomp_com, nfc_com
  INTEGER :: Ndisk, Ntape, Ntemp
  REAL(8) :: U(:), V(:), Yh(:), Yp(:)
  INTEGER :: j, nctnf
  !     BEGIN BLOCK PERMITTING ...EXITS TO 80
  !* FIRST EXECUTABLE STATEMENT  DSTOR1
  nctnf = ncomp_com*nfc_com
  DO j = 1, nctnf
    U(j) = Yh(j)
  END DO
  IF ( inhomo_com==1 ) THEN
    !
    !           NONZERO PARTICULAR SOLUTION
    !
    IF ( Ntemp==0 ) THEN
      !
      DO j = 1, ncomp_com
        V(j) = c_com*Yp(j)
      END DO
      !
      !        IS OUTPUT INFORMATION TO BE WRITTEN TO DISK
      !
      IF ( Ndisk==1 ) WRITE (Ntape) (V(j),j=1,ncomp_com), (U(j),j=1,nctnf)
    ELSE
      !
      DO j = 1, ncomp_com
        V(j) = Yp(j)
        !     .........EXIT
      END DO
    END IF
    !
    !           ZERO PARTICULAR SOLUTION
    !
    !     ......EXIT
  ELSEIF ( Ntemp/=1 ) THEN
    DO j = 1, ncomp_com
      V(j) = 0.0D0
    END DO
    IF ( Ndisk==1 ) WRITE (Ntape) (V(j),j=1,ncomp_com), (U(j),j=1,nctnf)
  END IF
  !
END SUBROUTINE DSTOR1
