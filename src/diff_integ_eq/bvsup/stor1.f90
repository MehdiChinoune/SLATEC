!** STOR1
SUBROUTINE STOR1(U,Yh,V,Yp,Ntemp,Ndisk,Ntape)
  !> Subsidiary to BVSUP
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (STOR1-S, DSTOR1-D)
  !***
  ! **Author:**  Watts, H. A., (SNLA)
  !***
  ! **Description:**
  !
  !- *********************************************************************
  !             0 -- Storage at output points.
  !     NTEMP =
  !             1 -- Temporary storage
  !- *********************************************************************
  !
  !***
  ! **See also:**  BVSUP
  !***
  ! **Routines called:**  (NONE)
  !***
  ! COMMON BLOCKS    ML8SZ

  !* REVISION HISTORY  (YYMMDD)
  !   750601  DATE WRITTEN
  !   890921  Realigned order of variables in certain COMMON blocks.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910722  Updated AUTHOR section.  (ALS)
  USE ML, ONLY : c_com, inhomo_com, ncomp_com, nfc_com
  !
  INTEGER, INTENT(IN) :: Ndisk, Ntape, Ntemp
  REAL(SP), INTENT(IN) :: Yh(:), Yp(:)
  REAL(SP), INTENT(OUT) :: U(:), V(:)
  !
  INTEGER :: j, nctnf
  !
  !- *********************************************************************
  !
  !* FIRST EXECUTABLE STATEMENT  STOR1
  nctnf = ncomp_com*nfc_com
  DO j = 1, nctnf
    U(j) = Yh(j)
  END DO
  IF( inhomo_com/=1 ) THEN
    !
    !   ZERO PARTICULAR SOLUTION
    !
    IF( Ntemp==1 ) RETURN
    DO j = 1, ncomp_com
      V(j) = 0._SP
    END DO
    IF( Ndisk==1 ) WRITE (Ntape) (V(j),j=1,ncomp_com), (U(j),j=1,nctnf)
    !
    !   NONZERO PARTICULAR SOLUTION
    !
  ELSEIF( Ntemp==0 ) THEN
    !
    DO j = 1, ncomp_com
      V(j) = c_com*Yp(j)
    END DO
    !
    !  IS OUTPUT INFORMATION TO BE WRITTEN TO DISK
    !
    IF( Ndisk==1 ) WRITE (Ntape) (V(j),j=1,ncomp_com), (U(j),j=1,nctnf)
  ELSE
    !
    DO j = 1, ncomp_com
      V(j) = Yp(j)
    END DO
    RETURN
  END IF
  !
END SUBROUTINE STOR1