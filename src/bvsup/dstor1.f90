!** DSTOR1
SUBROUTINE DSTOR1(U,Yh,V,Yp,Ntemp,Ndisk,Ntape)
  USE DML, ONLY : C, INHomo, NCOmp, NFC
  IMPLICIT NONE
  !>
  !***
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

  INTEGER j, nctnf, Ndisk, Ntape, Ntemp
  REAL(8) :: U(*), V(*), Yh(*), Yp(*)
  !     BEGIN BLOCK PERMITTING ...EXITS TO 80
  !* FIRST EXECUTABLE STATEMENT  DSTOR1
  nctnf = NCOmp*NFC
  DO j = 1, nctnf
    U(j) = Yh(j)
  ENDDO
  IF ( INHomo==1 ) THEN
    !
    !           NONZERO PARTICULAR SOLUTION
    !
    IF ( Ntemp==0 ) THEN
      !
      DO j = 1, NCOmp
        V(j) = C*Yp(j)
      ENDDO
      !
      !        IS OUTPUT INFORMATION TO BE WRITTEN TO DISK
      !
      IF ( Ndisk==1 ) WRITE (Ntape) (V(j),j=1,NCOmp), (U(j),j=1,nctnf)
    ELSE
      !
      DO j = 1, NCOmp
        V(j) = Yp(j)
        !     .........EXIT
      ENDDO
    ENDIF
    !
    !           ZERO PARTICULAR SOLUTION
    !
    !     ......EXIT
  ELSEIF ( Ntemp/=1 ) THEN
    DO j = 1, NCOmp
      V(j) = 0.0D0
    ENDDO
    IF ( Ndisk==1 ) WRITE (Ntape) (V(j),j=1,NCOmp), (U(j),j=1,nctnf)
  ENDIF
  !
END SUBROUTINE DSTOR1
