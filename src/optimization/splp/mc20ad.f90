!** MC20AD
SUBROUTINE MC20AD(Nc,Maxa,A,Inum,Jptr,Jnum,Jdisp)
  !>
  !  Subsidiary to DSPLP
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      DOUBLE PRECISION (MC20AS-S, MC20AD-D)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !     THIS SUBPROGRAM IS A SLIGHT MODIFICATION OF A SUBPROGRAM
  !     FROM THE C. 1979 AERE HARWELL LIBRARY.  THE NAME OF THE
  !     CORRESPONDING HARWELL CODE CAN BE OBTAINED BY DELETING
  !     THE FINAL LETTER =D= IN THE NAMES USED HERE.
  !     REVISED SEP. 13, 1979.
  !
  !     ROYALTIES HAVE BEEN PAID TO AERE-UK FOR USE OF THEIR CODES
  !     IN THE PACKAGE GIVEN HERE.  ANY PRIMARY USAGE OF THE HARWELL
  !     SUBROUTINES REQUIRES A ROYALTY AGREEMENT AND PAYMENT BETWEEN
  !     THE USER AND AERE-UK.  ANY USAGE OF THE SANDIA WRITTEN CODES
  !     DSPLP( ) (WHICH USES THE HARWELL SUBROUTINES) IS PERMITTED.
  !
  !***
  ! **See also:**  DSPLP
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   811215  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)

  INTEGER :: Nc, Jdisp, Jptr(Nc), Maxa
  INTEGER :: Inum(Maxa), Jnum(Maxa)
  REAL(8) :: A(Maxa)
  INTEGER :: i, ice, icep, j, ja, jb, jce, jcep, k, kr, locc, nul
  REAL(8) :: ace, acep
  !* FIRST EXECUTABLE STATEMENT  MC20AD
  nul = -Jdisp
  !**      CLEAR JPTR
  DO j = 1, Nc
    Jptr(j) = 0
  END DO
  !**      COUNT THE NUMBER OF ELEMENTS IN EACH COLUMN.
  DO k = 1, Maxa
    j = Jnum(k) + Jdisp
    Jptr(j) = Jptr(j) + 1
  END DO
  !**      SET THE JPTR ARRAY
  k = 1
  DO j = 1, Nc
    kr = k + Jptr(j)
    Jptr(j) = k
    k = kr
  END DO
  !
  !**      REORDER THE ELEMENTS INTO COLUMN ORDER.  THE ALGORITHM IS AN
  !        IN-PLACE SORT AND IS OF ORDER MAXA.
  DO i = 1, Maxa
    !        ESTABLISH THE CURRENT ENTRY.
    jce = Jnum(i) + Jdisp
    IF ( jce/=0 ) THEN
      ace = A(i)
      ice = Inum(i)
      !        CLEAR THE LOCATION VACATED.
      Jnum(i) = nul
      !        CHAIN FROM CURRENT ENTRY TO STORE ITEMS.
      DO j = 1, Maxa
        !        CURRENT ENTRY NOT IN CORRECT POSITION.  DETERMINE CORRECT
        !        POSITION TO STORE ENTRY.
        locc = Jptr(jce)
        Jptr(jce) = Jptr(jce) + 1
        !        SAVE CONTENTS OF THAT LOCATION.
        acep = A(locc)
        icep = Inum(locc)
        jcep = Jnum(locc)
        !        STORE CURRENT ENTRY.
        A(locc) = ace
        Inum(locc) = ice
        Jnum(locc) = nul
        !        CHECK IF NEXT CURRENT ENTRY NEEDS TO BE PROCESSED.
        IF ( jcep==nul ) EXIT
        !        IT DOES.  COPY INTO CURRENT ENTRY.
        ace = acep
        ice = icep
        jce = jcep + Jdisp
      END DO
    END IF
    !
  END DO
  !
  !**      RESET JPTR VECTOR.
  ja = 1
  DO j = 1, Nc
    jb = Jptr(j)
    Jptr(j) = ja
    ja = jb
  END DO
END SUBROUTINE MC20AD
