!** MC20AS
SUBROUTINE MC20AS(Nc,Maxa,A,Inum,Jptr,Jnum,Jdisp)
  IMPLICIT NONE
  !>
  !***
  !  Subsidiary to SPLP
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (MC20AS-S, MC20AD-D)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !     THIS SUBPROGRAM IS A SLIGHT MODIFICATION OF A SUBPROGRAM
  !     FROM THE C. 1979 AERE HARWELL LIBRARY.  THE NAME OF THE
  !     CORRESPONDING HARWELL CODE CAN BE OBTAINED BY DELETING
  !     THE FINAL LETTER =S= IN THE NAMES USED HERE.
  !     REVISED SEP. 13, 1979.
  !
  !     ROYALTIES HAVE BEEN PAID TO AERE-UK FOR USE OF THEIR CODES
  !     IN THE PACKAGE GIVEN HERE.  ANY PRIMARY USAGE OF THE HARWELL
  !     SUBROUTINES REQUIRES A ROYALTY AGREEMENT AND PAYMENT BETWEEN
  !     THE USER AND AERE-UK.  ANY USAGE OF THE SANDIA WRITTEN CODES
  !     SPLP( ) (WHICH USES THE HARWELL SUBROUTINES) IS PERMITTED.
  !
  !***
  ! **See also:**  SPLP
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   811215  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  
  REAL ace, acep
  INTEGER i, ice, icep, j, ja, jb, jce, jcep, Jdisp, Jptr, k, &
    kr, loc, Maxa, Nc, null
  INTEGER Inum(*), Jnum(*)
  REAL A(*)
  DIMENSION Jptr(Nc)
  !* FIRST EXECUTABLE STATEMENT  MC20AS
  null = -Jdisp
  !**      CLEAR JPTR
  DO j = 1, Nc
    Jptr(j) = 0
  ENDDO
  !**      COUNT THE NUMBER OF ELEMENTS IN EACH COLUMN.
  DO k = 1, Maxa
    j = Jnum(k) + Jdisp
    Jptr(j) = Jptr(j) + 1
  ENDDO
  !**      SET THE JPTR ARRAY
  k = 1
  DO j = 1, Nc
    kr = k + Jptr(j)
    Jptr(j) = k
    k = kr
  ENDDO
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
      Jnum(i) = null
      !        CHAIN FROM CURRENT ENTRY TO STORE ITEMS.
      DO j = 1, Maxa
        !        CURRENT ENTRY NOT IN CORRECT POSITION.  DETERMINE CORRECT
        !        POSITION TO STORE ENTRY.
        loc = Jptr(jce)
        Jptr(jce) = Jptr(jce) + 1
        !        SAVE CONTENTS OF THAT LOCATION.
        acep = A(loc)
        icep = Inum(loc)
        jcep = Jnum(loc)
        !        STORE CURRENT ENTRY.
        A(loc) = ace
        Inum(loc) = ice
        Jnum(loc) = null
        !        CHECK IF NEXT CURRENT ENTRY NEEDS TO BE PROCESSED.
        IF ( jcep==null ) EXIT
        !        IT DOES.  COPY INTO CURRENT ENTRY.
        ace = acep
        ice = icep
        jce = jcep + Jdisp
      ENDDO
    ENDIF
    !
  ENDDO
  !
  !**      RESET JPTR VECTOR.
  ja = 1
  DO j = 1, Nc
    jb = Jptr(j)
    Jptr(j) = ja
    ja = jb
  ENDDO
END SUBROUTINE MC20AS
