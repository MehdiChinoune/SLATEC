!** PRWPGE
SUBROUTINE PRWPGE(Key,Ipage,Lpg,Sx,Ix)
  !> Subsidiary to SPLP
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (PRWPGE-S, DPRWPG-D)
  !***
  ! **Author:**  Hanson, R. J., (SNLA)
  !           Wisniewski, J. A., (SNLA)
  !***
  ! **Description:**
  !
  !     PRWPGE LIMITS THE TYPE OF STORAGE TO A SEQUENTIAL SCHEME.
  !     VIRTUAL MEMORY PAGE READ/WRITE SUBROUTINE.
  !
  !     DEPENDING ON THE VALUE OF KEY, SUBROUTINE PRWPGE() PERFORMS A PAGE
  !     READ OR WRITE OF PAGE IPAGE. THE PAGE HAS LENGTH LPG.
  !
  !     KEY       IS A FLAG INDICATING WHETHER A PAGE READ OR WRITE IS
  !               TO BE PERFORMED.
  !               IF KEY = 1 DATA IS READ.
  !               IF KEY = 2 DATA IS WRITTEN.
  !     IPAGE     IS THE PAGE NUMBER OF THE MATRIX TO BE ACCESSED.
  !     LPG       IS THE LENGTH OF THE PAGE OF THE MATRIX TO BE ACCESSED.
  !   SX(*),IX(*) IS THE MATRIX TO BE ACCESSED.
  !
  !     THIS SUBROUTINE IS A MODIFICATION OF THE SUBROUTINE LRWPGE,
  !     SANDIA LABS. REPT. SAND78-0785.
  !     MODIFICATIONS BY K.L. HIEBERT AND R.J. HANSON
  !     REVISED 811130-1000
  !     REVISED YYMMDD-HHMM
  !
  !***
  ! **See also:**  SPLP
  !***
  ! **Routines called:**  PRWVIR, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   811215  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900328  Added TYPE section.  (WRB)
  !   900510  Fixed error messages and replaced GOTOs with
  !           IF-THEN-ELSE.  (RWC)
  !   910403  Updated AUTHOR and DESCRIPTION sections.  (WRB)
  USE service, ONLY : XERMSG
  INTEGER :: Ipage, Lpg, Ix(Lpg), Key
  REAL(SP) :: Sx(Lpg)
  !* FIRST EXECUTABLE STATEMENT  PRWPGE
  !
  !     CHECK IF IPAGE IS IN RANGE.
  !
  IF( Ipage<1 ) CALL XERMSG('PRWPGE',&
    'THE VALUE OF IPAGE (PAGE NUMBER) WAS NOT IN THE RANGE1<=IPAGE<=MAXPGE.',55,1)
  !
  !     CHECK IF LPG IS POSITIVE.
  !
  IF( Lpg<=0 ) CALL XERMSG('PRWPGE',&
    'THE VALUE OF LPG (PAGE LENGTH) WAS NONPOSITIVE.',55,1)
  !
  !     DECIDE IF WE ARE READING OR WRITING.
  !
  IF( Key==1 ) THEN
    !
    !        CODE TO DO A PAGE READ.
    !
    CALL PRWVIR(Key,Ipage,Lpg,Sx,Ix)
  ELSEIF( Key==2 ) THEN
    !
    !        CODE TO DO A PAGE WRITE.
    !
    CALL PRWVIR(Key,Ipage,Lpg,Sx,Ix)
  ELSE
    CALL XERMSG('PRWPGE',&
      'THE VALUE OF KEY (READ-WRITE FLAG) WAS NOT 1 OR 2.',55,1)
  END IF
END SUBROUTINE PRWPGE
